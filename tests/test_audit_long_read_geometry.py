"""Focused controls for the long-read geometry CLI and exact anchor scan."""

from collections import Counter, defaultdict
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from audit_long_read_geometry import canonical, exact_word_counts, packed_words, reverse
from discover_graph_frontier import DEFAULT_TRAVERSAL_STATES, evaluate

SCRIPT = SCRIPTS / "audit_long_read_geometry.py"


def random_dna(seed, size):
    rng = random.Random(seed)
    return bytes(rng.choice(b"ACGT") for _ in range(size))


def fasta(name, sequence):
    return b">" + name.encode() + b"\n" + sequence + b"\n"


def fastq(name, movie, sequence):
    return b"@" + name.encode() + b" " + movie.encode() + b"/ccs\n" + sequence + b"\n+\n" + b"I" * len(sequence) + b"\n"


class ExactScanTests(unittest.TestCase):
    def test_overlap_reverse_palindrome_unknown_bases_and_widths(self):
        seq = b"AAAATATATNTTTANCTAG" + random_dna(421, 209)
        queries = [b"AAAA", b"AAA", b"ATAT", b"TAG", b"CTAG", b"GT", b"T", seq[60:75]]
        words = defaultdict(set)
        for word in queries:
            words[len(word)].add(canonical(word))
        expected = Counter()
        for width, keys in words.items():
            for start in range(len(seq) - width + 1):
                word = seq[start:start + width]
                if canonical(word) in keys:
                    expected[canonical(word)] += 1
        self.assertEqual(exact_word_counts(seq, packed_words(words)), expected)
        self.assertEqual(exact_word_counts(b"ATATAT", packed_words({4: {b"ATAT"}})), {b"ATAT": 2})
        self.assertEqual(exact_word_counts(b"AAANAAA", packed_words({4: {b"AAAA"}})), {})


class GeometryCliTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.mt = random_dna(100, 128)
        self.nuclear = random_dna(200, 130)
        self.plastid = random_dna(300, 102)
        self.query = reverse(self.mt[10:110])
        self.original = fastq("original", "m123456/17", self.query)

    def file(self, relative_path, content):
        (self.root / relative_path).write_bytes(content)
        return dict(path=relative_path, sha256=sha256(content).hexdigest())

    def plan(self, *, shared_guard=False, same_scope_competitor=False, donor=True,
             query_mismatch=False, empty=False, graph_contains_original=True, extra_paf="", skip=False):
        nucleus = self.nuclear
        if shared_guard:
            # Exact obsolete right anchor on the reverse strand, absent from the PAF.
            nucleus += reverse(self.mt[74:85])
        if same_scope_competitor:
            close = bytearray(self.mt[10:110])
            for i in range(4, len(close), 7):
                close[i] = b"ACGT"[(b"ACGT".index(close[i]) + 1) % 4]
            nucleus += bytes(close)
        records = {"mt": self.mt, "plastid": self.plastid, "chr": nucleus, "chr2": random_dna(432, 104)}
        mt_component = self.file("mt.fa", fasta("mt_orig", self.mt))
        pt_component = self.file("pt.fa", fasta("pt_orig", self.plastid))
        nuclear_component = self.file("nuc.fa", fasta("c1", nucleus) + fasta("c2", records["chr2"]))
        combined = self.file("panel.fa", b"".join(fasta(k, v) for k, v in records.items()))
        panel = dict(fasta=combined, records={
            "mt": dict(role="mitochondrion", source=mt_component, source_record="mt_orig"),
            "plastid": dict(role="plastid", source=pt_component, source_record="pt_orig"),
            "chr": dict(role="nuclear", source=nuclear_component, source_record="c1"),
            "chr2": dict(role="nuclear", source=nuclear_component, source_record="c2"),
        })
        matches = sum(a == b for a, b in zip(reverse(self.query), self.mt[10:110]))
        paf = f"original\t100\t0\t100\t-\tmt\t128\t10\t110\t{matches}\t100\t60\tcg:Z:100M\n"
        if skip:
            paf = f"original\t100\t0\t100\t-\tmt\t128\t9\t110\t{matches}\t101\t60\tcg:Z:1N100M\n"
        paf += extra_paf
        if same_scope_competitor:
            paf += (f"original\t100\t0\t100\t-\tchr\t{len(nucleus)}\t{len(self.nuclear)}\t"
                    f"{len(nucleus)}\t86\t100\t2\tcg:Z:100M\n")
        query = self.query if not query_mismatch else b"A" + self.query
        p = dict(schema="ovasm.long-read-geometry.v1", scope="Synthetic capped mapper evidence",
                 original_fastqs=[self.file("reads.fastq", b"" if empty else self.original)],
                 mapper_queries=self.file("queries.fa", b"" if empty else fasta("original", query)),
                 paf=dict(data=self.file("reads.paf", b"" if empty else paf.encode("ascii")),
                          provenance=self.file("paf-proof.txt", b"synthetic mapper with a secondary cap\n")),
                 panel=panel,
                 windows=[dict(id="window", reference="mt", interval=[50, 74], anchor_bp=11,
                               mode="exact_context", exclusive_sides=["left", "right"],
                               graph_required=True, physical_boundaries=[], require_full_query=True)],
                 graph=dict(gfa=self.file("graph.gfa", b"S\tnode\t" + self.mt + b"\n"),
                            reads=[self.file("graph.fastq", (self.original if graph_contains_original else b"") +
                                             (fastq("donor", "m123456/22", self.query) if donor else b""))],
                            provenance=self.file("graph-proof.txt", b"synthetic graph input relation\n")),
                 window_provenance=self.file("window-proof.txt", b"registered independent window location\n"))
        (self.root / "plan.json").write_text(json.dumps(p))

    def run_cli(self, expected_success=True):
        out = self.root / "geometry.json"
        before = set(self.root.iterdir())
        p = subprocess.run([sys.executable, "-B", str(SCRIPT), str(self.root / "plan.json"),
                            "--out", str(out)], capture_output=True, text=True, check=False)
        self.assertEqual(p.returncode == 0, expected_success, p.stderr)
        if expected_success:
            self.assertEqual(set(self.root.iterdir()) - before, {out})
            result = json.loads(out.read_bytes())
            self.assertFalse(result["native_admission_created"])
            self.assertFalse(result["competition_complete"])
            return result
        self.assertEqual(set(self.root.iterdir()), before)
        return p.stderr

    def test_capped_paf_keeps_only_geometry_and_an_independent_graph_donor(self):
        self.plan()
        r = self.run_cli()
        w = r["records"][0]["windows"]["window"]
        self.assertEqual(w["competition_state"], "unresolved_reported_paf_only")
        self.assertEqual(w["other_graph_input_exact_donors"], ["m123456/22"])
        self.assertTrue(w["graph_corroborated_observation"])
        self.assertEqual(w["witnesses"][0]["placement"]["strand"], "-")
        self.assertEqual(w["witnesses"][0]["right"]["reference"], [74, 85])
        self.assertEqual(r["panel_records_by_role"]["nuclear"]["records"], 2)

    def test_hidden_reverse_complement_nuclear_match_moves_the_anchor(self):
        self.plan(shared_guard=True)
        r = self.run_cli()
        w = r["records"][0]["windows"]["window"]
        self.assertEqual(len(w["witnesses"]), 1)
        self.assertGreater(w["witnesses"][0]["right"]["reference"][0], 74)
        self.assertEqual(w["witnesses"][0]["reported_paf_roles_covering_same_scope"], ["mitochondrion"])

    def test_reported_same_scope_competitor_removes_unopposed_graph_observation(self):
        self.plan(same_scope_competitor=True)
        r = self.run_cli()
        w = r["records"][0]["windows"]["window"]
        self.assertEqual(w["witnesses"][0]["reported_paf_roles_covering_same_scope"], ["mitochondrion", "nuclear"])
        self.assertFalse(w["graph_corroborated_observation"])

    def test_candidate_does_not_supply_its_own_independent_graph_donor(self):
        self.plan(donor=False)
        r = self.run_cli()
        w = r["records"][0]["windows"]["window"]
        self.assertTrue(w["graph_spelling_unique"])
        self.assertEqual(w["other_graph_input_exact_donors"], [])
        self.assertFalse(w["graph_corroborated_observation"])

    def test_query_sequence_mismatch_fails_closed(self):
        self.plan(query_mismatch=True)
        self.assertIn("mapper query bytes or enrollment differ", self.run_cli(False))

    def test_panel_byte_mismatch_fails_closed(self):
        self.plan()
        self.file("panel.fa", fasta("mt", self.mt))
        self.assertIn("input SHA-256 mismatch", self.run_cli(False))

    def test_empty_enrollment_is_a_valid_unresolved_report(self):
        self.plan(empty=True)
        r = self.run_cli()
        self.assertEqual(r["original_molecules_enrolled"], 0)
        self.assertEqual(r["records"], [])


    def test_full_edit_export_binds_exact_report_and_original_record(self):
        self.plan(graph_contains_original=False)
        report, manifest = self.run_full_edit()
        self.assertTrue(report["native_admission_created"])
        self.assertFalse(report["competition_complete"])
        row = report["records"][0]
        model = row["full_edit_competition"]
        self.assertEqual(model["status"], "qualified_within_edit_model")
        self.assertEqual(model["radius"], 0)
        self.assertEqual(set(model["competing_roles"]), {"nuclear", "plastid"})
        self.assertTrue(all(v["status"] == "complete" and v["best_hit"] is None
                            for v in model["competing_roles"].values()))
        self.assertEqual(manifest["schema"], "ovasm.read-admission.v1")
        record, = manifest["records"]
        for key in ("read_id_hex", "sequence_sha256", "quality_sha256"):
            self.assertEqual(record[key], row[key])
        out = self.root / "full.geometry.json"
        self.assertEqual(record["provenance"], {"path": str(out), "sha256": sha256(out.read_bytes()).hexdigest()})

    def adjacent_boundary_plan(self, query, cigar, *, mode="adjacent_boundary", mismatches=0):
        self.plan(graph_contains_original=False)
        path = self.root / "plan.json"
        plan = json.loads(path.read_bytes())
        raw = reverse(query)
        length = len(raw)
        plan["windows"][0]["mode"] = mode
        plan["original_fastqs"] = [self.file("reads.fastq", fastq("original", "m123456/17", raw))]
        plan["mapper_queries"] = self.file("queries.fa", fasta("original", raw))
        paf = (f"original\t{length}\t0\t{length}\t-\tmt\t128\t10\t110\t"
               f"{min(100, length) - mismatches}\t{max(100, length)}\t60\tcg:Z:{cigar}\n")
        plan["paf"]["data"] = self.file("reads.paf", paf.encode("ascii"))
        path.write_text(json.dumps(plan))

    def test_interior_indel_with_two_exact_junctions_qualifies(self):
        reference = self.mt[10:110]
        self.adjacent_boundary_plan(reference[:48] + b"A" + reference[48:], "48M1I52M")
        report, manifest = self.run_full_edit()
        row = report["records"][0]
        model = row["full_edit_competition"]
        self.assertEqual(model["status"], "qualified_within_edit_model")
        self.assertEqual(model["radius"], 1)
        self.assertEqual(model["verified_target_upper_bound"]["cigar"], "48M1I52M")
        self.assertFalse(model["verified_target_upper_bound"]["witness"]["exact_boundary"])
        self.assertEqual(row["windows"]["window"]["other_graph_input_exact_donors"], ["m123456/22"])
        self.assertEqual(manifest["records"][0]["sequence_sha256"], row["sequence_sha256"])

    def test_indel_at_junction_remains_geometry_without_admission(self):
        reference = self.mt[10:110]
        self.adjacent_boundary_plan(reference[:40] + reference[41:], "40M1D59M")
        report, manifest = self.run_full_edit()
        self.assertIsNone(manifest)
        row = report["records"][0]
        window = row["windows"]["window"]
        self.assertTrue(window["witnesses"])
        self.assertEqual(window["excluded_alignment_counts"], {"exact_junction_missing": 1})
        self.assertEqual(row["full_edit_competition"]["status"], "unresolved_no_full_target_geometry")

    def test_junction_insertion_blocks_admission_even_if_core_is_exact(self):
        reference = self.mt[10:110]
        self.adjacent_boundary_plan(reference[:40] + b"A" + reference[40:], "40M1I60M")
        report, manifest = self.run_full_edit()
        self.assertIsNone(manifest)
        row = report["records"][0]
        window = row["windows"]["window"]
        self.assertTrue(window["witnesses"][0]["exact_boundary"])
        self.assertEqual(window["excluded_alignment_counts"], {"exact_junction_missing": 1})
        self.assertEqual(row["full_edit_competition"]["status"], "unresolved_no_full_target_geometry")

    def test_isolated_graph_core_does_not_corroborate_adjacent_junctions(self):
        reference = self.mt[10:110]
        self.adjacent_boundary_plan(reference[:48] + b"A" + reference[48:], "48M1I52M")
        path = self.root / "plan.json"
        plan = json.loads(path.read_bytes())
        plan["graph"]["gfa"] = self.file("graph.gfa", b"S\tnode\t" + self.mt[50:74] + b"\n")
        plan["graph"]["reads"] = [self.file("graph.fastq", fastq("donor", "m123456/22", self.mt[50:74]))]
        path.write_text(json.dumps(plan))
        report, manifest = self.run_full_edit()
        self.assertIsNone(manifest)
        row = report["records"][0]
        window = row["windows"]["window"]
        self.assertTrue(window["witnesses"])
        self.assertFalse(window["graph_spelling_unique"])
        self.assertEqual(window["other_graph_input_exact_donors"], [])
        self.assertEqual(row["full_edit_competition"]["status"], "unresolved_no_independent_graph_support")

    def test_interior_indel_with_hidden_equal_competitor_is_withheld(self):
        reference = self.mutate(self.mt[10:110], [0])
        oriented = reference[:48] + b"A" + reference[48:]
        self.nuclear += self.mutate(oriented, [36, 68])
        self.adjacent_boundary_plan(oriented, "48M1I52M", mismatches=1)
        report, manifest = self.run_full_edit()
        self.assertIsNone(manifest)
        row = report["records"][0]
        self.assertTrue(row["windows"]["window"]["graph_corroborated_observation"])
        model = row["full_edit_competition"]
        self.assertEqual(model["radius"], 2)
        self.assertEqual(model["status"], "competing_within_edit_bound")
        self.assertEqual(model["competing_roles"]["nuclear"]["best_hit"]["edit_distance"], 2)

    def test_exact_context_still_requires_an_exact_candidate_core(self):
        reference = self.mt[10:110]
        self.adjacent_boundary_plan(reference[:48] + b"A" + reference[48:], "48M1I52M", mode="exact_context")
        report, manifest = self.run_full_edit()
        self.assertIsNone(manifest)
        row = report["records"][0]
        self.assertEqual(row["windows"]["window"]["excluded_alignment_counts"], {"exact_core_missing": 1})
        self.assertEqual(row["full_edit_competition"]["status"], "unresolved_no_full_target_geometry")

    def test_selected_full_alignment_is_auditable_even_if_geometry_deduplicates(self):
        alternate = "original\t100\t0\t100\t-\tmt\t128\t10\t110\t99\t101\t25\tcg:Z:1I1D99M\n"
        self.plan(graph_contains_original=False, extra_paf=alternate)
        report, manifest = self.run_full_edit()
        self.assertIsNotNone(manifest)
        row = report["records"][0]
        self.assertEqual(len(row["windows"]["window"]["witnesses"]), 1)
        upper = row["full_edit_competition"]["verified_target_upper_bound"]
        self.assertEqual(upper["unit_edit_count"], 0)
        self.assertEqual(upper["cigar"], "100M")
        self.assertEqual(upper["witness"]["scope"]["read_sha256"], row["sequence_sha256"])
        self.assertEqual(upper["witness"]["boundary"]["reference"], [50, 74])

    def test_terminal_skip_keeps_geometry_without_full_edit_qualification(self):
        self.plan(graph_contains_original=False, skip=True)
        report, manifest = self.run_full_edit()
        self.assertIsNone(manifest)
        row = report["records"][0]
        self.assertTrue(row["windows"]["window"]["witnesses"])
        self.assertEqual(row["full_edit_competition"]["status"], "unresolved_no_full_target_geometry")

    def test_unreported_continuous_edit_competitor_prevents_export(self):
        oriented = self.mutate(self.mt[10:110], [0, 1, 2, 3])
        self.query = reverse(oriented)
        self.original = fastq("original", "m123456/17", self.query)
        self.nuclear += self.mutate(oriented, [36, 68])
        self.plan(graph_contains_original=False)
        report, manifest = self.run_full_edit()
        self.assertIsNone(manifest)
        self.assertFalse(report["native_admission_created"])
        row = report["records"][0]
        window = row["windows"]["window"]
        self.assertTrue(window["graph_corroborated_observation"])
        self.assertEqual(window["witnesses"][0]["reported_paf_roles_covering_same_scope"], ["mitochondrion"])
        model = row["full_edit_competition"]
        self.assertEqual(model["radius"], 4)
        self.assertEqual(model["status"], "competing_within_edit_bound")
        self.assertEqual(model["competing_roles"]["nuclear"]["best_hit"]["edit_distance"], 2)

    def frontier_plan(self, *, explicit=True, negative=False, core_only_donor=False,
                      bind_gapped=False, duplicate_geometry=False, nuclear_diagnostic=None):
        left, right = random_dna(17, 60), random_dna(23, 60)
        ordinary = left + right
        target = left + random_dna(34, 21) + right
        first = reverse(target) if negative else target
        strand = "-" if negative else "+"
        originals = self.file("originals.fastq", fastq("A", "m123456/17", first) +
                              fastq("B", "m123456/19", ordinary))
        queries = self.file("queries.fa", fasta("A", first) + fasta("B", ordinary))
        gfa = self.file("graph.gfa", b"S\tcore\t" + ordinary + b"\n")
        discovery = self.root / "discovery-plan.json"
        discovery.write_text(json.dumps({"schema": "ovasm.graph-frontier-discovery.v1",
                                         "scope": "Synthetic original differences", "graph": gfa,
                                         "raw_reads": [originals], "scan_k": 9}))
        evidence = BytesIO()
        report = evaluate(discovery, evidence, DEFAULT_TRAVERSAL_STATES, False)
        row, = (json.loads(line) for line in evidence.getvalue().splitlines())
        self.assertEqual(row["record_index0"], 0)
        event, = row["events"]
        start = event["affected_raw_interval0"][0]
        raw_core = start - 9, start
        core = (len(first) - raw_core[1], len(first) - raw_core[0]) if negative else raw_core
        evidence_file = self.file("discovery.jsonl", evidence.getvalue())
        report["evidence"] = evidence_file | {"path": str(self.root / evidence_file["path"])}
        report_file = self.file("discovery.json", json.dumps(report).encode())
        donor = target[slice(*core)] if core_only_donor else ordinary
        fixed = f"A\t141\t0\t141\t{strand}\tmt\t141\t0\t141\t"
        matches = 141 if nuclear_diagnostic is None else 139
        exact = (fixed + f"{matches}\t141\t60\tcg:Z:141M\n").encode()
        gapped = (fixed + "140\t142\t20\tcg:Z:1I1D140M\n").encode()
        alias = gapped.replace(b"\t20\tcg:", b"\t23\tcg:")
        second = b"B\t120\t0\t120\t+\tmt\t141\t0\t141\t120\t141\t60\tcg:Z:60M21D60M\n"
        diagnostic_row = b"A\t141\t0\t141\t+\tchr\t141\t0\t141\t138\t141\t60\tcg:Z:141M\n"
        lines = exact + gapped + (alias if duplicate_geometry else b"") + second
        if nuclear_diagnostic:
            lines += diagnostic_row
        p = dict(schema="ovasm.long-read-geometry.v1", scope="One real graph event and one matching graph-spelled read",
                 original_fastqs=[originals], mapper_queries=queries,
                 windows=[dict(id="window", reference="mt", interval=core, anchor_bp=9, mode="exact_context",
                               exclusive_sides=["left", "right"], graph_required=True,
                               physical_boundaries=[], require_full_query=True)],
                 paf=dict(data=self.file("reads.paf", lines), provenance=self.file("paf-proof.txt", b"fixture\n")),
                 graph=dict(gfa=gfa, reads=[self.file("graph.fastq", fastq("donor", "m123456/22", donor))],
                            provenance=self.file("graph-proof.txt", b"fixture\n")),
                 window_provenance=self.file("window-proof.txt", b"derived from a frozen event\n"))
        mitochondrial, nuclear = target, random_dna(27, 250)
        if nuclear_diagnostic is not None:
            mitochondrial = self.mutate(target, [15, 122])
            nuclear = self.mutate(target, [35, 96, 108])
            p["windows"].append(dict(p["windows"][0], id="nuclear_diagnostic", reference="chr",
                                     interval=[70, 79], graph_required=False))
        references = {"mt": ("mitochondrion", mitochondrial), "pt": ("plastid", random_dna(19, 220)),
                      "chr": ("nuclear", nuclear)}
        p["panel"] = dict(fasta=self.file("panel.fa", b"".join(fasta(n, v) for n, (_, v) in references.items())),
                          records={n: dict(role=role, source=self.file(n + ".fa", fasta(n, seq)), source_record=n)
                                   for n, (role, seq) in references.items()})
        fields = ("read_id_hex", "sequence_sha256", "quality_sha256")
        raw_identity = {k: row[k] for k in fields} | {"molecule": "m123456/17"}
        event_hash = sha256(json.dumps(raw_identity | {"event_index0": 0, "event": event},
                                      sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        selected_paf = alias if bind_gapped and duplicate_geometry else gapped if bind_gapped else exact
        registration = raw_identity | {"input_index0": 0, "record_index0": 0, "event_index0": 0,
                                       "event_sha256": event_hash, "window_id": "window", "side": "left",
                                       "paf_row_sha256": sha256(selected_paf).hexdigest()}
        binding = {"schema": "ovasm.frontier-event-bindings.v1", "discovery_report": report_file,
                   "discovery_evidence": evidence_file, "records": [registration]}
        if explicit:
            p["frontier_events"] = self.file("frontier.json", json.dumps(binding).encode())
        (self.root / "plan.json").write_text(json.dumps(p))
        return p, binding

    def test_window_only_can_admit_an_original_without_a_frontier_event(self):
        self.frontier_plan(explicit=False)
        report, manifest = self.run_full_edit()
        self.assertEqual(report["full_edit_competition"]["eligible_records"], 2)
        self.assertEqual([bytes.fromhex(r["read_id_hex"]).split()[0] for r in manifest["records"]], [b"A", b"B"])

    def test_frontier_binding_excludes_graph_spelled_original_and_preserves_reverse_strand(self):
        self.frontier_plan(negative=True)
        report, manifest = self.run_full_edit()
        self.assertEqual(report["original_molecules_scanned"], 2)
        self.assertEqual(report["original_molecules_enrolled"], 1)
        self.assertEqual(report["frontier_event_binding"]["ignored_paf_rows"], 1)
        row, = report["records"]
        self.assertEqual(row["query"], "A")
        upper = row["full_edit_competition"]["verified_target_upper_bound"]
        self.assertEqual(upper["strand"], "-")
        self.assertEqual(upper["frontier_events"][0]["independent_donor_molecules"], ["m123456/22"])
        self.assertEqual(manifest["records"][0]["read_id_hex"], row["read_id_hex"])

    def test_only_bound_paf_can_supply_the_edit_bound_even_with_identical_geometry_rows(self):
        _, binding = self.frontier_plan(bind_gapped=True, duplicate_geometry=True)
        report, manifest = self.run_full_edit()
        self.assertIsNotNone(manifest)
        upper = report["records"][0]["full_edit_competition"]["verified_target_upper_bound"]
        self.assertEqual(upper["unit_edit_count"], 2)
        self.assertEqual(upper["cigar"], "1I1D140M")
        event, = upper["frontier_events"]
        self.assertEqual(event["paf_row_sha256"], binding["records"][0]["paf_row_sha256"])

    def test_core_only_donor_cannot_support_a_frontier_event(self):
        self.frontier_plan(core_only_donor=True)
        report, manifest = self.run_full_edit()
        self.assertIsNone(manifest)
        row = report["records"][0]
        self.assertTrue(row["windows"]["window"]["other_graph_input_exact_donors"])
        self.assertEqual(row["full_edit_competition"]["status"], "unresolved_no_verified_frontier_target_geometry")
        self.assertEqual(report["frontier_event_binding"]["verified_without_independent_donor"], 1)

    def test_unregistered_other_role_diagnostic_can_veto_the_registered_target(self):
        self.frontier_plan(nuclear_diagnostic=True)
        report, manifest = self.run_full_edit()
        self.assertIsNone(manifest)
        row = report["records"][0]
        self.assertTrue(row["windows"]["nuclear_diagnostic"]["witnesses"])
        self.assertEqual(row["full_edit_competition"]["status"], "unresolved_conflicting_geometry")
        self.assertEqual(report["frontier_event_binding"]["verified_with_independent_donor"], 1)

    def test_unobserved_diagnostic_preserves_the_bound_target_search(self):
        self.frontier_plan(nuclear_diagnostic=False)
        report, manifest = self.run_full_edit()
        self.assertIsNotNone(manifest)
        row = report["records"][0]
        self.assertEqual(row["windows"]["nuclear_diagnostic"]["witnesses"], [])
        model = row["full_edit_competition"]
        self.assertEqual(model["status"], "qualified_within_edit_model")
        self.assertEqual(model["radius"], 2)
        self.assertIsNone(model["competing_roles"]["nuclear"]["best_hit"])

    def test_wrong_event_or_paf_binding_publishes_no_artifacts(self):
        plan, binding = self.frontier_plan()
        original = binding["records"][0]
        for key in ("event_sha256", "paf_row_sha256", "quality_sha256", "record_index0"):
            with self.subTest(key=key):
                binding["records"] = [original | {key: 1 if key == "record_index0" else "0" * 64}]
                plan["frontier_events"] = self.file("frontier.json", json.dumps(binding).encode())
                plan_path = self.root / "plan.json"
                plan_path.write_text(json.dumps(plan))
                out = self.root / "audit.json"
                admission = self.root / "admission.json"
                args = [sys.executable, "-B", str(SCRIPT), str(plan_path), "--out", str(out),
                        "--full-edit-competition", "--out-admission", str(admission)]
                result = subprocess.run(args, capture_output=True, text=True, check=False)
                self.assertEqual(result.returncode, 1, result.stderr)
                self.assertIn("frontier", result.stderr)
                self.assertFalse(out.exists() or admission.exists())

    def competing_source_plan(self, *, partial=False, close_mitochondrion=False):
        self.plan(graph_contains_original=False)
        path = self.root / "plan.json"
        plan = json.loads(path.read_bytes())
        oriented = self.mutate(self.nuclear[10:110], [5, 90])
        query = reverse(oriented)
        plan["original_fastqs"] = [self.file("reads.fastq", fastq("original", "m123456/17", query))]
        plan["mapper_queries"] = self.file("queries.fa", fasta("original", query))
        query_start, ref_end = (1, 109) if partial else (0, 110)
        aligned = ref_end - 10
        paf = (f"original\t100\t{query_start}\t100\t-\tchr\t130\t10\t{ref_end}\t"
               f"{aligned - 2}\t{aligned}\t60\tcg:Z:{aligned}M\n")
        plan["paf"]["data"] = self.file("reads.paf", paf.encode("ascii"))
        plan["windows"][0].update(reference="chr", graph_required=False, require_full_query=False)
        if close_mitochondrion:
            mt = self.mt + self.mutate(oriented, [30, 75])
            plan["panel"]["records"]["mt"]["source"] = self.file("mt.fa", fasta("mt_orig", mt))
            with (self.root / "panel.fa").open("rb") as original_panel:
                other_references = original_panel.read().split(b"\n", 2)[2]
            plan["panel"]["fasta"] = self.file("panel.fa", fasta("mt", mt) + other_references)
        path.write_text(json.dumps(plan))

    def run_source_exclusions(self, *args):
        out = self.root / "full.geometry.json"
        exclusion = self.root / "source-exclusions.json"
        command = [sys.executable, "-B", str(SCRIPT), str(self.root / "plan.json"),
                   "--out", str(out), "--full-edit-competition", "--exclude-target", "mitochondrion",
                   "--out-exclusions", str(exclusion), *args]
        process = subprocess.run(command, capture_output=True, text=True, check=False)
        self.assertEqual(process.returncode, 0, process.stderr)
        return json.loads(out.read_bytes()), json.loads(exclusion.read_bytes()) if exclusion.exists() else None

    def test_full_competing_source_excludes_only_the_named_destination(self):
        self.competing_source_plan()
        report, manifest = self.run_source_exclusions()
        row, = report["records"]
        decision = row["full_edit_source_exclusion"]
        self.assertEqual(decision["status"], "qualified_within_edit_model")
        self.assertEqual(decision["source_role"], "nuclear")
        self.assertEqual(decision["radius"], 2)
        witness = decision["verified_source_upper_bound"]["witness"]
        self.assertEqual(witness["scope"]["interval"], [0, 100])
        self.assertEqual(witness["placement"]["strand"], "-")
        self.assertEqual(set(decision["competing_roles"]), {"mitochondrion", "plastid"})
        self.assertTrue(report["native_exclusions_created"])
        self.assertFalse(report["native_admission_created"])
        self.assertEqual(manifest["schema"], "ovasm.read-source-exclusions.v1")
        original, = manifest["records"]
        self.assertEqual(original["target"], "mitochondrion")
        for key in ("read_id_hex", "sequence_sha256", "quality_sha256"):
            self.assertEqual(original[key], row[key])
        out = self.root / "full.geometry.json"
        self.assertEqual(original["provenance"],
                         {"path": str(out), "sha256": sha256(out.read_bytes()).hexdigest()})

    def test_partial_competing_witness_remains_unresolved(self):
        self.competing_source_plan(partial=True)
        report, manifest = self.run_source_exclusions()
        row, = report["records"]
        self.assertTrue(row["windows"]["window"]["witnesses"])
        self.assertEqual(row["full_edit_source_exclusion"]["status"], "unresolved_no_full_competing_geometry")
        self.assertIsNone(manifest)

    def test_unreported_equal_destination_placement_withholds_exclusion(self):
        self.competing_source_plan(close_mitochondrion=True)
        report, manifest = self.run_source_exclusions()
        row, = report["records"]
        decision = row["full_edit_source_exclusion"]
        self.assertEqual(decision["status"], "competing_within_edit_bound")
        self.assertEqual(decision["competing_roles"]["mitochondrion"]["best_hit"]["edit_distance"], 2)
        self.assertIsNone(manifest)

    def test_registered_partial_other_role_withholds_exclusion(self):
        self.competing_source_plan()
        path = self.root / "plan.json"
        plan = json.loads(path.read_bytes())
        oriented = self.mutate(self.nuclear[10:110], [5, 14, 23, 90])
        query = reverse(oriented)
        plan["original_fastqs"] = [self.file("reads.fastq", fastq("original", "m123456/17", query))]
        plan["mapper_queries"] = self.file("queries.fa", fasta("original", query))
        mt = self.mt + oriented[:30]
        plan["panel"]["records"]["mt"]["source"] = self.file("mt.fa", fasta("mt_orig", mt))
        other_references = (self.root / "panel.fa").read_bytes().split(b"\n", 2)[2]
        plan["panel"]["fasta"] = self.file("panel.fa", fasta("mt", mt) + other_references)
        plan["windows"].append(dict(id="mt_diagnostic", reference="mt", interval=[140, 148], anchor_bp=8,
                                    mode="exact_context", exclusive_sides=["left", "right"],
                                    graph_required=True, physical_boundaries=[], require_full_query=True))
        paf = (b"original\t100\t0\t100\t-\tchr\t130\t10\t110\t96\t100\t60\tcg:Z:100M\n"
               b"original\t100\t70\t100\t-\tmt\t158\t128\t158\t30\t30\t60\tcg:Z:30M\n")
        plan["paf"]["data"] = self.file("reads.paf", paf)
        path.write_text(json.dumps(plan))
        report, manifest = self.run_source_exclusions()
        row, = report["records"]
        self.assertEqual(row["full_edit_source_exclusion"]["status"], "unresolved_conflicting_geometry")
        self.assertFalse(row["windows"]["mt_diagnostic"]["witnesses"])
        self.assertTrue(row["windows"]["mt_diagnostic"]["partial_witnesses_for_source_exclusion"])
        self.assertIsNone(manifest)

    def test_competition_budget_exhaustion_withholds_exclusion(self):
        self.competing_source_plan(close_mitochondrion=True)
        report, manifest = self.run_source_exclusions("--max-edit-text-bases", "1")
        row, = report["records"]
        self.assertEqual(row["full_edit_source_exclusion"]["status"], "unresolved_edit_budget")
        self.assertIsNone(manifest)

    def run_full_edit(self):
        out = self.root / "full.geometry.json"
        admission = self.root / "read-admission.json"
        p = subprocess.run([sys.executable, "-B", str(SCRIPT), str(self.root / "plan.json"),
                            "--out", str(out), "--full-edit-competition", "--out-admission", str(admission)],
                           capture_output=True, text=True, check=False)
        self.assertEqual(p.returncode, 0, p.stderr)
        return json.loads(out.read_bytes()), json.loads(admission.read_bytes()) if admission.exists() else None

    @staticmethod
    def mutate(sequence, offsets):
        mutated = bytearray(sequence)
        for pos in offsets:
            mutated[pos] = ord("C") if mutated[pos] == ord("A") else ord("A")
        return bytes(mutated)



if __name__ == "__main__":
    unittest.main()
