"""Focused controls for the long-read geometry CLI and exact anchor scan."""

from collections import Counter, defaultdict
from hashlib import sha256
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
             query_mismatch=False, empty=False, graph_contains_original=True):
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
