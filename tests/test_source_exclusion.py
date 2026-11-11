import hashlib
import json
from pathlib import Path
import random
import runpy
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "validate_source_exclusion.py"
V = runpy.run_path(str(SCRIPT))
G = V["Graph"]


def receipt(path):
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


class SourceExclusion(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.graph = self.root / "frozen.gfa"
        self.proof = self.root / "read-proof.json"
        self.proof.write_text('{"scope":"synthetic fixed graph and independent source truth"}\n')
        rng = random.Random(70327)
        self.bases = {n: "".join(rng.choices("ACGT", k=35))
                      for n in ("NL", "NR", "TL", "TR", "X")}
        self.prefix = ["H\tVN:Z:1.0", "# metadata untouched"]
        self.prefix += [f"S\t{n}\t{s}\tdp:f:50" for n, s in self.bases.items()]
        self.links = [
            "L\tNL\t+\tX\t+\t0M",
            "L\tX\t-\tNL\t-\t0M",
            "L\tX\t+\tNR\t+\t0M",
            "L\tNR\t-\tX\t-\t0M",
            "L\tTL\t+\tX\t+\t0M",
            "L\tX\t+\tTR\t+\t0M",
        ]
        self.target = self.bases["TL"] + self.bases["X"] + self.bases["TR"]
        self.competing = self.bases["NL"] + self.bases["X"] + self.bases["NR"]
        self.write_graph()

    def write_graph(self, extra=()):
        self.graph.write_bytes("\r\n".join(self.prefix + self.links + list(extra)).encode())

    def manifest(self, *, competitor_state="certified", protected=None):
        graph = receipt(self.graph)
        windows = []
        for name, role, sequence, state in [
            ("nuclear", "competing", self.competing, competitor_state),
            ("target", "protected", protected or self.target, "certified"),
        ]:
            cp = self.root / f"{name}.cert.json"
            cp.write_text(json.dumps({
                "schema": V["CERTIFICATE"], "graph_sha256": graph["sha256"],
                "window_id": name, "role": role, "window_sha256": V["sha"](sequence.encode()),
                "full_window": True, "status": state, "provenance": [receipt(self.proof)],
            }))
            windows.append({"id": name, "role": role, "sequence": sequence, "certificate": receipt(cp)})
        path = self.root / "manifest.json"
        path.write_text(json.dumps({
            "schema": V["MANIFEST"], "scope": "synthetic", "graph": graph,
            "windows": windows, "protected_scope_state": "complete_for_declared_scope",
        }))
        return path

    def invoke(self, manifest, report, view=None, *extra):
        argv = [sys.executable, "-B", "-I", str(SCRIPT), str(manifest), "--out", str(report)]
        if view is not None:
            argv += ["--out-gfa", str(view)]
        return subprocess.run(argv + list(extra), cwd=SCRIPT.parents[1], capture_output=True, text=True)

    def test_materialized_view_preserves_every_unaffected_byte(self):
        manifest = self.manifest()
        report = self.root / "audit.json"
        view = self.root / "view.gfa"
        before = self.graph.read_bytes()
        completed = self.invoke(manifest, report, view)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        result = json.loads(report.read_text())
        self.assertEqual(result["state"], "eligible_scoped_graph_exclusion")
        self.assertEqual(result["safe_edges"], [["NL+", "X+", 0], ["NR-", "X-", 0]])
        self.assertEqual(result["graph_view"]["removed_link_records"], 4)
        self.assertEqual(result["graph_view"]["sha256"], receipt(view)["sha256"])
        self.assertEqual(view.read_bytes(), "\r\n".join(self.prefix + self.links[4:]).encode())
        self.assertEqual(self.graph.read_bytes(), before)

    def test_unresolved_duplicate_and_budget_outputs_have_no_view(self):
        for name in ("unresolved", "alias", "budget"):
            with self.subTest(name=name):
                self.write_graph([f"S\tALT\t{self.competing}"] if name == "alias" else [])
                manifest = self.manifest(competitor_state="unresolved" if name == "unresolved" else "certified")
                audit = self.root / f"{name}.json"
                view = self.root / f"{name}.gfa"
                extra = ("--state-cap", "1") if name == "budget" else ()
                completed = self.invoke(manifest, audit, view, *extra)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                result = json.loads(audit.read_text())
                self.assertEqual(result["state"], "abstained")
                self.assertEqual(result["safe_edges"], [])
                self.assertNotIn("graph_view", result)
                self.assertFalse(view.exists())

    def test_annotated_paths_cannot_be_mistaken_for_a_materialized_view(self):
        self.write_graph(["P\tnuclear\tNL+,X+,NR+\t0M,0M"])
        manifest = self.manifest()
        verified = V["validate"](manifest)
        self.assertEqual(verified["state"], "eligible_scoped_graph_exclusion")
        report = self.root / "audit.json"
        view = self.root / "view.gfa"
        completed = self.invoke(manifest, report, view)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        result = json.loads(report.read_text())
        self.assertEqual(result["state"], "abstained")
        self.assertIn("P or W", result["reason"])
        self.assertFalse(view.exists())

    def test_graph_tags_changed_after_audit_prevent_materialization(self):
        result = V["validate"](self.manifest())
        self.graph.write_bytes(self.graph.read_bytes().replace(b"dp:f:50", b"dp:f:51", 1))
        with self.assertRaisesRegex(V["Abstain"], "SHA-256 mismatch"):
            V["materialize_view"](result)

    def test_existing_and_aliased_output_paths_are_protected(self):
        manifest = self.manifest()
        original = self.graph.read_bytes()
        audit = self.root / "report.json"
        completed = self.invoke(manifest, audit, self.graph)
        self.assertNotEqual(completed.returncode, 0)
        self.assertFalse(audit.exists())
        self.assertEqual(self.graph.read_bytes(), original)
        completed = self.invoke(manifest, audit, audit)
        self.assertNotEqual(completed.returncode, 0)
        self.assertFalse(audit.exists())

    def test_second_output_failure_rolls_back_only_the_new_view(self):
        existing = self.root / "existing.json"
        existing.write_bytes(b"keep")
        view = self.root / "candidate.gfa"
        with self.assertRaises(FileExistsError):
            V["publish_outputs"]([(view, b"view"), (existing, b"audit")])
        self.assertFalse(view.exists())
        self.assertEqual(existing.read_bytes(), b"keep")
        self.assertEqual(list(self.root.glob(".*")), [])


def brute(graph, query):
    """Independent full-offset, whole-spelling enumeration for tiny graphs."""
    found = set()
    outgoing = graph.outgoing()
    for first, sequence in graph.oriented.items():
        for start in range(len(sequence)):
            pending = [((first,), (), sequence)]
            while pending:
                path, overlaps, spelled = pending.pop()
                available = spelled[start:]
                if len(available) >= len(query):
                    if available[:len(query)] == query:
                        stop = len(graph.oriented[path[-1]]) - len(available) + len(query)
                        found.add(V["route_signature"](graph, path, overlaps, start, stop))
                elif query.startswith(available):
                    for following, overlap in outgoing.get(path[-1], ()):
                        pending.append((path + (following,), overlaps + (overlap,),
                                        spelled + graph.oriented[following][overlap:]))
    return found


class ExactWalks(unittest.TestCase):
    def test_both_orientations_match_the_independent_offset_oracle(self):
        rng = random.Random(48428)
        comparisons = 0
        for i in range(36):
            nodes = {n: "".join(rng.choices("ACGT", k=3 + i % 3)) for n in ("A", "B", "C")}
            ports = {k + sign: seq if sign == "+" else V["rc"](seq)
                     for k, seq in nodes.items() for sign in "+-"}
            links = {}
            for _ in range(4):
                a, b = rng.choices(list(ports), k=2)
                overlap = int(ports[a][-1:] == ports[b][:1] and rng.randrange(2) == 0)
                key = V["edge_key"](a, b, overlap)
                links.setdefault(key[:2], key)
            graph = G(ports, set(links.values()))
            queries = ["".join(rng.choices("ACGT", k=j)) for j in range(1, 8)]
            queries.extend(graph.oriented.values())
            queries.extend(graph.oriented[a] + graph.oriented[b][o:] for a, b, o in graph.edges)
            for query in queries:
                obtained = V["enumerate_routes"](graph, query, 200_000, 4096)
                self.assertTrue(obtained["complete"])
                projected = {
                    V["route_signature"](graph, r["walk"], r["overlaps"], r["entry_offset0"], r["exit_offset0"])
                    for r in obtained["routes"]
                }
                self.assertEqual(projected, brute(graph, query) | brute(graph, V["rc"](query)))
                comparisons += 1
        self.assertEqual(comparisons, 603)


if __name__ == "__main__":
    unittest.main()
