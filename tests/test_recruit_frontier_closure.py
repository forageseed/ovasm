import importlib.util
from itertools import product
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/recruit_frontier_closure.py"
SPEC = importlib.util.spec_from_file_location("ovasm_recruit_frontier_closure", SCRIPT)
TOOL = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = TOOL
SPEC.loader.exec_module(TOOL)


class FrontierClosureTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def fixture(self, label, edges, evidence, *, seeds=("a",), coverage="complete", duplicate=False):
        directory = self.root / label
        directory.mkdir()

        def save(filename, data):
            path = directory / filename
            path.write_text(json.dumps(data) + "\n")
            return {"path": filename, "sha256": TOOL.sha(path.read_bytes())}

        graph = save("graph.json", {"schema": TOOL.GRAPH, "nodes": sorted(set(seeds) | {
            n for _, source, destination in edges for n in (source, destination)
        }), "edges": [{"id": edge, "from": source, "to": destination}
                        for edge, source, destination in edges]})
        eligibility = save("eligibility.json", {"rule": "locked fixture policy"})
        read_pool = save("read-pool.json", {"read_pool": "immutable fixture fragments"})
        provenance = save("provenance.json", {"claim": "synthetic observations with known labels"})
        binding = {"scope": "fixture", "graph_sha256": graph["sha256"],
                   "eligibility_sha256": eligibility["sha256"], "read_pool_sha256": read_pool["sha256"],
                   "seed_set_sha256": TOOL.canonical_sha(sorted(set(seeds)))}
        scope = save("scope.json", {"schema": TOOL.SCOPE, **binding, "coverage": coverage,
                                   "provenance": [provenance], "explanation": "fixture declaration"})

        records = []
        for index, (edge, decision, basis) in enumerate(evidence):
            opp = [] if basis in ("no_opportunity", "unmeasured") else ["frag-b", "frag-a"]
            support = ["frag-a"] if basis in ("target_exact_bridge", "competing_context") else []
            if duplicate:
                opp = opp + opp
                support = support + support
            records.append(save(f"edge-{index}.json", {
                "schema": TOOL.EDGE, **binding, "edge_id": edge, "decision": decision, "basis": basis,
                "full_edge": decision != "unresolved", "opportunity_fragment_ids": opp,
                "support_fragment_ids": support, "provenance": [provenance],
                "explanation": "synthetic fixed opportunity",
            }))
        if duplicate:
            records += records
        manifest = save("manifest.json", {
            "schema": TOOL.MANIFEST, "scope": "fixture", "graph": graph, "eligibility": eligibility,
            "read_pool": read_pool, "scope_receipt": scope, "seeds": list(seeds), "edge_evidence": records,
        })
        return directory / manifest["path"]

    def test_order_and_duplicate_evidence_leave_the_same_closed_reachability(self):
        edges = [("ab", "a", "b"), ("bc", "b", "c"), ("ca", "c", "a"),
                 ("bd", "b", "d"), ("parallel", "a", "b"), ("hidden", "d", "e")]
        evidence = [(edge, "admitted", "target_exact_bridge") for edge in ("ab", "bc", "ca")]
        evidence += [("bd", "excluded", "competing_context"),
                     ("parallel", "excluded", "calibrated_absence")]
        first = TOOL.verify(self.fixture("one", edges, evidence, seeds=("isolated", "a")))
        second = TOOL.verify(self.fixture("two", edges[::-1], evidence[::-1],
                                          seeds=("a", "a", "isolated"), duplicate=True))
        self.assertNotEqual(first["manifest_sha256"], second["manifest_sha256"])
        self.assertEqual(first["decision_sha256"], second["decision_sha256"])
        self.assertEqual(second["seed_nodes"], ["a", "isolated"])
        self.assertEqual(second["reachable_nodes"], ["a", "b", "c", "isolated"])
        self.assertEqual(second["admitted_edges"], ["ab", "bc", "ca"])
        self.assertEqual(second["excluded_edges"], ["bd", "parallel"])
        self.assertEqual(second["state"], "converged")
        hidden = next(edge for edge in second["edges"] if edge["id"] == "hidden")
        self.assertEqual(hidden["state"], "not_reached")
        ab = next(edge for edge in second["edges"] if edge["id"] == "ab")
        record, = ab["evidence"]
        self.assertEqual(record["opportunity_fragment_ids"], ["frag-a", "frag-b"])
        self.assertEqual(record["support_fragment_ids"], ["frag-a"])

    def test_held_frontier_and_unexamined_budget_have_separate_states(self):
        edges = [("admit", "a", "b"), ("hold", "a", "h"), ("next", "b", "c"),
                 ("after", "c", "d")]
        evidence = [("admit", "admitted", "target_exact_bridge"),
                    ("hold", "unresolved", "opportunity_only"),
                    ("next", "admitted", "target_exact_bridge"),
                    ("after", "unresolved", "no_opportunity")]
        path = self.fixture("held", edges, evidence)
        early = TOOL.verify(path, round_budget=1, edge_budget=2)
        self.assertEqual(early["state"], "truncated")
        self.assertEqual(early["unresolved_edges"], ["hold"])
        self.assertEqual(early["unevaluated_frontier_edges"], ["next"])
        self.assertEqual(early["budget_reasons"], ["round_budget", "edge_budget"])
        self.assertEqual(early["reachable_nodes"], ["a", "b"])
        self.assertEqual(next(e["state"] for e in early["edges"] if e["id"] == "after"), "not_reached")

        wave = TOOL.verify(path, edge_budget=1)
        self.assertEqual(wave["state"], "truncated")
        self.assertEqual(wave["reachable_nodes"], ["a"])
        self.assertEqual(wave["unevaluated_frontier_edges"], ["admit", "hold"])
        self.assertEqual(wave["rounds"], [])

        done = TOOL.verify(path, round_budget=3, edge_budget=4)
        self.assertEqual(done["state"], "unresolved")
        self.assertTrue(done["frontier_exhausted"])
        self.assertEqual(done["unresolved_edges"], ["after", "hold"])
        self.assertEqual(done["budget_reasons"], [])
        self.assertEqual(done["reachable_nodes"], ["a", "b", "c"])

    def test_conflicting_certificates_and_missing_evidence_never_admit(self):
        edges = [("conflict", "a", "b"), ("missing", "a", "c")]
        evidence = [("conflict", "admitted", "target_exact_bridge"),
                    ("conflict", "excluded", "competing_context")]
        report = TOOL.verify(self.fixture("conflict", edges, evidence, duplicate=True))
        self.assertEqual(report["state"], "unresolved")
        self.assertEqual(report["reachable_nodes"], ["a"])
        self.assertEqual({e["reason"] for e in report["edges"]},
                         {"conflicting_edge_evidence", "no_edge_evidence"})

    def test_exact_budget_boundary_and_declared_scope(self):
        edges = [("one", "a", "b")]
        evidence = [("one", "admitted", "target_exact_bridge")]
        exact = TOOL.verify(self.fixture("exact", edges, evidence), round_budget=1, edge_budget=1)
        self.assertEqual(exact["state"], "converged")
        partial = TOOL.verify(self.fixture("partial", [], [], coverage="partial"), edge_budget=0)
        self.assertEqual(partial["state"], "unresolved")
        self.assertEqual(partial["unresolved_edges"], [])
        empty = TOOL.verify(self.fixture("empty", [], []), round_budget=0, edge_budget=0)
        self.assertEqual(empty["state"], "converged")

    def test_stale_read_pool_binding_and_invalid_support_fail_closed(self):
        path = self.fixture("stale", [("one", "a", "b")],
                            [("one", "admitted", "target_exact_bridge")])
        m = json.loads(path.read_text())
        pool = path.parent / "read-pool.json"
        pool.write_text("different reads\n")
        m["read_pool"]["sha256"] = TOOL.sha(pool.read_bytes())
        path.write_text(json.dumps(m))
        with self.assertRaisesRegex(ValueError, "binding mismatch"):
            TOOL.verify(path)

        path = self.fixture("bad", [("one", "a", "b")],
                            [("one", "admitted", "target_exact_bridge")])
        m = json.loads(path.read_text())
        cert = path.parent / m["edge_evidence"][0]["path"]
        body = json.loads(cert.read_text())
        body["support_fragment_ids"] = ["not_an_opportunity"]
        cert.write_text(json.dumps(body))
        m["edge_evidence"][0]["sha256"] = TOOL.sha(cert.read_bytes())
        path.write_text(json.dumps(m))
        with self.assertRaisesRegex(ValueError, "outside the opportunity set"):
            TOOL.verify(path)

    def test_cli_refuses_replacing_an_existing_report(self):
        path = self.fixture("cli", [], [])
        out = self.root / "report.json"
        command = [sys.executable, "-B", "-I", str(SCRIPT), str(path), "--out", str(out)]
        subprocess.run(command, check=True, capture_output=True)
        before = out.read_bytes()
        result = subprocess.run(command, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(out.read_bytes(), before)

    def test_small_graphs_match_independent_reachability_oracle(self):
        edges = {edge.id: edge for edge in [
            TOOL.Edge("ab", "a", "b"), TOOL.Edge("bc", "b", "c"),
            TOOL.Edge("ca", "c", "a"), TOOL.Edge("bb", "b", "b"),
        ]}
        for states in product(("admitted", "excluded", "unresolved"), repeat=len(edges)):
            frozen = dict(zip(edges, states))
            for seeds in (["a"], ["b"], ["a", "c"]):
                reached = set(seeds)
                for _ in range(len(edges)):
                    previous = set(reached)
                    reached.update(edge.destination for key, edge in edges.items()
                                   if frozen[key] == "admitted" and edge.source in previous)
                expected_held = sorted(key for key, edge in edges.items()
                                       if frozen[key] == "unresolved" and edge.source in reached)
                decision, _ = TOOL.closure(seeds, edges, {
                    key: [{"decision": state}] for key, state in frozen.items()
                }, True)
                self.assertEqual(decision["reachable_nodes"], sorted(reached))
                self.assertEqual(decision["unresolved_edges"], expected_held)
                self.assertEqual(decision["state"], "unresolved" if expected_held else "converged")


if __name__ == "__main__":
    unittest.main()
