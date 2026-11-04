import importlib.util
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/assembly_benchmarks/audit_missing_graph_links.py"
spec = importlib.util.spec_from_file_location("audit_missing_graph_links", SCRIPT)
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)

def paf(read, qstart, qend, strand, node, tstart=0, tend=600, identity=600):
    return "\t".join(map(str, [read, 1226, qstart, qend, strand, node, 600, tstart, tend, identity, 600, 60])) + "\n"

class MissingLinks(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.graph = self.root / "graph.gfa"
        self.graph.write_text("".join(f"S\t{name}\t{'A'*600}\n" for name in ("a", "b", "c")) + "L\ta\t+\tb\t+\t0M\n")
        self.paf = self.root / "reads.paf"

    def tearDown(self):
        self.tmp.cleanup()

    def test_missing_link_reverse_complement_and_duplicates(self):
        lines = [paf("r1", 0, 600, "+", "a"), paf("r1", 626, 1226, "+", "c"),
                 paf("r1", 626, 1226, "+", "c"),
                 paf("r2", 0, 600, "-", "c"), paf("r2", 626, 1226, "-", "a"),
                 paf("r3", 0, 600, "+", "a"), paf("r3", 626, 1226, "+", "b")]
        self.paf.write_text("".join(lines))
        out = model.audit(self.graph, self.paf)
        self.assertEqual(len(out["candidates"]), 1)
        c = out["candidates"][0]
        self.assertEqual((c["from"], c["from_orientation"], c["to"], c["to_orientation"]), ("a", "+", "c", "+"))
        self.assertEqual(c["distinct_reads"], 2)
        self.assertEqual((c["gap_min"], c["gap_max"]), (26, 26))
        self.assertFalse(c["near_gap_limit"])

    def test_nonterminal_and_low_identity_alignments_excluded(self):
        self.paf.write_text(paf("r1", 0, 400, "+", "a", 0, 400, 400) + paf("r1", 426, 1026, "+", "c") +
                            paf("r2", 0, 600, "+", "a", identity=400) + paf("r2", 626, 1226, "+", "c"))
        self.assertEqual(model.audit(self.graph, self.paf)["candidates"], [])

    def test_gap_limit_warning_and_conflicting_query_length(self):
        left = "r\t1200\t0\t500\t+\ta\t600\t100\t600\t500\t500\t60\n"
        right = "r\t1200\t600\t1100\t+\tc\t600\t0\t500\t500\t500\t60\n"
        self.paf.write_text(left + right)
        candidate = model.audit(self.graph, self.paf)["candidates"][0]
        self.assertEqual(candidate["gap_max"], 100)
        self.assertTrue(candidate["near_gap_limit"])
        self.paf.write_text(left + right.replace("r\t1200", "r\t1300", 1))
        with self.assertRaisesRegex(ValueError, "inconsistent lengths"):
            model.audit(self.graph, self.paf)

    def test_rejects_inconsistent_graph_and_paf(self):
        self.paf.write_text(paf("r1", 0, 600, "+", "missing"))
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            model.audit(self.graph, self.paf)

if __name__ == "__main__":
    unittest.main()
