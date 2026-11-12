import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/audit_repeat_span_opportunity.py"
GFA = (
    b"S\tl\tTCACGT\nS\tr\tACGTACGT\nS\tt\tGTACGTT\nS\tx\tGTCCCT\nS\ty\tGTTAAA\n"
    b"L\tl\t+\tr\t+\t4M\nL\tr\t+\tr\t+\t4M\nL\tr\t+\tt\t+\t2M\n"
    b"L\tr\t+\tx\t+\t2M\nL\tt\t+\ty\t+\t3M\n"
)


class RepeatSpanOpportunityTest(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.plan = {
            "schema": "ovasm.repeat-span-opportunity-plan.v1",
            "scope": "synthetic repeat routes",
            "graph": self.file("graph.gfa", GFA),
            "raw_reads": [self.file("reads.fa", b">a\nAAAAAAAA\n>b\nACGTA\nCGTACGTA\n")],
            "routes": [{"id": "single", "path": ["l+", "r+", "x+"]}],
        }

    def file(self, name, data):
        (self.root / name).write_bytes(data)
        return {"path": name, "sha256": sha256(data).hexdigest()}

    def invoke(self):
        path = self.root / "plan.json"
        path.write_text(json.dumps(self.plan))
        return subprocess.run([sys.executable, "-B", str(SCRIPT), str(path), "--out", str(self.root / "audit.json")],
                              capture_output=True, text=True)

    def test_literal_span_in_both_orientations_and_with_repeated_visits(self):
        self.plan["routes"].extend([
            {"id": "tandem", "path": ["l+", "r+", "r+", "t+", "y+"]},
            {"id": "reverse", "path": ["y-", "t-", "r-", "r-", "l-"]},
        ])
        fastq = b"".join(b"@" + str(n).encode() + b"\r\n" + b"T" * n + b"\r\n+\r\n" + b"I" * n + b"\r\n"
                         for n in (9, 10, 18, 19))
        self.plan["raw_reads"].append(self.file("reads.fq.gz", gzip.compress(fastq)))
        run = self.invoke()
        self.assertEqual(run.returncode, 0, run.stderr)
        audit = json.loads((self.root / "audit.json").read_text())
        self.assertEqual((audit["raw_records"], audit["raw_bases"], audit["maximum_read_length"]), (6, 77, 19))
        single = audit["routes"]["single"]
        self.assertEqual(single["overlaps_bp"], [4, 2])
        self.assertEqual((single["interior_spelling_bp"], single["minimum_literal_read_length"]), (8, 10))
        self.assertEqual(single["raw_records_at_least_minimum"], 4)
        for name in ("tandem", "reverse"):
            route = audit["routes"][name]
            self.assertEqual((route["interior_spelling_bp"], route["minimum_literal_read_length"]), (17, 19))
            self.assertEqual(route["raw_records_at_least_minimum"], 1)
        self.assertEqual(audit["routes"]["reverse"]["overlaps_bp"], [3, 2, 4, 4])

    def test_unknown_route_has_no_length_claim(self):
        self.plan["routes"][0]["path"] = ["l+", "r+", "y+"]
        run = self.invoke()
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("missing or ambiguous link", run.stderr)
        self.assertFalse((self.root / "audit.json").exists())

    def test_parallel_graph_overlap_fails_closed(self):
        self.plan["graph"] = self.file("graph.gfa", GFA + b"L\tl\t+\tr\t+\t0M\n")
        run = self.invoke()
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("distinct parallel overlaps", run.stderr)
        self.assertFalse((self.root / "audit.json").exists())

    def test_changed_raw_data_fails_closed(self):
        (self.root / "reads.fa").write_bytes(b">different\nAAAA\n")
        run = self.invoke()
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("input SHA-256 mismatch", run.stderr)
        self.assertFalse((self.root / "audit.json").exists())


if __name__ == "__main__":
    unittest.main()
