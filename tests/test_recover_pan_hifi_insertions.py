import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/assembly_benchmarks/recover_pan_hifi_insertions.py"


class PanRecoveryTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = Path(self.tmp.name)
        self.reference = "A" * 500 + "C" * 500
        self.expected = self.reference[:500] + "G" * 30 + self.reference[500:]
        (self.d / "pan.fa").write_text(">other\nAAAA\n>chosen description\n" + self.reference + "\n")
        (self.d / "reads.fq").write_text("".join("@r%d\n%s\n+\n%s\n" % (i, self.expected, "I" * len(self.expected)) for i in range(3)))
        self.mapper = self.d / "mapper"
        self.mapper.write_text("#!/usr/bin/env python3\n" + "for i in range(3): print('r%d\t1030\t0\t1030\t+\tselected_pan_representative\t1000\t0\t1000\t1000\t1030\t60\ttp:A:P\tcg:Z:500M30I500M' % i)\n")
        self.mapper.chmod(0o755)

    def tearDown(self):
        self.tmp.cleanup()

    def cmd(self, record="chosen"):
        return [sys.executable, "-I", str(SCRIPT), "--pan-fasta", str(self.d / "pan.fa"), "--record", record, "--hifi", str(self.d / "reads.fq"), "--out-dir", str(self.d / "out"), "--minimap2", str(self.mapper), "--min-reads", "3"]

    def test_explicit_record_recovers_and_refuses_overwrite(self):
        z = subprocess.run(self.cmd(), capture_output=True, text=True)
        self.assertEqual(z.returncode, 0, z.stderr)
        result = json.loads((self.d / "out/workflow.json").read_text())
        candidate = "".join((self.d / "out/candidate.fasta").read_text().splitlines()[1:])
        self.assertEqual(candidate, self.expected)
        self.assertEqual(result["event_count"], 1)
        self.assertIn("not independent", result["scope"])
        self.assertNotEqual(subprocess.run(self.cmd(), capture_output=True).returncode, 0)
        self.assertEqual("".join((self.d / "out/candidate.fasta").read_text().splitlines()[1:]), self.expected)

    def test_missing_or_duplicate_record_fails_before_output(self):
        self.assertNotEqual(subprocess.run(self.cmd("missing"), capture_output=True).returncode, 0)
        self.assertFalse((self.d / "out").exists())
        with (self.d / "pan.fa").open("a") as f:
            f.write(">chosen\nACGT\n")
        self.assertNotEqual(subprocess.run(self.cmd(), capture_output=True).returncode, 0)
        self.assertFalse((self.d / "out").exists())


if __name__ == "__main__":
    unittest.main()
