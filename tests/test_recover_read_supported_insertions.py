import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/assembly_benchmarks/recover_read_supported_insertions.py"


class RecoveryTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.d = Path(self.tmp.name)
        ref = "A" * 100 + "C" * 100
        insert = "G" * 30
        forward = ref[:100] + insert + ref[100:]
        reverse = forward.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        (self.d / "ref.fa").write_text(">target\n" + ref + "\n")
        with (self.d / "reads.fq").open("w") as f:
            for name, seq in (("a", forward), ("b", forward), ("c", reverse)):
                f.write("@" + name + "\n" + seq + "\n+\n" + "I" * len(seq) + "\n")
        with (self.d / "map.paf").open("w") as f:
            for name, strand in (("a", "+"), ("b", "+"), ("c", "-")):
                f.write("\t".join([name,"230","0","230",strand,"target","200","0","200","200","230","60","tp:A:P","cg:Z:100M30I100M"]) + "\n")
        self.ref, self.expected = ref, forward

    def tearDown(self):
        self.tmp.cleanup()

    def command(self, out, minimum=3):
        return [sys.executable, str(SCRIPT), "--reference", str(self.d/"ref.fa"), "--reads", str(self.d/"reads.fq"), "--paf", str(self.d/"map.paf"), "--out", str(out), "--min-reads", str(minimum), "--min-identity", "0.8", "--min-flank", "50"]

    def test_recovers_forward_and_reverse_strand_without_overwrite(self):
        out = self.d / "candidate"
        cmd = self.command(out)
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        sequence = "".join((self.d/"candidate.fasta").read_text().splitlines()[1:])
        report = json.loads((self.d/"candidate.json").read_text())
        self.assertEqual(sequence, self.expected)
        self.assertEqual(report["events"][0]["read_support"], 3)
        self.assertEqual(report["events"][0]["modal_sequence_reads"], 3)
        self.assertNotEqual(subprocess.run(cmd, capture_output=True).returncode, 0)

    def test_support_threshold_prevents_correction(self):
        out = self.d / "below"
        subprocess.run(self.command(out, minimum=4), check=True, capture_output=True)
        report = json.loads((self.d/"below.json").read_text())
        self.assertEqual(report["candidate_length"], len(self.ref))
        self.assertEqual(report["events"], [])

    def test_conflicting_alignment_is_excluded_from_support(self):
        p = self.d / "map.paf"
        line = p.read_text().splitlines()[0].replace("\t+\t", "\t-\t")
        with p.open("a") as f:
            f.write(line + "\n")
        subprocess.run(self.command(self.d / "ambiguous"), check=True, capture_output=True)
        report = json.loads((self.d / "ambiguous.json").read_text())
        self.assertEqual(report["ambiguous_event_read_pairs_excluded"], 1)
        self.assertEqual(report["events"], [])

    def test_low_modal_agreement_is_reported_without_inserting(self):
        p = self.d / "reads.fq"
        content = p.read_text()
        original = "@b\n" + self.expected + "\n"
        altered = self.expected[:110] + "T" + self.expected[111:]
        self.assertIn(original, content)
        p.write_text(content.replace(original, "@b\n" + altered + "\n"))
        cmd = self.command(self.d / "lowmodal") + ["--min-modal-fraction", "0.8"]
        subprocess.run(cmd, check=True, capture_output=True)
        report = json.loads((self.d / "lowmodal.json").read_text())
        self.assertEqual(report["candidate_length"], len(self.ref))
        self.assertEqual(report["events"], [])
        self.assertEqual(report["rejected_low_modal_fraction"][0]["modal_reads"], 2)

    def test_rejects_mismatched_paf_reference(self):
        p = self.d/"map.paf"
        p.write_text(p.read_text().replace("\ttarget\t", "\twrong\t"))
        self.assertNotEqual(subprocess.run(self.command(self.d/"bad"), capture_output=True).returncode, 0)


if __name__ == "__main__":
    unittest.main()
