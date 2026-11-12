import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/extract_seed_crop_originals.py"


def fq(header, seq, qual):
    return b"@" + header + b"\r\n" + seq + b"\r\n+ extra\r\n" + qual + b"\r\n"


class CropOriginalsTest(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.report = self.root / "report.json"
        self.originals = self.root / "originals.fq"
        self.plan = {
            "schema": "ovasm.seed-crop-originals-plan.v1", "scope": "synthetic originals",
            "candidate_sets": [], "raw_reads": [],
        }

    def file(self, name, data):
        (self.root / name).write_bytes(data)
        return {"path": name, "sha256": sha256(data).hexdigest()}

    def add_set(self, name, fasta):
        self.plan["candidate_sets"].append({"id": name} | self.file(name + ".fa", fasta))

    def run_cli(self):
        plan = self.root / "plan.json"
        plan.write_text(json.dumps(self.plan))
        return subprocess.run([sys.executable, "-B", str(SCRIPT), str(plan), "--out", str(self.report),
                              "--out-reads", str(self.originals)], capture_output=True, text=True)

    def test_exact_full_records_and_union_with_late_duplicate(self):
        self.add_set("old", b">a cluster=1 span=1-4\nCgT\n>b span=1-4\nGGC\n")
        self.add_set("new", b">a cluster=1 span=1-4\nCgT\n>c span=0-4\nTGAA\n"
                     b">d span=0-4\nTAGG\n>e span=7-9\nTC\n>f span=0-2\nCA\n")
        header = b"a m123/45/ccs full header"
        record = fq(header, b"aCgTAA", b"#$%&'(")
        self.plan["raw_reads"] = [
            self.file("orig.fq", record + fq(b"b movie/1", b"AGGC", b"1234")
                      + fq(b"d movie/2", b"TAAA", b"4567") + fq(b"e movie/3", b"GGC", b"789")),
            self.file("second.fq.gz", gzip.compress(
                fq(b"archive b unrelated", b"AGGC", b"0000")
                + fq(b"archive c another", b"TGAA", b"IIII")
                + fq(b"archive a same candidate name other bases", b"CGTT", b"3333")
            )) | {"id_field": 2},
        ]
        run = self.run_cli()
        self.assertEqual(run.returncode, 0, run.stderr)
        out = json.loads(self.report.read_bytes())
        expected = (b"@a m123/45/ccs full header\naCgTAA\n+\n#$%&'(\n"
                    b"@archive c another\nTGAA\n+\nIIII\n")
        self.assertEqual(self.originals.read_bytes(), expected)
        self.assertEqual(out["output_fastq"]["sha256"], sha256(expected).hexdigest())
        self.assertEqual(out["candidate_counts"], {
            "bound": 3, "multiple_exact_originals": 1, "crop_bytes_mismatch": 1,
            "span_outside_originals": 1, "identifier_not_observed": 1,
        })
        self.assertEqual([row["records"] for row in out["raw_reads"]], [4, 3])
        self.assertEqual([row["original_output_record_index0"] for row in out["candidates"]], [0, None, 0, 1, None, None, None])
        self.assertEqual(out["originals"][0]["read_id_hex"], header.hex())
        self.assertEqual(out["originals"][0]["quality_sha256"], sha256(b"#$%&'(").hexdigest())
        self.assertEqual(out["originals"][1]["input_index0"], 1)
        self.assertEqual(out["originals"][1]["record_index0"], 1)

    def test_trailing_carriage_returns_round_trip(self):
        self.add_set("candidate", b">a span=0-4\nTGAA\n")
        raw = fq(b"a note\r", b"TGAA\r", b"IIII\r")
        self.plan["raw_reads"] = [self.file("reads.fq", raw)]
        run = self.run_cli()
        self.assertEqual(run.returncode, 0, run.stderr)
        output = self.originals.read_bytes()
        self.assertEqual(output, b"@a note\r\r\nTGAA\r\r\n+\nIIII\r\r\n")
        audit = json.loads(self.report.read_bytes())
        row, = audit["originals"]
        self.assertEqual(row["read_id_hex"], b"a note\r".hex())
        self.assertEqual(row["sequence_sha256"], sha256(b"TGAA\r").hexdigest())
        self.assertEqual(row["quality_sha256"], sha256(b"IIII\r").hexdigest())
        sys.path.insert(0, str(SCRIPT.parent))
        from prepare_read_admission import LockedFile, fastq
        source = LockedFile.lock(self.root, {"path": "originals.fq", "sha256": sha256(output).hexdigest()})
        self.assertEqual(list(fastq(source)), [(b"a note\r", b"TGAA\r", b"IIII\r")])

    def test_wrong_compressed_digest_publishes_neither_output(self):
        self.add_set("candidate", b">a span=0-4\nTGAA\n")
        binding = self.file("reads.fq.gz", gzip.compress(fq(b"a", b"TGAA", b"IIII")))
        self.plan["raw_reads"] = [binding | {"sha256": "0" * 64}]
        run = self.run_cli()
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("original SHA-256 or byte count mismatch", run.stderr)
        self.assertFalse(self.report.exists())
        self.assertFalse(self.originals.exists())

    def test_invalid_crop_span_publishes_neither_output(self):
        self.add_set("candidate", b">a span=0-5\nTGAA\n")
        self.plan["raw_reads"] = [self.file("reads.fq", fq(b"a", b"TGAA", b"IIII"))]
        run = self.run_cli()
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("candidate span disagrees", run.stderr)
        self.assertFalse(self.report.exists())
        self.assertFalse(self.originals.exists())

    def test_without_matches_publishes_report_only(self):
        self.add_set("candidate", b">a span=0-4\nTGAA\n")
        self.plan["raw_reads"] = [self.file("reads.fq", fq(b"b", b"TGAA", b"IIII"))]
        run = self.run_cli()
        self.assertEqual(run.returncode, 0, run.stderr)
        self.assertEqual(json.loads(self.report.read_bytes())["output_records"], 0)
        self.assertFalse(self.originals.exists())


if __name__ == "__main__":
    unittest.main()
