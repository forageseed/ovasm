"""Regression controls for exact, provenance-bound read source screening."""
import gzip
import hashlib
import json
import random
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "prepare_read_admission.py"
ROLES = ("mitochondrion", "plastid", "nuclear")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def fasta(name, seq):
    return b">" + name.encode() + b"\n" + seq + b"\n"


def fastq(rid, seq, qual):
    return b"@" + rid + b"\n" + seq + b"\n+\n" + qual + b"\n"


def rc(seq):
    return seq.translate(bytes.maketrans(b"ACGT", b"TGCA"))[::-1]


class PrepareReadAdmission(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        rand = random.Random(81)
        self.seq = bytes(rand.choice(b"ACGT") for _ in range(96))
        self.window = {"id": "junction", "reference": "mt", "start0": 31, "boundary0": 46, "end0": 63}
        self.rid = b"original/2 provenance tag"
        self.qual = b"*" * len(self.seq)
        self.origin = self.file("geometry.json", b'{"scope":"independent synthetic paired geometry"}\n')
        self.panel = {
            "mitochondrion": self.file("mt.fa", fasta("mt", self.seq)),
            "plastid": self.file("pt.fa.gz", gzip.compress(fasta("pt", b"T" * 160))),
            "nuclear": self.file("nuclear.fa", fasta("n", b"C" * 160)),
        }
        self.case_no = 0

    def file(self, name, data):
        (self.root / name).write_bytes(data)
        return {"path": name, "sha256": sha(data)}

    def qualify(self, rid, seq, qual, **overrides):
        row = {"read_id_hex": rid.hex(), "sequence_sha256": sha(seq),
               "quality_sha256": sha(qual), "window": self.window,
               "verdict": "target_supported"}
        return row | overrides

    def screen(self, *, rows=None, raw=None, cap=None, proof_updates=None, plan_updates=None):
        self.case_no += 1
        stem = f"case{self.case_no}"
        if raw is None:
            raw = fastq(self.rid, self.seq, self.qual)
        source = self.file(f"{stem}.reads.gz", gzip.compress(raw))
        plan = {"schema": "ovasm.read-source-screen.v1", "scope": "synthetic source screen",
                "target": "mitochondrion", "reads": [source], "panel": self.panel,
                "windows": [self.window]}
        if rows is not None:
            proof = {"schema": "ovasm.read-source-qualification.v1",
                     "scope": "same original paired geometry", "target": "mitochondrion",
                     "panel_sha256": {role: self.panel[role]["sha256"] for role in ROLES},
                     "provenance": self.origin, "records": rows}
            if proof_updates:
                proof.update(proof_updates)
            plan["qualification"] = self.file(f"{stem}.qual.json", json.dumps(proof).encode())
        if plan_updates:
            plan.update(plan_updates)
        spec = self.file(f"{stem}.json", json.dumps(plan).encode())
        out = self.root / f"{stem}.out"
        args = [sys.executable, "-B", "-I", str(SCRIPT), str(self.root / spec["path"]), "--out", str(out)]
        if cap is not None:
            args += ["--max-candidates", str(cap)]
        run = subprocess.run(args, capture_output=True, text=True, cwd=self.root, check=False)
        return run, out

    def result(self, run, out):
        self.assertEqual(run.returncode, 0, run.stderr)
        data = json.loads((out / "source-screen.json").read_text())
        manifest_path = out / "read-admission.json"
        return data, json.loads(manifest_path.read_text()) if manifest_path.exists() else None

    def test_native_manifest_binds_forward_reverse_and_full_header(self):
        rid2 = b"reversed/1 instrument metadata"
        raw = fastq(self.rid, self.seq.lower(), self.qual) + fastq(rid2, rc(self.seq), self.qual)
        rows = [self.qualify(self.rid, self.seq.lower(), self.qual), self.qualify(rid2, rc(self.seq), self.qual)]
        run, out = self.screen(rows=rows, raw=raw)
        report, native = self.result(run, out)
        self.assertEqual(report["decision_counts"], {"target_supported": 2})
        self.assertEqual(native["schema"], "ovasm.read-admission.v1")
        self.assertEqual([r["read_id_hex"] for r in native["records"]], sorted(r["read_id_hex"] for r in rows))
        for row in native["records"]:
            proof = Path(row["provenance"]["path"])
            self.assertEqual(proof, out / "source-screen.json")
            self.assertEqual(row["provenance"]["sha256"], sha(proof.read_bytes()))
        self.assertEqual([r["full_record_placements"]["mitochondrion"][0]["strand"] for r in report["records"]],
                         ["+", "-"])

    def test_missing_stale_or_explicit_unresolved_qualification(self):
        base = self.qualify(self.rid, self.seq, self.qual)
        variants = [(None, "missing"), ([base | {"quality_sha256": "1" * 64}], "record_mismatch"),
                    ([base | {"verdict": "unresolved"}], "unresolved")]
        for rows, expected in variants:
            with self.subTest(expected=expected):
                report, native = self.result(*self.screen(rows=rows))
                self.assertIsNone(native)
                decision = report["records"][0]
                self.assertEqual(decision["status"], "unresolved_qualification")
                self.assertEqual(decision["windows"]["junction"]["qualification"], expected)

    def test_shared_exact_window_vetoes_a_unique_whole_record(self):
        word = self.seq[self.window["start0"]:self.window["end0"]]
        self.panel["plastid"] = self.file("pt.fa", fasta("p", rc(word)))
        report, native = self.result(*self.screen(rows=[self.qualify(self.rid, self.seq, self.qual)]))
        self.assertIsNone(native)
        row = report["records"][0]
        self.assertEqual(row["status"], "unresolved_exact_placement")
        self.assertFalse(row["full_record_placements"]["plastid"])
        self.assertEqual(report["windows_exact_placements"]["junction"]["plastid"][0]["strand"], "-")

    def test_shared_full_record_and_same_role_duplicate_are_withheld(self):
        proof = [self.qualify(self.rid, self.seq, self.qual)]
        self.panel["nuclear"] = self.file("nuclear.fa", fasta("n", self.seq))
        report, native = self.result(*self.screen(rows=proof))
        self.assertIsNone(native)
        self.assertEqual(report["records"][0]["status"], "competing_exact_record")
        self.panel["nuclear"] = self.file("nuclear.fa", fasta("n", b"C" * 160))
        self.panel["mitochondrion"] = self.file("mt.fa", fasta("mt", self.seq) + fasta("mt2", self.seq))
        report, native = self.result(*self.screen(rows=proof))
        self.assertIsNone(native)
        self.assertEqual(report["records"][0]["status"], "unresolved_exact_placement")
        self.assertEqual(len(report["records"][0]["full_record_placements"]["mitochondrion"]), 2)

    def test_cap_truncation_publishes_no_manifest(self):
        raw = fastq(self.rid, self.seq, self.qual) + fastq(b"second", self.seq, self.qual)
        proof = [self.qualify(self.rid, self.seq, self.qual)]
        report, native = self.result(*self.screen(rows=proof, raw=raw, cap=1))
        self.assertIsNone(native)
        self.assertEqual(report["state"], "truncated")
        self.assertEqual(report["reason"], "candidate_cap")

    def test_duplicate_original_header_rejected_even_when_one_record_misses_window(self):
        raw = fastq(self.rid, self.seq, self.qual) + fastq(self.rid, b"C" * len(self.seq), self.qual)
        run, out = self.screen(rows=[self.qualify(self.rid, self.seq, self.qual)], raw=raw)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("duplicate candidate FASTQ header", run.stderr)
        self.assertFalse(out.exists())

    def test_mutated_provenance_and_coordinate_mismatch_are_errors(self):
        self.root.joinpath(self.origin["path"]).write_text("mutated")
        run, out = self.screen(rows=[self.qualify(self.rid, self.seq, self.qual)])
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("SHA-256 mismatch", run.stderr)
        self.assertFalse(out.exists())
        self.origin = self.file("geometry.json", b"bound")
        stale = self.qualify(self.rid, self.seq, self.qual,
                             window=self.window | {"boundary0": self.window["boundary0"] + 1})
        run, out = self.screen(rows=[stale])
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("qualification junction window differs", run.stderr)
        self.assertFalse(out.exists())


if __name__ == "__main__":
    unittest.main()
