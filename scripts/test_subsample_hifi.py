#!/usr/bin/env python3
import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from types import SimpleNamespace

spec = importlib.util.spec_from_file_location("subsample_hifi", Path(__file__).with_name("subsample_hifi.py"))
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class SubsamplingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.a = self.root / "a.fastq"
        self.b = self.root / "b.fastq.gz"
        # Repeated IDs across and within sources must still have independent record keys.
        payload = b"".join(b"@duplicate\nACGT\nAC\n+\nIII\nIII\n" for _ in range(200))
        self.a.write_bytes(payload)
        with gzip.open(self.b, "wb") as f:
            f.write(payload)

    def tearDown(self):
        self.temp.cleanup()

    def run_subsets(self, name, reads=None, compresslevel=1):
        out = self.root / name
        argv = []
        for path in reads or [self.a, self.b]:
            argv += ["--reads", str(path)]
        argv += ["--source-id", "runA", "--source-id", "runB",
                 "--fraction", ".2", "--fraction", ".6", "--fraction", "1",
                 "--seed", "7", "--seed", "8", "--out", str(out),
                 "--compresslevel", str(compresslevel)]
        mod.main(argv)
        return out, json.loads((out / "manifest.json").read_text())

    def test_nested_deterministic_relocated_and_complete_records(self):
        first, report = self.run_subsets("first")
        second, again = self.run_subsets("second")
        moved_a = self.root / "different-name.fq"
        moved_b = self.root / "different-name.fq.gz"
        moved_a.write_bytes(self.a.read_bytes())
        moved_b.write_bytes(self.b.read_bytes())
        moved, _ = self.run_subsets("moved", [moved_a, moved_b])
        self.assertEqual(report["status"], "complete")
        self.assertFalse((first / "INCOMPLETE").exists())
        by_seed = {}
        for entry in report["subsets"]:
            name = entry["path"]
            self.assertEqual((first / name).read_bytes(), (second / name).read_bytes())
            self.assertEqual((first / name).read_bytes(), (moved / name).read_bytes())
            with gzip.open(first / name, "rb") as handle:
                rows = list(mod.records(handle))
            self.assertEqual(len(rows), entry["reads"])
            self.assertEqual(sum(r[2] for r in rows), entry["bases"])
            self.assertTrue(all(r[2] == 6 for r in rows))
            by_seed.setdefault(entry["seed"], []).append(entry["reads"])
        for sizes in by_seed.values():
            self.assertEqual(sizes, sorted(sizes))
            self.assertEqual(sizes[-1], 400)
        self.assertNotEqual(by_seed[7], by_seed[8])
        # Independently check membership nesting with distinguishable read contents.
        self.a.write_bytes(b"".join(f"@r{i}\nACGT\n+\nIIII\n".encode() for i in range(200)))
        unique, unique_report = self.run_subsets("unique")
        sets = []
        for entry in unique_report["subsets"][:3]:
            with gzip.open(unique / entry["path"], "rb") as handle:
                sets.append({r[1] for r in mod.records(handle)})
        self.assertTrue(sets[0] <= sets[1] <= sets[2])

    def test_missing_final_newline_does_not_fuse_sources(self):
        self.a.write_bytes(b"@last\nACGT\n+\nIIII")
        out, report = self.run_subsets("no-newline")
        entry = next(e for e in report["subsets"] if e["seed"] == 7 and e["fraction"] == "1")
        with gzip.open(out / entry["path"], "rb") as handle:
            rows = list(mod.records(handle))
        self.assertEqual(len(rows), 201)
        self.assertEqual(rows[0][1], b"last")

    def test_low_disk_space_leaves_incomplete_manifest(self):
        with patch.object(mod.shutil, "disk_usage", return_value=SimpleNamespace(free=0)):
            with self.assertRaisesRegex(OSError, "below --min-free-gib"):
                self.run_subsets("no-space")
        out = self.root / "no-space"
        self.assertTrue((out / "INCOMPLETE").exists())
        report = json.loads((out / "manifest.json").read_text())
        self.assertEqual(report["status"], "incomplete")
        self.assertEqual(report["min_free_gib"], 30)
        self.assertIn("outputs incomplete", report["error"])

    def test_disk_space_is_rechecked_during_input(self):
        self.a.write_bytes(b"@a\nACGT\n+\nIIII\n" * 1100)
        with patch.object(mod.shutil, "disk_usage", side_effect=[
                SimpleNamespace(free=100 * (1 << 30)), SimpleNamespace(free=0)]) as disk:
            with self.assertRaisesRegex(OSError, "below --min-free-gib"):
                self.run_subsets("space-ran-out")
            self.assertEqual(disk.call_count, 2)
        self.assertTrue((self.root / "space-ran-out" / "INCOMPLETE").exists())

    def test_compression_level_does_not_change_selected_records(self):
        fast, fast_report = self.run_subsets("fast", compresslevel=1)
        dense, dense_report = self.run_subsets("dense", compresslevel=9)
        self.assertEqual(fast_report["compresslevel"], 1)
        self.assertEqual(dense_report["compresslevel"], 9)
        for a, b in zip(fast_report["subsets"], dense_report["subsets"]):
            self.assertEqual(a["uncompressed_fastq_sha256"], b["uncompressed_fastq_sha256"])
            self.assertEqual(a["reads"], b["reads"])
            self.assertEqual(gzip.decompress((fast / a["path"]).read_bytes()),
                             gzip.decompress((dense / b["path"]).read_bytes()))

    def test_refuses_overwrite(self):
        self.run_subsets("same")
        with self.assertRaises(FileExistsError):
            self.run_subsets("same")

    def test_invalid_quality_marks_incomplete(self):
        self.a.write_bytes(b"@bad\nACGT\n+\nII\n")
        with self.assertRaises(ValueError):
            self.run_subsets("invalid")
        out = self.root / "invalid"
        self.assertTrue((out / "INCOMPLETE").exists())
        report = json.loads((out / "manifest.json").read_text())
        self.assertEqual(report["status"], "incomplete")
        self.assertIn("truncated qualities", report["error"])


if __name__ == "__main__":
    unittest.main()
