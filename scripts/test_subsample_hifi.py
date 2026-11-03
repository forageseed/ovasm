#!/usr/bin/env python3
import gzip
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

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

    def run_subsets(self, name, reads=None):
        out = self.root / name
        argv = []
        for path in reads or [self.a, self.b]:
            argv += ["--reads", str(path)]
        argv += ["--source-id", "runA", "--source-id", "runB",
                 "--fraction", ".2", "--fraction", ".6", "--fraction", "1",
                 "--seed", "7", "--seed", "8", "--out", str(out)]
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
