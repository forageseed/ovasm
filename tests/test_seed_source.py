import gzip
from hashlib import sha256
import json
from pathlib import Path
import random
import runpy
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_seed_source.py"
TOOL = runpy.run_path(str(SCRIPT))


class SeedSourceAudit(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        rng = random.Random(412)
        self.target = bytes(rng.choices(b"ACGT", k=64))
        self.competing = bytes(rng.choices(b"ACGT", k=64))
        self.novel = bytes(rng.choices(b"ACGT", k=64))
        self.proof = self.write("proof.json", b'{"source_evidence": "control"}\n')
        self.target_ref = self.write("mt.fa", b">mt\n" + self.target + b"\n")
        self.competing_ref = self.write("nuclear.fa", b">chr2\n" + self.competing + b"\n")
        self.manifest = {
            "schema": TOOL["MANIFEST"], "scope": "frozen synthetic source comparison",
            "references": [dict(self.ref(self.target_ref), id="mt", role="target"),
                           dict(self.ref(self.competing_ref), id="chr2", role="competing")],
            "windows": [self.window("t", "mt", "mt"), self.window("c", "chr2", "chr2")],
        }

    def write(self, name, data):
        path = self.root / name
        path.write_bytes(data)
        return path

    def ref(self, path):
        return {"path": path.name, "sha256": sha256(path.read_bytes()).hexdigest()}

    def window(self, name, reference, record):
        return {"id": name, "reference": reference, "record": record,
                "interval": [0, 64], "boundary": 32, "provenance": self.ref(self.proof)}

    def assemble(self, originals, seeds, as_gzip=False):
        if as_gzip:
            data = b"".join(
                b"@" + name + b" notes\r\n" + seq + b"\r\n+\r\n" + b"I" * len(seq) + b"\r\n"
                for name, seq in originals
            )
            raw = self.write("raw.fq.gz", gzip.compress(data, mtime=0))
        else:
            raw = self.write("raw.fa", b"".join(
                b">" + name + b" original\n" + seq + b"\n" for name, seq in originals))
        seed = self.write("seeds.fa", b"".join(seeds))
        self.manifest.update(raw_reads=[self.ref(raw)], seeds=self.ref(seed))
        return self.write("manifest.json", json.dumps(self.manifest).encode())

    @staticmethod
    def crop(name, sequence, start=0, width=None):
        width = len(sequence) if width is None else width
        return (b">" + name + b" whole_read_scan span=" + f"{start}-{start + width}".encode()
                + b"\r\n" + sequence[:32] + b"\r\n" + sequence[32:] + b"\r\n")

    def test_certified_competitor_only_is_removed(self):
        reverse = self.competing.translate(TOOL["COMPLEMENT"])[::-1]
        mixed = self.competing + self.target
        originals = [(b"target", self.target), (b"competing", self.competing),
                     (b"reverse", reverse), (b"unknown", self.novel), (b"mixed", mixed),
                     (b"outside", self.competing + self.novel), (b"changed", self.target),
                     (b"repeated", self.competing), (b"repeated", self.competing),
                     (b"multiple", self.competing * 2)]
        crops = [self.crop(name, seq) for name, seq in originals[:5]]
        crops += [self.crop(b"outside", self.novel, start=64),
                  self.crop(b"changed", self.competing), self.crop(b"repeated", self.competing),
                  self.crop(b"multiple", self.competing * 2), self.crop(b"missing", self.novel)]
        path = self.assemble(originals, crops, as_gzip=True)
        report, retained = TOOL["audit"](path)
        self.assertEqual(report["counts"], {"target_supported": 1, "competing_supported": 2,
                                            "unresolved": 7})
        self.assertEqual(retained, [crops[0], *crops[3:]])
        rows = {row["read_id"]: row for row in report["records"]}
        self.assertEqual(rows["reverse"]["observations"][0]["original_positions"], [(0, 64, "-")])
        self.assertEqual(rows["mixed"]["reason"], "mixed_source_contexts")
        self.assertEqual(rows["outside"]["reason"], "source_window_ambiguous_or_outside_crop")
        self.assertEqual(rows["changed"]["reason"], "crop_does_not_bind_to_original")
        self.assertEqual(rows["repeated"]["reason"], "original_identifier_repeated")
        self.assertEqual(rows["multiple"]["reason"], "source_window_ambiguous_or_outside_crop")
        self.assertEqual(rows["missing"]["reason"], "original_missing")
        out, view = self.root / "audit.json", self.root / "retained.fa"
        subprocess.run([sys.executable, "-B", "-I", str(SCRIPT), str(path),
                        "--out", str(out), "--out-fasta", str(view)], check=True)
        self.assertEqual(view.read_bytes(), b"".join(retained))
        self.assertEqual(json.loads(out.read_text())["retained_fasta_sha256"], self.ref(view)["sha256"])

    def test_explicit_original_identifier_field_supports_archive_prefix(self):
        path = self.assemble([(b"SRR123.1 movie/1", self.competing)],
                             [self.crop(b"movie/1", self.competing)])
        report, _ = TOOL["audit"](path)
        self.assertEqual(report["records"][0]["reason"], "original_missing")
        self.manifest["raw_reads"][0]["id_field"] = 2
        path.write_text(json.dumps(self.manifest))
        report, retained = TOOL["audit"](path)
        self.assertEqual(report["records"][0]["state"], "competing_supported")
        self.assertEqual(retained, [])

    def test_ambiguous_panel_window_cannot_remove_seed(self):
        alias = self.write("plastid.fa", b">pt\n" + self.competing + b"\n")
        self.manifest["references"].append(dict(self.ref(alias), id="pt", role="competing"))
        crop = self.crop(b"seen", self.competing)
        path = self.assemble([(b"seen", self.competing)], [crop])
        report, retained = TOOL["audit"](path)
        self.assertFalse(report["windows"]["c"]["unique_in_named_panel"])
        self.assertEqual(report["records"][0]["state"], "unresolved")
        self.assertEqual(retained, [crop])

    def test_provenance_mismatch_leaves_no_outputs(self):
        path = self.assemble([(b"seen", self.competing)], [self.crop(b"seen", self.competing)])
        self.proof.write_bytes(b"modified\n")
        out, view = self.root / "audit.json", self.root / "retained.fa"
        result = subprocess.run([sys.executable, "-B", "-I", str(SCRIPT), str(path),
                                 "--out", str(out), "--out-fasta", str(view)], capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"SHA-256 mismatch", result.stderr)
        self.assertFalse(out.exists())
        self.assertFalse(view.exists())

    def test_publish_rolls_back_if_second_output_appears(self):
        out = self.root / "audit.json"
        occupied = self.write("existing.fa", b"keep me\n")
        with self.assertRaises(FileExistsError):
            TOOL["publish"]([(out, [b"new"]), (occupied, [b"bad"])])
        self.assertFalse(out.exists())
        self.assertEqual(occupied.read_bytes(), b"keep me\n")

    def test_source_overlap_requires_original_even_when_seed_matches(self):
        seed = self.crop(b"seen", self.competing)
        path = self.assemble([(b"seen", self.competing), (b"extra", self.target)], [seed])
        extra = self.write("second.fa", b">seen\n" + self.competing + b"\n")
        self.manifest["raw_reads"].append(self.ref(extra))
        path.write_text(json.dumps(self.manifest))
        report, retained = TOOL["audit"](path)
        self.assertEqual(report["records"][0]["reason"], "original_identifier_repeated")
        self.assertEqual(retained, [seed])


if __name__ == "__main__":
    unittest.main()
