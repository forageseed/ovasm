import gzip
from hashlib import sha256
import json
import os
from pathlib import Path
import random
import runpy
import shutil
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
            "target": "mitochondrion",
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
                        "--out", str(out), "--out-fasta", str(view), "--python-placements"], check=True)
        self.assertEqual(view.read_bytes(), b"".join(retained))
        result = json.loads(out.read_text())
        self.assertEqual(result["target"], "mitochondrion")
        self.assertEqual(result["retained_fasta_sha256"], self.ref(view)["sha256"])

    def test_strict_binds_distinct_crops_with_reused_ids_and_excludes_outside_evidence(self):
        crops = [self.crop(b"same", self.target), self.crop(b"same", self.competing),
                 self.crop(b"ambiguous", self.target), self.crop(b"outside", self.novel, start=64),
                 self.crop(b"mixed", self.target)]
        path = self.assemble([(b"same", self.target), (b"same", self.competing),
                             (b"ambiguous", self.target + self.novel),
                             (b"ambiguous", self.target + self.competing),
                             (b"outside", self.target + self.novel),
                             (b"mixed", self.target + self.competing)], crops, as_gzip=True)
        legacy, _ = TOOL["audit"](path)
        self.assertNotIn("retention_policy", legacy)
        self.assertEqual(legacy["records"][0]["reason"], "original_identifier_repeated")
        report, retained = TOOL["audit"](path, TOOL["TARGET_ONLY"])
        self.assertEqual(report["retention_policy"], "target_supported_only")
        self.assertEqual(report["counts"], {"target_supported": 1, "competing_supported": 1, "unresolved": 3})
        self.assertEqual(retained, [crops[0]])
        self.assertEqual(report["records"][0]["original_record_index"], 0)
        self.assertEqual(report["records"][1]["original_record_index"], 1)
        self.assertEqual(report["records"][2]["reason"], "crop_matches_multiple_original_records")
        self.assertEqual(report["records"][3]["reason"], "source_window_ambiguous_or_outside_crop")
        self.assertEqual(report["records"][4]["reason"], "mixed_source_contexts")
        out, view = self.root / "strict.json", self.root / "strict.fa"
        subprocess.run([sys.executable, "-B", "-I", str(SCRIPT), str(path),
                        "--retention-policy", "target-supported-only", "--out", str(out),
                        "--out-fasta", str(view), "--python-placements"], check=True)
        self.assertEqual(view.read_bytes(), crops[0])

    def test_strict_empty_is_reportable_and_repeated_input_alias_fails(self):
        path = self.assemble([(b"unknown", self.novel)], [self.crop(b"unknown", self.novel)])
        report, retained = TOOL["audit"](path, TOOL["TARGET_ONLY"])
        self.assertEqual(report["retained_fasta_sha256"], sha256(b"").hexdigest())
        self.assertEqual(retained, [])
        self.manifest["raw_reads"].append(dict(self.manifest["raw_reads"][0]))
        path.write_text(json.dumps(self.manifest))
        with self.assertRaisesRegex(ValueError, "repeats an original input file"):
            TOOL["audit"](path, TOOL["TARGET_ONLY"])

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
                                 "--out", str(out), "--out-fasta", str(view), "--python-placements"],
                                capture_output=True)
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


    def test_default_requires_an_ovasm_executable_and_does_not_fall_back(self):
        path = self.assemble([(b"seen", self.target)], [self.crop(b"seen", self.target)])
        bin_dir = self.root / "bin"
        bin_dir.mkdir()
        out, view = self.root / "audit.json", self.root / "retained.fa"
        command = [sys.executable, "-B", "-I", str(SCRIPT), str(path), "--out", str(out), "--out-fasta", str(view)]
        env = {**os.environ, "PATH": str(bin_dir)}
        cases = (
            ([], b"OVASM executable not found on PATH"),
            (["--exact-placement-bin", str(bin_dir / "missing")], b"cannot resolve OVASM executable"),
            (["--exact-placement-bin", str(bin_dir)], b"not an executable file"),
            (["--exact-placement-bin", str(bin_dir), "--python-placements"], b"not allowed with argument"),
        )
        for extra, message in cases:
            with self.subTest(extra=extra):
                failed = subprocess.run([*command, *extra], env=env, capture_output=True, timeout=5)
                self.assertNotEqual(failed.returncode, 0)
                self.assertIn(message, failed.stderr)
                self.assertFalse(out.exists())
                self.assertFalse(view.exists())

        if os.name != "nt":
            fake = bin_dir / "ovasm"
            fake.write_text(f"#!{sys.executable}\nraise SystemExit(9)\n")
            fake.chmod(0o700)
            failed = subprocess.run(command, env=env, capture_output=True, timeout=5)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn(b"OVASM exact placement backend failed (panel", failed.stderr)
            self.assertFalse(out.exists())
            self.assertFalse(view.exists())

        subprocess.run([*command, "--python-placements"], env=env, check=True, timeout=5)
        self.assertTrue(out.exists())
        self.assertTrue(view.exists())

    def test_requested_batch_failure_leaves_no_outputs(self):
        path = self.assemble([(b"seen", self.target)], [self.crop(b"seen", self.target)])
        prefix = b"OVXPOS1\tpanel\n"
        hit = b"H\t0\t0\t0\n"
        terminal = b"E\t2\t4\n"
        cases = (
            (b"invalid\n", 0, False),
            (prefix + hit, 0, False),
            (prefix + b"E\t2\t5\n", 0, False),
            (prefix + hit * 2 + terminal, 0, False),
            (prefix + b"H\t2\t0\t0\n" + terminal, 0, False),
            (prefix + b"H\t0\t4\t0\n" + terminal, 0, False),
            (prefix + b"H\t0\t0\t99999999999999999999999999\n" + terminal, 0, False),
            (prefix + b"H\t0\t2\t0\n" + terminal, 0, False),
            (prefix + hit + terminal + b"unexpected\n", 0, False),
            (prefix + hit + terminal, 7, False),
            (b"OVXPOS1\trecords\nH\t0\t0\t0\n", 0, True),
        )
        out, view = self.root / "bad.json", self.root / "bad.fa"
        for index, (response, status, originals_only) in enumerate(cases):
            with self.subTest(index=index):
                source = (
                    f"#!{sys.executable}\n"
                    "import sys\n"
                    "sys.stdin.buffer.read()\n"
                    f"if {originals_only!r} and sys.argv[-1] == 'panel':\n"
                    f"    sys.stdout.buffer.write({(prefix + hit + terminal)!r})\n"
                    "    raise SystemExit(0)\n"
                    f"sys.stdout.buffer.write({response!r})\n"
                    f"raise SystemExit({status})\n"
                ).encode()
                binary = self.write(f"bad-backend-{index}", source)
                binary.chmod(0o700)
                result = subprocess.run([sys.executable, "-B", "-I", str(SCRIPT), str(path),
                                         "--out", str(out), "--out-fasta", str(view),
                                         "--exact-placement-bin", str(binary)],
                                        capture_output=True, timeout=5)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(b"exact placement", result.stderr)
                self.assertFalse(out.exists())
                self.assertFalse(view.exists())

    @unittest.skipUnless(os.environ.get("OVASM_EXACT_PLACEMENT_BIN"), "requires built OVASM exact-placement binary")
    def test_native_batch_preserves_capped_orientation_and_record_order(self):
        binary = Path(os.environ["OVASM_EXACT_PLACEMENT_BIN"]).resolve(strict=True)
        rng = random.Random(742)
        words = [b"A", b"AA", b"AAA", b"ACGT", b"AACC", b"GGTT", b"TA", b"AT"]
        words += [bytes(rng.choices(b"ACGT", k=rng.randrange(3, 9))) for _ in range(31)]
        records = [b"GGTTNGAACCAAACGT", b"AACCNGGTTACGT", b"", b"N" * 65534 + b"ACGT"]
        records += [bytes(rng.choices(b"ACGTN", k=370)) for _ in range(170)]
        batch = TOOL["batch_placement_lookup"]
        panel = batch(binary, "panel", records, words)
        by_record = batch(binary, "records", records, words)
        for word in words:
            expected = []
            for index, seq in enumerate(records):
                single = [(index, *match) for match in TOOL["placements"](seq, word)]
                self.assertEqual(by_record(word, index), single)
                expected.extend(single)
            self.assertEqual(panel(word), expected[:2])
        self.assertEqual(panel(b"AACC"), [(0, 6, 10, "+"), (0, 0, 4, "-")])
        self.assertEqual(by_record(b"ACGT", 3), [(3, 65534, 65538, "+"), (3, 65534, 65538, "-")])

    @unittest.skipUnless(os.environ.get("OVASM_EXACT_PLACEMENT_BIN"), "requires built OVASM exact-placement binary")
    def test_native_batch_cli_outputs_match_both_python_policies(self):
        binary = Path(os.environ["OVASM_EXACT_PLACEMENT_BIN"]).resolve(strict=True)
        reverse = self.competing.translate(TOOL["COMPLEMENT"])[::-1]
        originals = [(b"target", self.target), (b"shared", self.target), (b"shared", reverse),
                     (b"inside", self.novel + self.target), (b"outside", self.target + self.novel),
                     (b"mixed", self.competing + self.target)]
        crops = [self.crop(name, seq) for name, seq in originals[:3]]
        crops += [self.crop(b"inside", self.target, start=64),
                  self.crop(b"outside", self.novel, start=64), self.crop(b"mixed", self.target, start=64),
                  self.crop(b"missing", self.novel)]
        path = self.assemble(originals, crops, as_gzip=True)
        bin_dir = self.root / "bin"
        bin_dir.mkdir()
        shutil.copy2(binary, bin_dir / ("ovasm.exe" if os.name == "nt" else "ovasm"))
        modes = (
            ("python", ["--python-placements"], self.root / "empty"),
            ("default", [], bin_dir),
            ("explicit", ["--exact-placement-bin", str(binary)], self.root / "empty"),
        )
        for policy in (TOOL["EXCLUDE_COMPETING"], TOOL["TARGET_ONLY"]):
            published = {}
            for label, extra, search_path in modes:
                out, view = [self.root / f"{policy}-{label}.{suffix}" for suffix in ("json", "fa")]
                subprocess.run([sys.executable, "-B", "-I", str(SCRIPT), str(path), "--retention-policy", policy,
                                "--out", str(out), "--out-fasta", str(view), *extra],
                               env={**os.environ, "PATH": str(search_path)}, check=True, timeout=5)
                published[label] = (out.read_bytes(), view.read_bytes())
            self.assertEqual(published["default"], published["python"])
            self.assertEqual(published["explicit"], published["python"])


if __name__ == "__main__":
    unittest.main()
