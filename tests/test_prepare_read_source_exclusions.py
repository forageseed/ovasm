import gzip
from hashlib import sha256
import json
from pathlib import Path
import random
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/prepare_read_source_exclusions.py"
TARGET = "mitochondrion"


def reverse(sequence):
    return sequence.translate(bytes.maketrans(b"ACGT", b"TGCA"))[::-1]


def record(name, sequence, quality=None):
    quality = b"I" * len(sequence) if quality is None else quality
    return b"@" + name + b"\n" + sequence + b"\n+\n" + quality + b"\n"


class SourceExclusionTest(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        randomizer = random.Random(9406)
        self.outside, self.inside, self.target, self.background = (
            bytes(randomizer.choices(b"ACGT", k=40)) for _ in range(4)
        )
        self.window = self.outside + self.inside
        proof = self.file("registered.txt", b"independently registered boundary\n")
        self.plan = {
            "schema": "ovasm.read-source-exclusion-plan.v1",
            "scope": "synthetic source control",
            "target": TARGET,
            "references": [
                self.reference("nuc", "nuclear", self.window),
                self.reference("mt", TARGET, self.target),
            ],
            "raw_reads": [],
            "windows": [
                dict(id="outer", reference="nuc", record="nuc", interval=[0, 80],
                     boundary=40, exterior_side="left", provenance=proof),
                dict(id="protected", reference="mt", record="mt", interval=[0, 40],
                     boundary=20, provenance=proof),
            ],
        }

    def file(self, name, data):
        (self.root / name).write_bytes(data)
        return {"path": name, "sha256": sha256(data).hexdigest()}

    def reference(self, name, role, sequence):
        return dict(self.file(name + ".fa", b">" + name.encode() + b"\n" + sequence + b"\n"),
                    id=name, role=role)

    def invoke(self):
        plan = self.root / "plan.json"
        plan.write_text(json.dumps(self.plan))
        return subprocess.run(
            [sys.executable, "-B", str(SCRIPT), str(plan), "--out", str(self.root / "audit.json"),
             "--out-exclusions", str(self.root / "exclusions.json")],
            capture_output=True, text=True,
        )

    def result(self, name):
        return json.loads((self.root / name).read_text())

    def test_original_identity_and_strand_with_abstentions(self):
        name = b"forward full header"
        forward = self.background + self.window.lower()
        quality = b"#" * len(forward)
        backward = self.background + reverse(self.window)
        self.plan["raw_reads"] = [
            self.file("a.fq", record(name, forward, quality) + record(b"duplicate", self.background)
                      + record(b"same-first-word another name", self.target)),
            self.file("b.fq.gz", gzip.compress(
                record(b"reverse", backward)
                + record(b"mixed", self.window + self.target)
                + record(b"repeated", self.window * 2)
                + record(b"duplicate", self.window)
            )),
        ]
        run = self.invoke()
        self.assertEqual(run.returncode, 0, run.stderr)
        audit = self.result("audit.json")
        manifest = self.result("exclusions.json")
        self.assertEqual(audit["counts"], {
            "no_registered_window": 1, "target_context_only": 1, "competing_context": 2,
            "mixed_registered_roles": 1, "ambiguous_or_shared_competing_context": 1,
            "duplicate_original_header": 1,
        })
        excluded = manifest["records"]
        self.assertEqual([row["read_id_hex"] for row in excluded], [name.hex(), b"reverse".hex()])
        self.assertEqual(excluded[0]["sequence_sha256"], sha256(forward).hexdigest())
        self.assertEqual(excluded[0]["quality_sha256"], sha256(quality).hexdigest())
        for row in excluded:
            self.assertEqual(row["target"], TARGET)
            self.assertEqual(row["provenance"], {
                "path": str(self.root / "audit.json"),
                "sha256": sha256((self.root / "audit.json").read_bytes()).hexdigest(),
            })
        forward_record, reverse_record = [r for r in audit["records"] if r["state"] == "competing_context"]
        self.assertEqual([forward_record["original_record_ordinal0"], reverse_record["original_record_ordinal0"]], [0, 0])
        self.assertEqual([forward_record["original_input_index0"], reverse_record["original_input_index0"]], [0, 1])
        for item, strand, coords in [(forward_record, "+", [40, 80]), (reverse_record, "-", [80, 120])]:
            observation = item["observations"]["outer"]
            self.assertEqual(observation["original_positions"][0]["strand"], strand)
            self.assertEqual(observation["original_exterior_intervals"], [coords])

    def test_shared_exterior_blocks_a_unique_full_window(self):
        self.plan["references"].append(self.reference("pt", "plastid", reverse(self.outside)))
        self.plan["raw_reads"] = [self.file("reads.fa", b">read\n" + self.window + b"\n")]
        run = self.invoke()
        self.assertEqual(run.returncode, 0, run.stderr)
        audit = self.result("audit.json")
        self.assertEqual(len(audit["windows"]["outer"]["reference_placements"]), 1)
        self.assertEqual(audit["windows"]["outer"]["exterior_other_roles"], ["plastid"])
        self.assertIsNone(audit["records"][0]["quality_sha256"])
        self.assertEqual(audit["eligible_records"], 0)
        self.assertFalse((self.root / "exclusions.json").exists())

    def test_same_role_duplicate_reference_placement_abstains(self):
        self.plan["references"].append(self.reference("alias", "nuclear", self.window))
        self.plan["raw_reads"] = [self.file("reads.fq", record(b"read", self.window))]
        run = self.invoke()
        self.assertEqual(run.returncode, 0, run.stderr)
        audit = self.result("audit.json")
        self.assertEqual(len(audit["windows"]["outer"]["reference_placements"]), 2)
        self.assertEqual(audit["eligible_records"], 0)
        self.assertFalse((self.root / "exclusions.json").exists())

    def test_original_fasta_keeps_non_newline_bytes(self):
        header = b"forward full header  "
        prefix = b"\t  "
        line_break = b"\r\n"
        bases = self.window.lower()
        self.plan["raw_reads"] = [self.file("reads.fa", (
            b">" + header + line_break + prefix + bases[:25] + line_break + bases[25:] + line_break
            + b">noncontiguous" + line_break + bases[:40] + b" \t" + line_break + bases[40:]
        ))]
        run = self.invoke()
        self.assertEqual(run.returncode, 0, run.stderr)
        audit = self.result("audit.json")
        self.assertEqual(audit["counts"], {"competing_context": 1, "no_registered_window": 1})
        row, = self.result("exclusions.json")["records"]
        self.assertEqual(row["read_id_hex"], header.hex())
        self.assertEqual(row["sequence_sha256"], sha256(prefix + bases).hexdigest())
        self.assertIsNone(row["quality_sha256"])
        observation = audit["records"][0]["observations"]["outer"]
        self.assertEqual(observation["original_exterior_intervals"], [[3, 43]])

    def test_changed_provenance_aborts_without_publishing(self):
        self.plan["raw_reads"] = [self.file("reads.fq", record(b"read", self.window))]
        (self.root / "registered.txt").write_text("changed")
        run = self.invoke()
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("input sha-256 mismatch", run.stderr.lower())
        self.assertFalse((self.root / "audit.json").exists())
        self.assertFalse((self.root / "exclusions.json").exists())


if __name__ == "__main__":
    unittest.main()
