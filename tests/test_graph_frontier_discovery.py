import gzip
import hashlib
import io
import json
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "discover_graph_frontier.py"
TOOL = runpy.run_path(str(SCRIPT))


class FrontierGraphWords(unittest.TestCase):
    def words(self, gfa, k, budget=100_000):
        sequences, links = TOOL["parse_graph"](io.BytesIO(gfa))
        index, states, complete = TOOL["graph_words"](sequences, links, k, budget)
        return sequences, links, index, states, complete

    def test_multiple_exact_overlaps_and_crossing_disqualifies_unique(self):
        graph = (b"S\ta\tAC\nS\tb\tCG\nS\tc\tGT\nS\td\tTA\nS\tduplicate\tACGTA\nS\tunique\tTCCAGC\n"
                 b"L\ta\t+\tb\t+\t1M\nL\tb\t+\tc\t+\t1M\nL\tc\t+\td\t+\t1M\nL\td\t-\tc\t-\t1M\n")
        sequences, links, words, _, complete = self.words(graph, 5)
        self.assertTrue(complete)
        self.assertEqual(links["d-"], [("c-", 1)])
        self.assertIn(b"ACGTA", words)
        self.assertIsNone(words[b"ACGTA"])
        self.assertIsNone(words[b"TACGT"])
        self.assertEqual(words[b"TCCAG"], ("unique+", 0))
        gaps, represented, unique, observable, unobservable = TOOL["events"](b"ACGTA", words, sequences, links, 5)
        self.assertEqual((gaps, represented, unique, observable, unobservable), ([], 1, 0, 1, 0))

    def test_fully_overlapped_short_cycle_has_finite_walk(self):
        graph = (b"S\ta\tTGA\nS\tb\tA\nS\tc\tAC\nS\td\tCG\n"
                 b"L\ta\t+\tb\t+\t1M\nL\tb\t+\tb\t+\t1M\nL\tb\t+\tc\t+\t1M\nL\tc\t+\td\t+\t1M\n")
        _, _, words, states, complete = self.words(graph, 5)
        self.assertTrue(complete)
        self.assertLess(states, 100)
        self.assertIn(b"TGACG", words)
        self.assertIn(b"CGTCA", words)

    def test_reverse_raw_and_forward_graph_offsets(self):
        graph = b"S\tu\tAAAGTCC\n"
        sequences, links, words, _, _ = self.words(graph, 3)
        gaps, *_ = TOOL["events"](b"AAATTGTCC", words, sequences, links, 3)
        self.assertEqual(len(gaps), 1)
        self.assertEqual(gaps[0]["absent_word_start_interval0"], [1, 5])
        self.assertEqual(gaps[0]["affected_raw_interval0"], [1, 7])
        self.assertTrue(gaps[0]["two_unique_noncollinear_oriented_anchors"])
        self.assertFalse(gaps[0]["terminal_or_branch_port"])
        reverse, *_ = TOOL["events"](b"GGACAATTT", words, sequences, links, 3)
        self.assertEqual(reverse[0]["left_anchor"]["orientation"], "-")
        self.assertEqual(reverse[0]["left_anchor"]["oriented_interval0"], [1, 4])
        self.assertEqual(reverse[0]["left_anchor"]["forward_segment_interval0"], [3, 6])
        end, *_ = TOOL["events"](b"AAAGTCCAAC", words, sequences, links, 3)
        self.assertTrue(end[0]["terminal_or_branch_port"])

    def test_ambiguous_raw_words_break_runs(self):
        sequences, links, words, _, _ = self.words(b"S\tu\tAAAGTCC\n", 3)
        for raw in (b"AAANAGT", b"aaaRagt"):
            with self.subTest(raw=raw):
                gaps, represented, unique, observable, unobservable = TOOL["events"](
                    raw, words, sequences, links, 3)
                self.assertEqual((gaps, represented, unique, observable, unobservable), ([], 2, 2, 2, 3))

    def test_exact_overlap_geometry_and_unsupported_records_fail(self):
        for graph in (b"S\ta\tAC\nS\tb\tGT\nL\ta\t+\tb\t+\t1M\n",
                      b"S\ta\tAC\nS\tb\tCG\nL\ta\t+\tb\t+\t1I\n",
                      b"S\ta\tAC\nJ\ta\t+\ta\t-\t1\n"):
            with self.subTest(graph=graph), self.assertRaises(ValueError):
                self.words(graph, 3)


class FrontierDiscoveryCli(unittest.TestCase):
    def setUp(self):
        work = tempfile.TemporaryDirectory()
        self.addCleanup(work.cleanup)
        self.root = Path(work.name)
        self.manifest = self.root / "manifest.json"
        self.report = self.root / "report.json"
        self.evidence = self.root / "rows.jsonl.gz"
        self.data = dict(schema=TOOL["MANIFEST"], scope="synthetic fixture",
                         graph=self.ref("graph.gfa", b"S\ta\tAAAGTCC\nS\tb\tCAAC\nL\ta\t+\tb\t+\t1M\n"),
                         raw_reads=[self.ref("records.fastq.gz", gzip.compress(
                             b"@found full header\xff\naaaTTgTCC\n+\nIIIIIIIII\n"
                             b"@no-word\nCCCCCCCCC\n+\nIIIIIIIII\n"
                             b"@unobservable\nNNNNNNNNN\n+\nIIIIIIIII\n"))],
                         recruit=self.ref("recruit.json", b'{"tool":"ovasm find recruit","k":3}'))

    def ref(self, name, raw):
        path = self.root / name
        path.write_bytes(raw)
        return {"path": name, "sha256": hashlib.sha256(raw).hexdigest()}

    def run_cli(self, *flags):
        self.manifest.write_text(json.dumps(self.data))
        return subprocess.run([sys.executable, "-B", "-I", str(SCRIPT), str(self.manifest),
                               "--out", str(self.report), "--out-evidence", str(self.evidence), *flags],
                              text=True, capture_output=True)

    def test_hashes_output_and_no_graph_word_opportunities(self):
        result = self.run_cli("--include-no-graph-words")
        self.assertEqual(result.returncode, 0, result.stderr)
        d = json.loads(self.report.read_text())
        rows = [json.loads(s) for s in gzip.decompress(self.evidence.read_bytes()).splitlines()]
        self.assertEqual(d["schema"], TOOL["REPORT"])
        self.assertEqual(d["scan_k"], 3)
        self.assertEqual(d["graph_index"]["status"], "complete")
        self.assertEqual(d["raw_reads"][0]["counts"]["no_graph_words"], 1)
        self.assertEqual(d["evidence"]["sha256"], hashlib.sha256(self.evidence.read_bytes()).hexdigest())
        self.assertEqual(rows[0]["read_id_hex"], b"found full header\xff".hex())
        self.assertEqual(rows[0]["quality_sha256"], hashlib.sha256(b"IIIIIIIII").hexdigest())
        self.assertEqual(rows[0]["sequence_sha256"], hashlib.sha256(b"aaaTTgTCC").hexdigest())
        self.assertEqual(rows[1]["status"], "no_graph_words")
        self.assertEqual(rows[1]["events"], [])
        self.assertEqual(rows[2]["status"], "no_observable_words")
        self.assertEqual(rows[2]["unobservable_word_positions"], 7)
        self.assertEqual(d["raw_reads"][0]["counts"]["no_observable_words"], 1)

    def test_incomplete_index_cannot_assert_absence_or_uniqueness(self):
        result = self.run_cli("--max-traversal-states", "1", "--include-no-graph-words")
        self.assertEqual(result.returncode, 0, result.stderr)
        d = json.loads(self.report.read_text())
        self.assertEqual(d["graph_index"]["status"], "truncated")
        self.assertEqual(d["graph_index"]["traversal_states"], 1)
        self.assertIsNone(d["graph_index"]["unique_oriented_node_words"])
        self.assertEqual(d["raw_scan_status"], "not_run_incomplete_graph_index")
        self.assertIsNone(d["raw_reads"][0]["counts"])
        self.assertEqual(gzip.decompress(self.evidence.read_bytes()), b"")

    def test_manifest_rejects_stale_source_and_output_alias(self):
        self.data["raw_reads"][0]["sha256"] = "00" * 32
        failed = self.run_cli()
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("SHA-256 mismatch", failed.stderr)
        self.assertFalse(self.evidence.exists())
        self.assertFalse(self.report.exists())
        self.report.write_text("existing\n")
        failed = self.run_cli()
        self.assertNotEqual(failed.returncode, 0)
        self.assertEqual(self.report.read_text(), "existing\n")
        self.assertFalse(self.evidence.exists())

    def test_explicit_k_or_recruit_and_distinct_input_files(self):
        self.data["scan_k"] = 3
        bad = self.run_cli()
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn("exactly one", bad.stderr)
        del self.data["recruit"]
        self.data["raw_reads"].append(self.data["raw_reads"][0])
        bad = self.run_cli()
        self.assertNotEqual(bad.returncode, 0)
        self.assertIn("duplicate raw input", bad.stderr)
        self.assertFalse(self.evidence.exists())

    def test_fasta_shorter_than_k_is_counted(self):
        del self.data["recruit"]
        self.data["scan_k"] = 3
        self.data["raw_reads"] = [self.ref("records.fa", b">short full ID\nAC\n>case\nAaa\nGtc\n")]
        result = self.run_cli("--include-no-graph-words")
        self.assertEqual(result.returncode, 0, result.stderr)
        stats = json.loads(self.report.read_text())["raw_reads"][0]["counts"]
        self.assertEqual(stats["shorter_than_scan_k"], 1)
        self.assertEqual(stats["any_graph_words"], 1)
        self.assertEqual(stats["no_graph_words"], 0)
        self.assertEqual(stats["no_observable_words"], 1)


if __name__ == "__main__":
    unittest.main()
