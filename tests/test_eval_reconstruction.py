#!/usr/bin/env python3
"""Behavioural controls for scripts/eval_reconstruction.py (PAF-level; minimap2 not needed)."""
import importlib.util
from pathlib import Path
import unittest

SOURCE = Path(__file__).resolve().parents[1] / "scripts/eval_reconstruction.py"
spec = importlib.util.spec_from_file_location("eval_reconstruction", SOURCE)
ev = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ev)


def paf(q, qs, qe, t, tl, ts, te, match, cols, primary=True):
    return "\t".join(map(str, [q, 0, qs, qe, "+", t, tl, ts, te, match, cols, 60, f"tp:A:{'P' if primary else 'S'}"]))


class EvalTests(unittest.TestCase):
    def test_wrap_around_alignment_on_a_doubled_circle_counts_once(self):
        # truth circle 10 kb doubled to 20 kb; the query starts at 7 kb of the circle
        rows = ev.parse_paf([paf("q", 0, 10000, "t", 20000, 7000, 17000, 10000, 10000)])
        s = ev.score_v2(rows, {"t": 10000}, {"q": 10000})
        self.assertAlmostEqual(s["v2_recall"], 1.0)
        self.assertAlmostEqual(s["v2_precision"], 1.0)
        self.assertEqual(s["v2_copy_deficit_bp"], 0)

    def test_low_identity_and_short_blocks_are_not_good(self):
        rows = ev.parse_paf([paf("q", 0, 5000, "t", 10000, 0, 5000, 4800, 5000),
                             paf("q", 5000, 5500, "t", 10000, 5000, 5500, 500, 500)])
        s = ev.score_v2(rows, {"t": 10000}, {"q": 5500})
        self.assertEqual(s["v2_recall"], 0.0)
        self.assertAlmostEqual(s["v2_identity"], 0.96)

    def test_a_repeat_held_once_is_a_copy_deficit(self):
        # truth has a 3 kb repeat twice; the query holds it once: primary to copy 1, secondary to copy 2
        rows = ev.parse_paf([paf("q", 0, 3000, "t", 20000, 1000, 4000, 3000, 3000),
                             paf("q", 0, 3000, "t", 20000, 6000, 9000, 3000, 3000, primary=False)])
        s = ev.score_v2(rows, {"t": 10000}, {"q": 3000})
        self.assertEqual(s["v2_copy_deficit_bp"], 3000)
        self.assertAlmostEqual(s["v2_recall"], 0.3)

    def test_secondaries_do_not_inflate_v2_but_v1_takes_every_alignment(self):
        rows = ev.parse_paf([paf("q", 0, 2000, "t", 4000, 0, 2000, 2000, 2000),
                             paf("q", 0, 2000, "t", 4000, 2000, 4000, 1990, 2000, primary=False)])
        self.assertAlmostEqual(ev.score_v1(rows, {"t": 4000}, {"q": 2000})[0], 1.0)

    def test_read_supported_bases_absent_from_truth(self):
        rows = ev.parse_paf([paf("q", 0, 4000, "t", 8000, 0, 4000, 4000, 4000)])
        depth = {"q": [30] * 4000 + [25] * 1000 + [0] * 1000}
        s = ev.score_v2(rows, {"t": 4000}, {"q": 6000}, depth)
        self.assertEqual(s["v2_read_supported_absent_bp"], 1000)
        self.assertAlmostEqual(s["v2_precision_supported"], 5000 / 6000)
        self.assertEqual(s["v2_query_unaligned_bp"], 2000)

    def test_fold(self):
        self.assertEqual(ev.fold(9000, 12000, 10000), [(9000, 10000), (0, 2000)])
        self.assertEqual(ev.fold(0, 25000, 10000), [(0, 10000)])


if __name__ == "__main__":
    unittest.main()
