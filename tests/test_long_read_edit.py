"""Independent scalar and exact-oracle controls for the long-read edit search."""

from pathlib import Path
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from long_read_edit import (
    EditBudget, candidate_interval, cigar_edit_upper_bound, disjoint_blocks,
    edit_scores, reverse_complement, search_competitors,
)
from source_competition import parse_paf

MT, PT, NUC = "mitochondrion", "plastid", "nuclear"
LARGE = EditBudget(max_seed_hits=100_000, max_text_bases=100_000)


def scalar_scores(pattern, sequence, free_prefix):
    previous = list(range(len(pattern) + 1))
    for j, value in enumerate(sequence, 1):
        row = [0 if free_prefix else j]
        for i, base in enumerate(pattern, 1):
            row.append(min(row[-1] + 1, previous[i] + 1, previous[i - 1] + (base != value)))
        previous = row
        yield row[-1]


def full_substring_best(query, sequence, radius):
    best = None
    for strand in (query, reverse_complement(query)):
        for start in range(len(sequence)):
            text = sequence[start:start + len(query) + radius]
            for size, distance in enumerate(scalar_scores(strand, text, False), 1):
                if size >= len(query) - radius and distance <= radius and (best is None or distance < best):
                    best = distance
    return best


def make_paf(query_len, ref_len, ref_start, ref_stop, cigar, strand="+",
             query_start=0, query_stop=None):
    return parse_paf("\t".join(map(str, (
        "q", query_len, query_start, query_len if query_stop is None else query_stop, strand,
        "reference", ref_len, ref_start, ref_stop, 0, max(ref_len, query_len), 20, f"cg:Z:{cigar}",
    ))))


class EditPrimitiveTests(unittest.TestCase):
    def test_cigar_literal_substitution_insertion_deletion_and_reverse(self):
        ref = b"GGACGTCCGTTAA"
        query = b"ACGTACTTT"
        forward = make_paf(9, len(ref), 2, 11, "4=1I2M1D2=")
        reverse = make_paf(9, len(ref), 2, 11, "4M1I2X1D2M", strand="-")
        self.assertEqual(cigar_edit_upper_bound(forward, query.lower(), ref), 3)
        self.assertEqual(cigar_edit_upper_bound(reverse, reverse_complement(query), ref), 3)

    def test_skipped_partial_and_stale_length_cigar_rejected(self):
        ref = b"AACGTCC"
        for aln, raw, text in (
            (make_paf(6, 7, 0, 7, "3M1N3M"), b"AACGTC", ref),
            (make_paf(6, 7, 0, 5, "5M", query_start=1), b"AACGTC", ref),
            (make_paf(6, 7, 0, 6, "6M"), b"AACGTC", ref[:-1]),
        ):
            with self.assertRaises(ValueError):
                cigar_edit_upper_bound(aln, raw, text)

    def test_blocks_partition_and_candidate_bounds(self):
        self.assertEqual(disjoint_blocks(b"ACGTACG", 2), ((0, b"AC"), (2, b"GT"), (4, b"ACG")))
        self.assertEqual(candidate_interval(7, 2, 2, 5, 14), (1, 12))
        self.assertEqual(candidate_interval(7, 2, 0, 2, 8), (0, 8))
        self.assertIsNone(candidate_interval(7, 1, 0, 5, 7))
        with self.assertRaises(ValueError):
            disjoint_blocks(b"AN", 0)

    def test_bit_parallel_matches_independent_scalar_prefix_oracle(self):
        rng = random.Random(194)
        for _ in range(175):
            query = bytes(rng.choices(b"ACGTN", k=rng.randrange(1, 77)))
            text = bytes(rng.choices(b"ACGTN", k=rng.randrange(0, 81)))
            for free_prefix in (False, True):
                self.assertEqual(list(edit_scores(query, text, free_prefix=free_prefix)),
                                 list(scalar_scores(query, text, free_prefix)))


class CompleteSearchTests(unittest.TestCase):
    def test_search_matches_independent_all_substrings_oracle(self):
        rng = random.Random(434)
        for _ in range(325):
            query = bytes(rng.choices(b"ACGT", k=rng.randrange(1, 12)))
            radius = rng.randrange(min(len(query), 5))
            texts = [bytes(rng.choices(b"ACGTN", k=rng.randrange(1, 20))) for _ in range(3)]
            panel = {str(i): seq for i, seq in enumerate(texts)}
            roles = {"0": NUC, "1": NUC, "2": PT}
            result = search_competitors({"q": (query, radius)}, panel, roles, LARGE)["q"]
            for role, row in result.items():
                with self.subTest(query=query, radius=radius, role=role):
                    self.assertEqual(row.status, "complete")
                    expected = [distance for name, text in panel.items() if roles[name] == role
                                if (distance := full_substring_best(query, text, radius)) is not None]
                    hit = row.best_hit
                    self.assertEqual(hit.edit_distance if hit else None, min(expected, default=None))
                    if hit:
                        self.assertEqual(roles[hit.reference], role)
                        oriented = query if hit.strand == "+" else reverse_complement(query)
                        text = panel[hit.reference][slice(*hit.interval)]
                        self.assertEqual(list(scalar_scores(oriented, text, False))[-1], hit.edit_distance)

    def test_full_seed_validation_and_reverse_hit_coordinates(self):
        query = b"ACGTTGACCTTGATAAGCATCGTACGGGAAAGTCT"
        panel = {"prefix_only": query[:31] + b"NNNN",
                 "reverse_full": b"NN" + reverse_complement(query) + b"NN"}
        result = search_competitors({"q": (query.lower(), 0)}, panel,
                                    {"prefix_only": NUC, "reverse_full": PT}, LARGE)["q"]
        self.assertEqual(result[NUC].status, "complete")
        self.assertEqual(result[NUC].seed_hits, 0)
        self.assertIsNone(result[NUC].best_hit)
        hit = result[PT].best_hit
        self.assertEqual((hit.reference, hit.strand, hit.interval, hit.edit_distance),
                         ("reverse_full", "-", (2, len(query) + 2), 0))

    def test_seed_hit_cap_is_isolated_by_query_and_role(self):
        panel = {"repeated": b"AAAAA", "pt": b"CGT"}
        result = search_competitors({"repeat": (b"AA", 0), "single": (b"CGT", 0)},
                                    panel, {"repeated": NUC, "pt": PT}, EditBudget(1, 100))
        stopped = result["repeat"][NUC]
        self.assertEqual((stopped.status, stopped.reason, stopped.seed_hits), ("budget_exhausted", "seed_hits", 1))
        self.assertIsNone(stopped.best_hit)
        self.assertEqual(result["repeat"][PT].status, "complete")
        self.assertEqual(result["single"][PT].status, "complete")
        self.assertEqual(result["single"][PT].best_hit.edit_distance, 0)

    def test_exact_cap_accepts_all_and_zero_cap_can_prove_no_hits(self):
        panel = {"chromosome": b"AAANAA"}
        complete = search_competitors({"q": (b"AA", 0)}, panel, {"chromosome": NUC}, EditBudget(3, 100))["q"][NUC]
        self.assertEqual((complete.status, complete.seed_hits), ("complete", 3))
        absent = search_competitors({"q": (b"CG", 0)}, panel, {"chromosome": NUC}, EditBudget(0, 0))["q"][NUC]
        self.assertEqual((absent.status, absent.seed_hits), ("complete", 0))

    def test_panel_and_text_budgets_and_inapplicable_query(self):
        panel = {"nuclear": b"CCCGTCCAGGG", "pt": b"CGT"}
        roles = {"nuclear": NUC, "pt": PT}
        result = search_competitors({"q": (b"CGT", 0), "N": (b"AN", 0), "short": (b"A", 1)},
                                    panel, roles, EditBudget(100, 100, max_panel_bases=3))
        self.assertEqual((result["q"][NUC].status, result["q"][NUC].reason),
                         ("budget_exhausted", "panel_bases"))
        self.assertEqual(result["q"][PT].status, "complete")
        self.assertEqual(result["q"][PT].best_hit.edit_distance, 0)
        for key in ("N", "short"):
            self.assertTrue(all(row.status == "inapplicable" for row in result[key].values()))
        blocked = search_competitors({"q": (b"CGT", 0)}, panel, roles, EditBudget(100, 2))["q"]
        self.assertEqual((blocked[PT].status, blocked[PT].reason),
                         ("budget_exhausted", "edit_text_bases"))

    def test_reverse_start_recovery_is_bounded(self):
        query = b"ACGTCGAT"
        panel = {"nuclear": b"NNACGTAGATNN"}
        roles = {"nuclear": NUC}
        complete = search_competitors({"q": (query, 1)}, panel, roles, LARGE)["q"][NUC]
        self.assertEqual((complete.best_hit.interval, complete.best_hit.edit_distance), ((2, 10), 1))
        limited = search_competitors({"q": (query, 1)}, panel, roles,
                                      EditBudget(100, complete.examined_bases - 1))["q"][NUC]
        self.assertEqual((limited.status, limited.reason), ("budget_exhausted", "edit_traceback_bases"))
        self.assertIsNone(limited.best_hit)

    def test_invalid_roles_and_negative_budget_fail_closed(self):
        with self.assertRaises(ValueError):
            search_competitors({"q": (b"ACGT", 0)}, {"chr": b"ACT"}, {"other": NUC}, LARGE)
        with self.assertRaises(ValueError):
            EditBudget(-1, 20)


if __name__ == "__main__":
    unittest.main()
