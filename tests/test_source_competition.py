import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from source_competition import (
    Anchor, Placement, QueryScope, RoleSearch, decide, parse_paf, verify_witness,
)

MT, PT, NUC = "mitochondrion", "plastid", "nuclear"
REF = b"AACCGTGGCATGAATCCCTG"


def reverse_complement(raw):
    return raw.translate(bytes.maketrans(b"ACGT", b"TGCA"))[::-1]


def paf(cigar, *, ref="mt", strand="+", qlen=20, qa=0, qb=20,
        tlen=20, ta=0, tb=20, matches=None, columns=None, extra=()):
    if matches is None:
        matches = min(qb - qa, tb - ta)
    if columns is None:
        columns = max(qb - qa, tb - ta)
    return "\t".join(map(str, ("original", qlen, qa, qb, strand, ref, tlen, ta, tb,
                              matches, columns, 0, "cg:Z:" + cigar, *extra)))


def witness(aln, raw, ref=REF, role=MT, boundary=Anchor((8, 12), (8, 12)),
            left=Anchor((4, 8), (4, 8)), right=Anchor((12, 16), (12, 16)),
            exact=True, span=None):
    scope = QueryScope.bind("original", raw, (0, len(raw)) if span is None else span)
    return verify_witness(scope, aln, raw, ref, role, boundary, left, right, exact_boundary=exact)


class GeometryTests(unittest.TestCase):
    def test_exact_forward_joins_adjacent_diagonal_cigar_runs(self):
        aln = parse_paf(paf("6M4=10M"))
        self.assertEqual(aln.project_exact((4, 12), REF.lower(), REF), (4, 12))
        self.assertIsNotNone(witness(aln, REF))

    def test_exact_reverse_unclipped_and_clipped(self):
        # Reference fragment [2,18) occurs in the reverse-complement of raw [3,19).
        raw = b"AAA" + reverse_complement(REF[2:18]) + b"C"
        aln = parse_paf(paf("16M", strand="-", qa=3, qb=19, ta=2, tb=18))
        self.assertEqual(aln.project_exact((13, 17), raw, REF), (4, 8))
        self.assertIsNone(aln.project_exact((2, 17), raw, REF))
        kwargs = dict(boundary=Anchor((9, 13), (8, 12)), left=Anchor((13, 17), (4, 8)),
                      right=Anchor((5, 9), (12, 16)))
        self.assertIsNone(witness(aln, raw, **kwargs))
        hit = witness(aln, raw, **kwargs, span=(3, 19))
        self.assertEqual(hit.placement, Placement(MT, "mt", (8, 12), "-"))
        self.assertEqual(hit.left.raw, (13, 17))

    def test_literal_projection_rejects_false_match_and_bad_anchor_offset(self):
        aln = parse_paf(paf("20M"))
        mismatch = REF[:6] + b"N" + REF[7:]
        self.assertIsNone(aln.project_exact((4, 8), mismatch, REF))
        self.assertIsNone(witness(aln, mismatch))
        self.assertIsNone(witness(aln, REF, left=Anchor((4, 8), (3, 7))))
        aln = parse_paf(paf("6M1X13M", matches=19))
        self.assertIsNone(aln.project_exact((4, 8), REF, REF))

    def test_gapped_deletion_between_exact_outer_anchors(self):
        raw = REF[:10] + REF[11:]
        aln = parse_paf(paf("10M1D9M", qlen=19, qb=19))
        kw = dict(boundary=Anchor((8, 11), (8, 12)), right=Anchor((11, 15), (12, 16)))
        self.assertIsNotNone(witness(aln, raw, **kw, exact=False))
        self.assertIsNone(witness(aln, raw, **kw, exact=True))
        self.assertIsNone(witness(aln, raw, **{**kw, "boundary": Anchor((8, 11), (9, 12))},
                                  exact=False))
        self.assertIsNone(witness(aln, raw, **{**kw, "left": Anchor((7, 11), (7, 12))},
                                  exact=False))
        self.assertIsNone(aln.project_exact((9, 11), raw, REF))

    def test_gapped_insertion_between_exact_outer_anchors(self):
        raw = REF[:10] + b"A" + REF[10:]
        aln = parse_paf(paf("10M1I10M", qlen=21, qb=21))
        self.assertIsNotNone(witness(aln, raw, exact=False,
                                     boundary=Anchor((8, 13), (8, 12)),
                                     right=Anchor((13, 17), (12, 16))))
        self.assertIsNone(aln.project_exact((9, 12), raw, REF))
        self.assertFalse(aln.projects_boundary(Anchor((10, 11), (10, 11))))

    def test_skip_and_disjoint_alignment_never_complete_a_bridge(self):
        raw = REF[:10] + REF[11:]
        skip = parse_paf(paf("10M1N9M", qlen=19, qb=19))
        self.assertIsNone(witness(skip, raw, exact=False,
                                  boundary=Anchor((8, 11), (8, 12)),
                                  right=Anchor((11, 15), (12, 16))))
        left = parse_paf(paf("10M", qb=10, tb=10))
        right = parse_paf(paf("10M", qa=10, ta=10))
        self.assertIsNone(witness(left, REF))
        self.assertIsNone(witness(right, REF))

    def test_duplicate_or_malformed_cigar_and_unbound_read(self):
        for line in (
            paf("9M10M"), paf("20M0I"), paf("6Mgarbage14M"), paf("10M10S"),
            paf("20M", strand="."), paf("20M", ta=-1), paf("20M", extra=("cg:Z:20M",)),
            paf("20M").rsplit("\t", 1)[0],
        ):
            with self.subTest(line=line), self.assertRaises(ValueError):
                parse_paf(line)
        aln = parse_paf(paf("20M"))
        scope = QueryScope.bind("different_original", REF, (0, 20))
        with self.assertRaises(ValueError):
            verify_witness(scope, aln, REF, REF, MT, Anchor((8, 12), (8, 12)),
                           Anchor((4, 8), (4, 8)), Anchor((12, 16), (12, 16)), exact_boundary=True)


class DecisionTests(unittest.TestCase):
    def setUp(self):
        self.w = witness(parse_paf(paf("20M")), REF)
        self.scope = self.w.scope
        self.placement = self.w.placement
        self.nuclear = Placement(NUC, "chr2", (8, 12), "+")

    def searches(self, mt=(), pt=(), nuc=(), incomplete=None, scope=None):
        scope = self.scope if scope is None else scope
        return tuple(RoleSearch(scope, role, role != incomplete, items)
                     for role, items in ((MT, mt), (PT, pt), (NUC, nuc)))

    def test_single_placement_and_exact_duplicate_evidence(self):
        searches = self.searches(mt=(self.placement, self.placement))
        got = decide(self.scope, MT, searches, (self.w, self.w))
        self.assertEqual(got.state, "target_supported_local")
        self.assertEqual(got.placement, self.placement)
        self.assertEqual(got, decide(self.scope, MT, searches[::-1], (self.w,)))

    def test_nuclear_and_plastid_are_distinct_competitor_roles(self):
        n = witness(parse_paf(paf("20M", ref="chr2")), REF, role=NUC)
        p = witness(parse_paf(paf("20M", ref="pt")), REF, role=PT)
        got = decide(self.scope, MT, self.searches(nuc=(n.placement,)), (n,))
        self.assertEqual(got.state, "competing_supported_local")
        self.assertEqual(got.placement.role, NUC)
        got = decide(self.scope, MT, self.searches(pt=(p.placement,)), (p,))
        self.assertEqual(got.state, "competing_supported_local")
        self.assertEqual(got.placement.role, PT)
        got = decide(self.scope, PT, self.searches(pt=(p.placement,)), (p,))
        self.assertEqual(got.state, "target_supported_local")

    def test_complete_search_with_credible_competitor_stays_unresolved(self):
        got = decide(self.scope, MT, self.searches(mt=(self.placement,), nuc=(self.nuclear,)), (self.w,))
        self.assertEqual((got.state, got.reason), ("unresolved", "competing_roles"))

    def test_unwitnessed_or_incomplete_stays_unresolved(self):
        for inputs, witnesses, reason in (
            (self.searches(mt=(self.placement,), incomplete=NUC), (self.w,), "incomplete_competition"),
            (self.searches(mt=(self.placement,))[:2], (self.w,), "incomplete_competition"),
            (self.searches(mt=(self.placement,)), (), "no_verified_geometry"),
            (self.searches(), (self.w,), "witness_outside_search"),
            (self.searches(mt=(self.placement,), scope=QueryScope.bind("original", REF, (1, 20))),
             (self.w,), "competition_scope_mismatch"),
            (self.searches(mt=(self.placement,), scope=QueryScope.bind("original", b"T" + REF[1:], (0, 20))),
             (self.w,), "competition_scope_mismatch"),
        ):
            with self.subTest(reason=reason):
                got = decide(self.scope, MT, inputs, witnesses)
                self.assertEqual((got.state, got.reason), ("unresolved", reason))

    def test_complete_requires_a_real_boolean_and_local_source_scope(self):
        searches = (RoleSearch(self.scope, MT, True, (self.placement,)),
                    RoleSearch(self.scope, PT, True, ()),
                    RoleSearch(self.scope, NUC, "unknown", ()))
        self.assertEqual(decide(self.scope, MT, searches, (self.w,)).reason,
                         "incomplete_competition")
        other = witness(parse_paf(paf("20M")), REF, span=(2, 18))
        self.assertEqual(decide(self.scope, MT, self.searches(mt=(self.placement,)),
                                (other,)).reason, "witness_scope_mismatch")

    def test_known_alias_is_deduplicated_but_unverified_alias_is_ambiguous(self):
        alias = Placement(MT, "alternate", (1, 5), "-")
        searches = self.searches(mt=(alias, self.placement))
        self.assertEqual(decide(self.scope, MT, searches, (self.w,)).reason, "multiple_placements")
        got = decide(self.scope, MT, searches, (self.w,), aliases={alias: self.placement})
        self.assertEqual(got.state, "target_supported_local")
        with self.assertRaises(ValueError):
            decide(self.scope, MT, searches, (self.w,), aliases={alias: self.nuclear})


if __name__ == "__main__":
    unittest.main()
