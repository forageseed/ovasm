"""Acceptance must detect chromosome errors hidden by aggregate coverage."""

import importlib.util
from pathlib import Path

import pytest

PATH = Path(__file__).resolve().parents[2] / "scripts/assembly_benchmarks/score_independent_multichrom.py"
SPEC = importlib.util.spec_from_file_location("independent_multichrom_score", PATH)
assert SPEC and SPEC.loader
scorer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(scorer)


def paf(q, qlen, qs, qe, t, tlen, ts, te, matches, block):
    return f"{q}\t{qlen}\t{qs}\t{qe}\t+\t{t}\t{tlen}\t{ts}\t{te}\t{matches}\t{block}\t60\n"


def test_circular_rotation_reverse_complement_and_panel_absent_recall():
    result = scorer.evaluate(
        {"circle": "ACGTTGCA", "line": "AAGTCC"},
        {"q1": "TGCAACGT", "q2": "GGACTT"}, {"circle"}, "", panel_absent={"line"},
    )
    assert result["accepted"]
    assert result["panel_absent_chromosomes"]["matched"] == 1
    assert not result["physical_topology_certified"]
    assert result["reference_record_scope"]["physical_chromosome_count"] is None
    assert result["reference_record_scope"]["matched_unit"] == "supplied reference assembly sequence record"


def test_linear_origin_cannot_be_rotated():
    result = scorer.evaluate({"line": "ACGTTGCC"}, {"q": "TGCCACGT"}, set(), "")
    assert not result["accepted"]


def test_merged_chromosomes_fail_despite_complete_aggregate_coverage():
    alignments = paf("merged", 20000, 0, 10000, "a", 10000, 0, 10000, 10000, 10000)
    alignments += paf("merged", 20000, 10000, 20000, "b", 10000, 0, 10000, 10000, 10000)
    result = scorer.evaluate({"a": "A" * 10000, "b": "C" * 10000}, {"merged": "A" * 10000 + "C" * 10000}, set(), alignments)
    assert all(c["aggregate_coverage"] == 1 for c in result["per_chromosome"].values())
    assert result["extra_sequence_bp"] == 0
    assert result["relative_length_difference"] == 0
    assert not result["one_to_one_chromosome_match"]
    assert not result["relaxed_accepted"]
    assert not result["accepted"]


def test_duplicate_one_chromosome_and_missing_other_fail_with_correct_length():
    result = scorer.evaluate({"a": "AAAACCCC", "b": "TTTTCCCC"}, {"q1": "AAAACCCC", "q2": "AAAACCCC"}, set(), "")
    assert result["relative_length_difference"] == 0
    assert result["matched_chromosomes"] == 1
    assert not result["accepted"]


def test_one_substitution_at_shifted_circular_origin_passes():
    line = paf("q", 10000, 0, 10000, "a", 20000, 7777, 17777, 9999, 10000)
    result = scorer.evaluate({"a": "A" * 10000}, {"q": "C" + "A" * 9999}, {"a"}, line)
    assert result["accepted"]
    assert result["per_chromosome"]["a"]["best_pair"]["target_coverage"] == 1


def test_low_identity_fails_even_when_length_and_coverage_pass():
    line = paf("q", 10000, 0, 10000, "a", 10000, 0, 10000, 9998, 10000)
    result = scorer.evaluate({"a": "A" * 10000}, {"q": "CC" + "A" * 9998}, set(), line)
    assert result["extra_sequence_bp"] == 0
    assert not result["accepted"]


def test_tandem_circular_traversal_does_not_become_a_bijective_pair():
    line = paf("q", 20000, 0, 20000, "a", 20000, 0, 20000, 20000, 20000)
    result = scorer.evaluate({"a": "A" * 10000}, {"q": "A" * 20000}, {"a"}, line)
    assert result["matched_chromosomes"] == 0
    assert not result["relaxed_accepted"]
    assert not result["accepted"]


def test_bijection_reassigns_ambiguous_identical_sequence():
    assert scorer.bijection({"q1": ["a", "b"], "q2": ["a"]}, ["q1", "q2"]) == {"q2": "a", "q1": "b"}


def test_empty_candidate_stays_a_failure_in_the_denominator():
    result = scorer.evaluate({"a": "ACGT"}, {}, set(), "")
    assert result["candidate_records"] == 0
    assert not result["accepted"]


@pytest.mark.parametrize("extra", ["GGGGTTTT", "AAAACCCC"])
def test_complete_reference_recovery_with_extra_records_has_two_scores(extra):
    result = scorer.evaluate(
        {"a": "AAAACCCC", "b": "ACGATGCA"},
        {"q1": "AAAACCCC", "q2": "ACGATGCA", "extra": extra}, set(), "",
    )
    assert result["relaxed_accepted"]
    assert result["scores"]["relaxed"]["score"] == 1
    assert not result["strict_accepted"]
    assert result["scores"]["strict"]["score"] == 0
    assert result["accepted"] is result["strict_accepted"]
    assert result["reference_record_match_fraction"] == 1
    assert len(result["unmatched_candidate_records"]) == 1
    assert "record_bijection" in result["scores"]["strict"]["failure_reasons"]
    assert result["candidate_bp"] == 24
    assert result["thresholds"]["length_tolerance"] == 0.001


@pytest.mark.parametrize("candidate,paf_text", [
    ({}, ""),
    ({"q1": "AAAACCCC", "q2": "AAAACCCC"}, ""),
    ({"merged": "AAAACCCCACGATGCA"}, ""),
])
def test_relaxed_score_does_not_accept_missing_duplicated_or_merged_references(candidate, paf_text):
    result = scorer.evaluate({"a": "AAAACCCC", "b": "ACGATGCA"}, candidate, set(), paf_text)
    assert not result["relaxed_accepted"]
    assert not result["strict_accepted"]
    assert result["scores"]["relaxed"]["score"] == 0


def test_relaxed_score_preserves_identity_and_single_record_coverage_requirements():
    alignment = paf("q", 10000, 0, 10000, "a", 10000, 0, 10000, 9998, 10000)
    result = scorer.evaluate({"a": "A" * 10000}, {"q": "CC" + "A" * 9998}, set(), alignment)
    assert not result["relaxed_accepted"]
    assert not result["strict_accepted"]


def test_both_protocols_pass_exact_complete_output():
    result = scorer.evaluate({"a": "ACGTACGA"}, {"q": "ACGTACGA"}, set(), "")
    assert result["accepted"] and result["relaxed_accepted"] and result["strict_accepted"]
    assert result["scores"]["strict"]["failure_reasons"] == []


@pytest.mark.parametrize("line", [
    paf("q", 5, 0, 4, "a", 4, 0, 4, 4, 4),
    paf("q", 4, 0, 4, "a", 4, 0, 5, 4, 4),
    paf("q", 4, 0, 4, "missing", 4, 0, 4, 4, 4),
])
def test_alignment_input_mismatch_fails_clearly(line):
    with pytest.raises(ValueError):
        scorer.evaluate({"a": "ACGT"}, {"q": "AAAA"}, set(), line)
