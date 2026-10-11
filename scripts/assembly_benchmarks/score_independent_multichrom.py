"""Post-assembly reference-record acceptance; no references in assembly selection.

Circular reference sequences are doubled only for alignment across their arbitrary
origin. A passing pair must cover both sequences in one alignment, traversing the
reference at most once. Aggregate coverage cannot substitute for a bijection.
Sequence acceptance does not establish a physical molecule's shape or count.
Both protocols use the same full candidate FASTA and per-reference pair criteria.
Relaxed acceptance requires every reference to have a distinct passing candidate;
strict acceptance also requires a bijection, total-length and extra-sequence limits.
The historical ``accepted`` field retains strict acceptance semantics.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in records:
                raise ValueError(f"Duplicate FASTA identifier: {name}")
            records[name] = []
        elif line.strip():
            if name is None:
                raise ValueError("Sequence before a FASTA header")
            records[name].append(line.strip().upper())
    sequences = {key: "".join(value) for key, value in records.items()}
    if any(not value for value in sequences.values()):
        raise ValueError("Empty FASTA record")
    return sequences


def union_length(intervals: list[tuple[int, int]]) -> int:
    total, end = 0, 0
    for left, right in sorted(intervals):
        if right > end:
            total += right - max(left, end)
            end = right
    return total


def circular_intervals(start: int, end: int, length: int) -> list[tuple[int, int]]:
    span = end - start
    if span >= length:
        return [(0, length)]
    left = start % length
    right = left + span
    return [(left, right)] if right <= length else [(left, length), (0, right - length)]


def bijection(edges: dict[str, list[str]], queries: list[str]) -> dict[str, str]:
    owners: dict[str, str] = {}

    def visit(query: str, seen: set[str]) -> bool:
        for target in edges.get(query, []):
            if target in seen:
                continue
            seen.add(target)
            if target not in owners or visit(owners[target], seen):
                owners[target] = query
                return True
        return False

    for query in queries:
        visit(query, set())
    return {query: target for target, query in owners.items()}


def evaluate(
    truth: dict[str, str], candidate: dict[str, str], circular: set[str],
    paf: str, *, identity_min: float = 0.9999, coverage_min: float = 0.999,
    length_tolerance: float = 0.001, extra_fraction_max: float = 0.001,
    panel_absent: set[str] | None = None,
) -> dict:
    if not truth or not circular <= truth.keys():
        raise ValueError("Truth must be nonempty and circular identifiers must exist")
    if panel_absent is not None and not panel_absent <= truth.keys():
        raise ValueError("Panel-absent identifiers must exist in truth")
    pairs: dict[tuple[str, str], dict] = {}
    query_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    target_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    complement = str.maketrans("ACGTN", "TGCAN")
    for query, sequence in candidate.items():
        reverse = sequence.translate(complement)[::-1]
        for target, reference in truth.items():
            if len(sequence) != len(reference):
                continue
            searchable = reference + reference if target in circular else reference
            if sequence not in searchable and reverse not in searchable:
                continue
            if target not in circular and sequence != reference and reverse != reference:
                continue
            pairs[query, target] = {
                "query": query, "target": target, "identity": 1.0,
                "query_coverage": 1.0, "target_coverage": 1.0,
                "reference_span": len(reference), "method": "exact_sequence",
            }
            query_intervals[query].append((0, len(sequence)))
            target_intervals[target].append((0, len(reference)))
    for line in paf.splitlines():
        row = line.split("\t")
        query, target = row[0], row[5]
        if query not in candidate or target not in truth:
            raise ValueError("PAF identifiers differ from the scoring inputs")
        qlen, qs, qe, tlen, ts, te, matches, block = map(
            int, (row[1], row[2], row[3], row[6], row[7], row[8], row[9], row[10])
        )
        length = len(truth[target])
        if qlen != len(candidate[query]) or tlen != length * (2 if target in circular else 1):
            raise ValueError("PAF lengths differ from the scoring inputs")
        if not (0 <= qs <= qe <= qlen and 0 <= ts <= te <= tlen and 0 <= matches <= block):
            raise ValueError("Invalid PAF coordinates or match count")
        identity = matches / block if block else 0.0
        intervals = circular_intervals(ts, te, length) if target in circular else [(ts, te)]
        if identity >= 0.98:
            query_intervals[query].append((qs, qe))
            target_intervals[target].extend(intervals)
        pair = {
            "query": query, "target": target, "identity": identity,
            "query_coverage": (qe - qs) / qlen,
            "target_coverage": union_length(intervals) / length,
            "reference_span": te - ts, "method": "single_alignment",
        }
        # Doubling a reference must never make a tandem duplication a full match.
        if te - ts > length + math.ceil(length * length_tolerance):
            continue
        def rank(value: dict) -> tuple:
            coverage = min(value["query_coverage"], value["target_coverage"])
            return (coverage >= coverage_min and value["identity"] >= identity_min,
                    coverage, value["identity"])
        if (query, target) not in pairs or rank(pair) > rank(pairs[query, target]):
            pairs[query, target] = pair
    edges: dict[str, list[str]] = defaultdict(list)
    for (query, target), pair in pairs.items():
        if pair["identity"] >= identity_min and min(pair["query_coverage"], pair["target_coverage"]) >= coverage_min:
            edges[query].append(target)
    mapping = bijection(edges, list(candidate))
    one_to_one = len(mapping) == len(truth) == len(candidate)
    truth_bp, candidate_bp = sum(map(len, truth.values())), sum(map(len, candidate.values()))
    extra = candidate_bp - sum(union_length(query_intervals[name]) for name in candidate)
    delta = abs(candidate_bp - truth_bp) / truth_bp
    reverse_mapping = {target: query for query, target in mapping.items()}
    reference_complete = len(mapping) == len(truth)
    strict_accepted = one_to_one and delta <= length_tolerance and extra / truth_bp <= extra_fraction_max
    strict_failures = []
    if not one_to_one:
        strict_failures.append("record_bijection")
    if delta > length_tolerance:
        strict_failures.append("total_length_difference")
    if extra / truth_bp > extra_fraction_max:
        strict_failures.append("extra_sequence_fraction")
    per_chromosome = {}
    for target, reference in truth.items():
        query = reverse_mapping.get(target)
        options = [v for (q, t), v in pairs.items() if t == target]
        best = pairs[query, target] if query is not None else max(
            options, key=lambda x: (min(x["query_coverage"], x["target_coverage"]), x["identity"]), default=None
        )
        per_chromosome[target] = {
            "length": len(reference), "reference_circular": target in circular,
            "matched_query": query, "accepted_pair": query is not None,
            "best_pair": best, "aggregate_coverage": union_length(target_intervals[target]) / len(reference),
        }
    return {
        "accepted": strict_accepted,
        "strict_accepted": strict_accepted,
        "relaxed_accepted": reference_complete,
        "scores": {
            "relaxed": {
                "score": float(reference_complete), "accepted": reference_complete,
                "criteria": "Every supplied reference has a distinct complete passing candidate; extra candidates allowed",
                "failure_reasons": [] if reference_complete else ["reference_record_match"],
            },
            "strict": {
                "score": float(strict_accepted), "accepted": strict_accepted,
                "criteria": "Passing reference/candidate bijection plus total-length and extra-sequence limits",
                "failure_reasons": strict_failures,
            },
            "scale": "0 or 1 acceptance score; not physical-origin or topology confidence",
        },
        "reference_record_match_fraction": len(mapping) / len(truth),
        "unmatched_candidate_records": sorted(set(candidate) - mapping.keys()),
        "thresholds": {
            "identity_min": identity_min, "coverage_min": coverage_min,
            "length_tolerance": length_tolerance, "extra_fraction_max": extra_fraction_max,
        },
        "one_to_one_chromosome_match": one_to_one, "matched_chromosomes": len(mapping),
        "truth_records": len(truth), "candidate_records": len(candidate),
        "truth_bp": truth_bp, "candidate_bp": candidate_bp,
        "relative_length_difference": delta, "extra_sequence_bp": extra,
        "extra_sequence_fraction": extra / truth_bp,
        "per_chromosome": per_chromosome,
        "panel_absent_chromosomes": None if panel_absent is None else {
            "scope": "No conserved organelle-encoded protein panel hit; not absence of all genes",
            "total": len(panel_absent), "matched": sum(name in reverse_mapping for name in panel_absent),
            "identifiers": sorted(panel_absent),
        },
        "physical_topology_certified": False,
        "reference_record_scope": {
            "matched_unit": "supplied reference assembly sequence record",
            "physical_chromosome_count": None,
            "physical_chromosome_count_source": "not established by FASTA records or sequence alignment",
        },
    }


def score(binary: Path, truth_path: Path, candidate_path: Path, circular: set[str], output: Path,
          *, panel_absent: set[str] | None = None) -> dict:
    truth, candidate = read_fasta(truth_path), read_fasta(candidate_path)
    output.mkdir(parents=True, exist_ok=True)
    doubled = output / "alignment-reference.fasta"
    doubled.write_text("".join(f">{name}\n{seq + seq if name in circular else seq}\n" for name, seq in truth.items()))
    if candidate:
        run = subprocess.run(
            [str(binary), "-x", "asm5", "-c", "--eqx", "--secondary=yes", "-N", "50", "-p", "0.1", "-t", "2", str(doubled), str(candidate_path)],
            capture_output=True, text=True, check=True,
        )
        paf = run.stdout
        (output / "alignment.stderr.log").write_text(run.stderr)
    else:
        paf = ""
    (output / "alignment.paf").write_text(paf)
    result = evaluate(truth, candidate, circular, paf, panel_absent=panel_absent)
    result["method"] = "Exact sequence or one minimap2 asm5 alignment spanning both records; circular reference doubled; no physical certification"
    (output / "scored.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--circular-identifiers", type=Path, required=True)
    parser.add_argument("--minimap2", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = score(args.minimap2, args.truth, args.candidate, set(json.loads(args.circular_identifiers.read_text())), args.output)
    print(json.dumps({key: result[key] for key in (
        "accepted", "relaxed_accepted", "strict_accepted", "scores",
        "matched_chromosomes", "relative_length_difference", "extra_sequence_bp",
    )}))
