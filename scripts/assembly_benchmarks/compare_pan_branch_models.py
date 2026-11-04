#!/usr/bin/env python3
"""Compare leave-one-run-out predictions of observed branch traversal fractions.

Input is a candidate/pairing JSON with each pairing's path and read_counts by
run. Groups are conditional on candidate paths and the fixed input graph. This
does not estimate a topology posterior or certify biological independence.
"""
import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path

MODELS = ("uniform", "binary_presence", "sample_mean_fraction", "pooled_fraction")

def flip(token):
    if len(token) < 2 or token[-1] not in "+-":
        raise ValueError(f"invalid oriented node: {token!r}")
    return token[:-1] + ("-" if token[-1] == "+" else "+")

def canonical(path):
    steps = path.split()
    if len(steps) != 3:
        raise ValueError(f"expected a three-node traversal: {path!r}")
    reverse = [flip(x) for x in reversed(steps)]
    return tuple(min(steps, reverse))

def evaluate(data, minimum_depth, pseudocount=0.5):
    if not math.isfinite(pseudocount) or pseudocount <= 0:
        raise ValueError("pseudocount must be positive and finite")
    if minimum_depth < 1:
        raise ValueError("minimum depth must be positive")
    paths = {}
    for row in data["candidates"]:
        for p in row["pairings"]:
            key = canonical(p["path"])
            counts = p["read_counts"]
            if any(not isinstance(v, int) or isinstance(v, bool) or v < 0 for v in counts.values()):
                raise ValueError("read counts must be nonnegative integers")
            if key in paths and paths[key] != counts:
                raise ValueError(f"inconsistent counts for {key}")
            paths[key] = counts
    if not paths:
        raise ValueError("no candidate pairings")
    samples = set(next(iter(paths.values())))
    if len(samples) < 3 or any(set(counts) != samples for counts in paths.values()):
        raise ValueError("every pairing must contain counts for the same three or more runs")
    groups = defaultdict(set)
    for a, middle, b in paths:
        groups[(a, middle)].add((a, middle, b))
        groups[(flip(b), flip(middle))].add((a, middle, b))
    events = []
    for held in sorted(samples):
        for group, pathset in sorted(groups.items()):
            ps = sorted(pathset)
            if len(ps) < 2:
                continue
            held_counts = [paths[p][held] for p in ps]
            held_n = sum(held_counts)
            if held_n < minimum_depth:
                continue
            training = []
            for sample in sorted(samples - {held}):
                counts = [paths[p][sample] for p in ps]
                if sum(counts) >= minimum_depth:
                    training.append(counts)
            if len(training) < 2:
                continue
            n = len(ps)
            observed = [x / held_n for x in held_counts]
            present = [sum(counts[i] >= 3 for counts in training) for i in range(n)]
            pooled = [sum(counts[i] for counts in training) for i in range(n)]
            probs = {
                "uniform": [1 / n] * n,
                "binary_presence": [(x + pseudocount) / (sum(present) + n * pseudocount) for x in present],
                "sample_mean_fraction": [sum((c[i] + pseudocount) / (sum(c) + n * pseudocount) for c in training) / len(training) for i in range(n)],
                "pooled_fraction": [(x + pseudocount) / (sum(pooled) + n * pseudocount) for x in pooled],
            }
            metrics = {}
            for model, pred in probs.items():
                metrics[model] = {"cross_entropy": -sum(q * math.log(p) for p, q in zip(pred, observed) if q),
                                  "brier": 1.0 - 2.0 * sum(p * q for p, q in zip(pred, observed)) + sum(p * p for p in pred),
                                  "squared_fraction_error": sum((p - q) ** 2 for p, q in zip(pred, observed)),
                                  "predicted": pred}
            events.append({"heldout": held, "group": group, "paths": [" ".join(p) for p in ps],
                           "heldout_counts": held_counts, "training_runs": len(training), "metrics": metrics})
    if not events:
        raise ValueError("no eligible run/branch groups")
    def summarize(subset):
        return {model: {metric: sum(e["metrics"][model][metric] for e in subset) / len(subset)
                        for metric in ("cross_entropy", "brier", "squared_fraction_error")} for model in MODELS}
    per_sample = {s: summarize([e for e in events if e["heldout"] == s])
                  for s in sorted(samples) if any(e["heldout"] == s for e in events)}
    return {"scope": "Held-out observed traversal fractions on fixed candidate-derived graph groups; not topology truth or a calibrated topology posterior. Runs, groups and reads may be correlated; unobserved and noncandidate paths are omitted.",
            "minimum_depth": minimum_depth, "minimum_training_runs": 2, "pseudocount": pseudocount,
            "event_count": len(events), "events_per_sample": dict(Counter(e["heldout"] for e in events)),
            "event_weighted_summary": summarize(events),
            "sample_weighted_summary": {model: {metric: sum(v[model][metric] for v in per_sample.values()) / len(per_sample)
                                                 for metric in ("cross_entropy", "brier", "squared_fraction_error")} for model in MODELS},
            "per_sample_summary": per_sample, "events": events}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--minimum-depth", type=int, default=10)
    parser.add_argument("--pseudocount", type=float, default=0.5)
    args = parser.parse_args()
    result = evaluate(json.loads(args.input.read_text()), args.minimum_depth, args.pseudocount)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ("events", "per_sample_summary")}, indent=2))

if __name__ == "__main__":
    main()
