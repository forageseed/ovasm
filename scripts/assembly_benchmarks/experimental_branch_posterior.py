#!/usr/bin/env python3
"""Experimental conditional Dirichlet model for panel and target traversal counts.

Input is a JSON object with candidates[].pairings[] containing a three-step path and
read_counts keyed by run name. This is not a calibrated physical-topology posterior.
"""
import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path


def flip(s):
    if len(s) < 2 or s[-1] not in "+-":
        raise ValueError("each path step needs an orientation")
    return s[:-1] + ("-" if s[-1] == "+" else "+")


def compute(data, target, strength, cap, minimum_depth):
    if not target or not math.isfinite(strength) or strength < 0 or not math.isfinite(cap) or cap < 0 or minimum_depth < 1:
        raise ValueError("invalid target, strength, cap or minimum panel depth")
    paths = {}
    for row in data["candidates"]:
        for pairing in row["pairings"]:
            path = pairing["path"]
            steps = path.split()
            if len(steps) != 3:
                raise ValueError("each pairing must contain three oriented steps")
            for step in steps:
                flip(step)
            counts = pairing["read_counts"]
            if not isinstance(counts, dict) or any(isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in counts.values()):
                raise ValueError("read counts must be nonnegative integers")
            if path in paths and paths[path] != counts:
                raise ValueError("inconsistent counts for duplicate pairing")
            paths[path] = counts
    for path in paths:
        reverse = " ".join(flip(step) for step in reversed(path.split()))
        if reverse in paths and reverse != path:
            raise ValueError("reverse-complement duplicate pairing may double-count read support")
    samples = sorted({s for counts in paths.values() for s in counts})
    if target not in samples:
        raise ValueError("target sample is absent from input")
    groups = defaultdict(set)
    for path in paths:
        a, mid, b = path.split()
        groups[(a, mid)].add(path)
        groups[(flip(b), flip(mid))].add(path)
    results = []
    for key, choices in sorted(groups.items()):
        choices = sorted(choices)
        if len(choices) < 2:
            continue
        if any(target not in paths[p] for p in choices):
            raise ValueError("target count missing for a branch alternative; missing is not zero")
        target_counts = [paths[p][target] for p in choices]
        donors = []
        for sample in samples:
            if sample == target or any(sample not in paths[p] for p in choices):
                continue
            counts = [paths[p][sample] for p in choices]
            if sum(counts) >= minimum_depth:
                donors.append((sample, counts))
        n = len(choices)
        if len(donors) >= 2:
            mean = [sum((c[i] + 0.5) / (sum(c) + 0.5 * n) for _, c in donors) / len(donors) for i in range(n)]
            applied_strength = strength
            prior_status = "panel_mean_fraction" if strength > 0 else "panel_mean_computed_zero_strength"
        else:
            mean = [1.0 / n] * n
            applied_strength = 0.0
            prior_status = "insufficient_panel_runs_base_prior_only"
        prior = [0.5 + applied_strength * p for p in mean]
        raw_n = sum(target_counts)
        factor = min(1.0, cap / raw_n) if cap > 0 and raw_n else 1.0
        effective = [c * factor for c in target_counts]
        alpha = [a + c for a, c in zip(prior, effective)]
        total = sum(alpha)
        rows = []
        for i, path in enumerate(choices):
            mean_i = alpha[i] / total
            rows.append({"path": path, "panel_mean_fraction": mean[i], "prior_alpha": prior[i], "target_raw_count": target_counts[i], "target_effective_count": effective[i], "posterior_alpha": alpha[i], "conditional_mean": mean_i, "conditional_variance": mean_i * (1.0 - mean_i) / (total + 1.0)})
        results.append({"oriented_entrance_and_repeat": list(key), "prior_status": prior_status, "applied_prior_strength": applied_strength, "panel_runs_eligible": [s for s, _ in donors], "panel_runs_used": [s for s, _ in donors] if len(donors) >= 2 and applied_strength > 0 else [], "target_total_raw": raw_n, "target_total_effective": sum(effective), "alternatives": rows})
    return {"scope": "Experimental conditional model of observed graph-traversal fractions, not a calibrated posterior over physical topologies, chromosome counts, organelle origin, or molecule abundance. Each group conditions on enumerated detected alternatives; unobserved paths and mapping/recruitment biases are omitted. Read assignments, groups and donor runs may be dependent. Panel prior is empirically estimated, not an independent biological truth. With cap > 0 this is a fractional-likelihood generalized posterior, not an ordinary iid-count posterior; without a cap it is a Dirichlet-multinomial conditional posterior only if the stated iid model holds. Marginal variances are model-based and not calibrated confidence intervals.", "target_sample": target, "panel_strength_requested": strength, "target_effective_cap": cap, "minimum_panel_depth": minimum_depth, "jeffreys_base_per_category": 0.5, "groups": results}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--target", required=True)
    p.add_argument("--prior-strength", type=float, required=True, help="sensitivity parameter; not automatically calibrated")
    p.add_argument("--target-effective-cap", type=float, required=True, help="0 uses all target counts; >0 downweights target counts to this total")
    p.add_argument("--minimum-panel-depth", type=int, default=10)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    if a.out.exists():
        raise FileExistsError("refusing to overwrite output")
    raw = a.input.read_bytes()
    result = compute(json.loads(raw), a.target, a.prior_strength, a.target_effective_cap, a.minimum_panel_depth)
    result["input"] = str(a.input)
    result["input_sha256"] = hashlib.sha256(raw).hexdigest()
    with a.out.open("x") as f:
        json.dump(result, f, indent=2)
        f.write("\n")
    print(json.dumps({"groups": len(result["groups"]), "scope": result["scope"]}))


if __name__ == "__main__":
    main()
