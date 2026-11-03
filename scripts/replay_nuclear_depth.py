"""Replay historical/current nuclear peak selection from a saved integer histogram.

Usage: python3 scripts/replay_nuclear_depth.py [fixture.json]
No reads, aligner, reference sequence, or third-party dependencies are needed.
The Rust regression test independently executes the actual current estimator.
"""
import json
import math
import sys
from pathlib import Path


def replay(rows, minimum, local_maxima):
    depth_bin = lambda depth: math.floor(math.log2(depth) / 0.25)
    bins = {}
    for depth, reads, _ in rows:
        if depth >= 3:
            index = depth_bin(depth)
            bins[index] = bins.get(index, 0) + reads
    raw = lambda b: bins.get(b, 0)
    key = lambda b: (raw(b - 1) + raw(b) + raw(b + 1), raw(b))
    candidates = sorted(bins, key=lambda b: (-key(b)[0], -key(b)[1], b))
    counts = {depth: reads for depth, reads, _ in rows}
    for peak in candidates:
        boundary = peak == depth_bin(3) and counts.get(3, 0) > counts.get(4, 0)
        if local_maxima and not (boundary or (key(peak) >= key(peak - 1) and key(peak) > key(peak + 1))):
            continue
        window = sorted(row for row in rows if row[0] >= 3 and abs(depth_bin(row[0]) - peak) <= 2)
        total = sum(row[1] for row in window)
        seen = 0
        for depth, reads, _ in window:
            seen += reads
            if 2 * seen >= total:
                break
        genome_size = sum(row[2] for row in window) // depth
        if genome_size >= minimum:
            return [depth, genome_size]
    return None


def main():
    default = Path(__file__).resolve().parents[1] / "tests/fixtures/lathyrus-nuclear-depth.json"
    fixture = json.loads(Path(sys.argv[1] if len(sys.argv) > 1 else default).read_text())
    result = {"sample": fixture["sample"]}
    for name, local in [("old_unfiltered", False), ("local_maxima", True)]:
        actual = replay(fixture["integer_read_depth_histogram"], fixture["min_nuclear_genome"], local)
        if actual != fixture[name + "_expected"]:
            raise AssertionError(f"{name}: {actual} != {fixture[name + '_expected']}")
        result[name] = {"depth": actual[0], "window_genome_bases": actual[1]}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
