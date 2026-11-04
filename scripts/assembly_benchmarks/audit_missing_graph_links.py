#!/usr/bin/env python3
"""Find candidate missing links from split read-to-node PAF alignments.

Candidates are exploratory: repeats, secondary alignments, chimeric reads,
missing references and nuclear insertions can all create apparent support.
Never use this output to certify organelle origin or modify an assembly.
"""
import argparse
from collections import defaultdict
import json
from pathlib import Path

def flip(sign):
    if sign not in ("+", "-"):
        raise ValueError("invalid strand")
    return "-" if sign == "+" else "+"

def canonical(a, ao, b, bo):
    return min((a, ao, b, bo), (b, flip(bo), a, flip(ao)))

def parse_graph(path):
    lengths = {}
    links = set()
    adjacency = defaultdict(list)
    with path.open() as file:
        for line in file:
            f = line.rstrip("\n").split("\t")
            if f[0] == "S":
                if len(f) < 3 or f[2] == "*" or f[1] in lengths:
                    raise ValueError("GFA segments need unique names and sequence")
                lengths[f[1]] = len(f[2])
            elif f[0] == "L":
                if len(f) < 6 or not f[5].endswith("M") or not f[5][:-1].isdigit():
                    raise ValueError("only GFA links with a single M overlap are supported")
                links.add(canonical(f[1], f[2], f[3], f[4]))
                a, b, overlap = (f[1], f[2]), (f[3], f[4]), int(f[5][:-1])
                adjacency[a].append((b, overlap))
                adjacency[(b[0], flip(b[1]))].append(((a[0], flip(a[1])), overlap))
    if not lengths:
        raise ValueError("empty graph")
    if any(a not in lengths or b not in lengths for a, _, b, _ in links):
        raise ValueError("link references missing segment")
    return lengths, links, adjacency

def possible_routes(edge, lengths, adjacency, gap_min, gap_max, tolerance, max_intermediates=4, max_states=5000):
    start, target = (edge[0], edge[1]), (edge[2], edge[3])
    stack = [(start, (start,), 0)]
    routes, explored, depth_limited = [], 0, False
    while stack and explored < max_states:
        node, path, distance = stack.pop()
        explored += 1
        for nxt, overlap in adjacency.get(node, []):
            if nxt in path:
                continue
            cost = distance - overlap + (0 if nxt == target else lengths[nxt[0]])
            next_path = path + (nxt,)
            if nxt == target:
                if len(next_path) > 2 and gap_min - tolerance <= cost <= gap_max + tolerance:
                    routes.append({"path": [n + o for n, o in next_path], "expected_gap": cost})
            elif len(next_path) <= max_intermediates + 1:
                stack.append((nxt, next_path, cost))
            else:
                depth_limited = True
    routes.sort(key=lambda x: (abs(x["expected_gap"] - (gap_min + gap_max) / 2), x["path"]))
    return routes[:5], bool(stack), depth_limited, len(routes)


def audit(gfa, paf, min_align=300, min_identity=0.98, end_slop=50, max_gap=100, max_overlap=50):
    if min_align < 1 or not 0 < min_identity <= 1 or min(end_slop, max_gap, max_overlap) < 0:
        raise ValueError("invalid alignment thresholds")
    lengths, links, adjacency = parse_graph(gfa)
    reads = defaultdict(list)
    query_lengths = {}
    with paf.open() as file:
        for line in file:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                raise ValueError("malformed PAF")
            name, strand, node = f[0], f[4], f[5]
            qlen, qs, qe, tlen, ts, te, matches, block = map(int, (f[1], f[2], f[3], f[6], f[7], f[8], f[9], f[10]))
            if name in query_lengths and query_lengths[name] != qlen:
                raise ValueError("inconsistent lengths for a read name")
            query_lengths[name] = qlen
            if node not in lengths or tlen != lengths[node] or strand not in ("+", "-"):
                raise ValueError("PAF target is inconsistent with graph")
            if not (0 <= qs <= qe <= qlen and 0 <= ts <= te <= tlen and 0 <= matches <= block):
                raise ValueError("invalid PAF coordinates or counts")
            if qe - qs >= min_align and block and matches / block >= min_identity:
                reads[name].append((qs, qe, node, strand, ts, te, tlen))
    support = defaultdict(dict)
    for name, aligns in reads.items():
        aligns.sort()
        for i, a in enumerate(aligns):
            for b in aligns[i + 1:]:
                gap = b[0] - a[1]
                if gap > max_gap:
                    break
                if gap < -max_overlap or a[2] == b[2]:
                    continue
                leaves_end = a[5] >= a[6] - end_slop if a[3] == "+" else a[4] <= end_slop
                enters_end = b[4] <= end_slop if b[3] == "+" else b[5] >= b[6] - end_slop
                if not leaves_end or not enters_end:
                    continue
                edge = canonical(a[2], a[3], b[2], b[3])
                if edge not in links:
                    old = support[edge].get(name)
                    if old is None or abs(gap) < abs(old):
                        support[edge][name] = gap
    candidates = [{"from": edge[0], "from_orientation": edge[1], "to": edge[2], "to_orientation": edge[3],
                   "distinct_reads": len(hits), "read_examples": sorted(hits)[:10],
                   "gap_min": min(hits.values()), "gap_max": max(hits.values()),
                   "reads_near_gap_limit": sum(g >= max_gap - 5 for g in hits.values()),
                   "near_gap_limit": sum(g >= max_gap - 5 for g in hits.values()) * 2 >= len(hits)}
                  for edge, hits in support.items()]
    for c in candidates:
        edge = (c["from"], c["from_orientation"], c["to"], c["to_orientation"])
        routes, truncated, depth_limited, route_count = possible_routes(edge, lengths, adjacency, c["gap_min"], c["gap_max"], max(10, 2 * end_slop))
        c["possible_existing_routes"] = routes
        c["matching_routes_found"] = route_count
        c["route_search_truncated"] = truncated
        c["route_depth_limited"] = depth_limited
    candidates.sort(key=lambda x: (-x["distinct_reads"], x["from"], x["from_orientation"], x["to"], x["to_orientation"]))
    return {"scope": "Candidate graph-discordant links from split alignments; not physical topology or certified organelle origin. Competing loci, repeat ambiguity, mapping errors, structural divergence and chimeras are not ruled out. Possible existing routes are bounded graph-only length checks, not proof that an alignment follows a route.",
            "parameters": {"min_align": min_align, "min_identity": min_identity, "end_slop": end_slop, "max_gap": max_gap, "max_overlap": max_overlap, "route_max_intermediates": 4, "route_max_states": 5000, "route_gap_tolerance": max(10, 2 * end_slop)},
            "reads_with_qualifying_alignments": len(reads), "candidates": candidates}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gfa", type=Path, required=True)
    p.add_argument("--paf", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--min-align", type=int, default=300)
    p.add_argument("--min-identity", type=float, default=.98)
    p.add_argument("--end-slop", type=int, default=50)
    p.add_argument("--max-gap", type=int, default=100)
    p.add_argument("--max-overlap", type=int, default=50)
    args = p.parse_args()
    out = audit(args.gfa, args.paf, args.min_align, args.min_identity, args.end_slop, args.max_gap, args.max_overlap)
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"reads_with_qualifying_alignments": out["reads_with_qualifying_alignments"], "candidates": len(out["candidates"]), "leading_candidates": out["candidates"][:5]}, indent=2))

if __name__ == "__main__":
    main()
