#!/usr/bin/env python3
"""Check an organelle assembly graph's links and repeat traversals against independent long reads.

Works on any GFA (OVASM, Oatk, HiMT, PMAT2 ...). The reads must not have been used to build
the graph (e.g. another sequencing run of the same specimen). Only minimap2 is required.

For every link a->b the reads that align contiguously across the junction (with an anchor on
both sides) support it. For every branching segment X (two or more entries or exits) each
entry x exit traversal a-X-b is tested the same way across the whole of X; a read that spans
several traversals of the same X is assigned to the best one, or is counted ambiguous. Reads
split between two segment ends that the graph does not join are reported as missed links.

Outputs (in --out): links.tsv, traversals.tsv, repeats.tsv, missed_links.tsv, summary.json.
"""

from argparse import ArgumentParser
from collections import defaultdict
import gzip
import json
import math
import os
import subprocess
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def flip(o):
    return "-" if o == "+" else "+"


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


class Graph:
    def __init__(self, path):
        self.seq = {}
        self.out = defaultdict(list)  # (node, orient) -> [((node, orient), overlap)]
        self.links = []  # canonical (a, ao, b, bo, ov, tags)
        seen = set()
        with opener(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] == "S":
                    s = f[2]
                    if s == "*":
                        raise SystemExit(f"{path}: segment {f[1]} has no sequence")
                    self.seq[f[1]] = s.upper()
                elif f[0] == "L":
                    a, ao, b, bo = f[1], f[2], f[3], f[4]
                    cig = f[5] if len(f) > 5 else "*"
                    ov = int(cig[:-1]) if cig.endswith("M") and cig[:-1].isdigit() else 0
                    key = canon_link(a, ao, b, bo)
                    if key in seen:
                        continue
                    seen.add(key)
                    tags = {t.split(":")[0]: t.split(":")[-1] for t in f[6:]}
                    self.links.append((a, ao, b, bo, ov, tags))
                    self.out[(a, ao)].append(((b, bo), ov))
                    self.out[(b, flip(bo))].append(((a, flip(ao)), ov))
        self.inn = defaultdict(list)
        for src, lst in self.out.items():
            for dst, ov in lst:
                self.inn[dst].append((src, ov))
        for k in self.out:
            self.out[k] = sorted(set(self.out[k]))
        for k in self.inn:
            self.inn[k] = sorted(set(self.inn[k]))
        self.names = list(self.seq)  # the GFA's own segment names
        # oriented segment end -> its name in the current graph, for leaving / entering it
        self.map_out = {(s, o): (s, o) for s in self.seq for o in "+-"}
        self.map_in = dict(self.map_out)

    def compact(self):
        """Merge non-branching chains into unitigs so that every tool's graph is compared alike
        (a repeat split into pieces is tested as one segment). Merged names join the parts with
        '|', e.g. 'u8.1+|u8.2+'."""
        def unique_next(e, chain):
            outs = self.out.get(e, [])
            if len(outs) != 1:
                return None
            nxt, ov = outs[0]
            if len(self.inn.get(nxt, [])) != 1 or nxt[0] in chain:
                return None
            return nxt, ov
        done, seq, map_out, map_in = set(), {}, {}, {}
        for s in list(self.seq):
            if s in done:
                continue
            chain = {s}
            fwd, e = [], (s, "+")
            while (n := unique_next(e, chain)) is not None:
                fwd.append(n)
                chain.add(n[0][0])
                e = n[0]
            back, e = [], (s, "-")
            while (n := unique_next(e, chain)) is not None:
                back.append(n)
                chain.add(n[0][0])
                e = n[0]
            # elements e0..ek in forward orientation, each with its overlap to the previous one;
            # back holds (s-)->(p1), (p1)->(p2) ... so p_i's overlap with s's side is back[i-1]
            rb = list(reversed(back))
            elems = [((n, flip(o)), rb[t - 1][1] if t else 0) for t, ((n, o), _) in enumerate(rb)]
            elems.append(((s, "+"), rb[-1][1] if rb else 0))
            elems += [(nxt, ov) for nxt, ov in fwd]
            name = s if len(elems) == 1 else "|".join(f"{n}{o}" for (n, o), _ in elems)
            parts = [self.oseq(n, o)[ov:] for (n, o), ov in elems]
            seq[name] = "".join(parts)
            done |= {n for (n, _), _ in elems}
            e0, ek = elems[0][0], elems[-1][0]
            map_out[ek], map_out[(e0[0], flip(e0[1]))] = (name, "+"), (name, "-")
            map_in[e0], map_in[(ek[0], flip(ek[1]))] = (name, "+"), (name, "-")
        links, seen = [], set()
        for a, ao, b, bo, ov, tags in self.links:
            src, dst = map_out.get((a, ao)), map_in.get((b, bo))
            if src is None or dst is None:
                continue  # internal to a chain
            key = canon_link(*src, *dst)
            if key not in seen:
                seen.add(key)
                links.append((*src, *dst, ov, tags))
        self.seq, self.links = seq, links
        self.map_out, self.map_in = map_out, map_in
        self.out, self.inn = defaultdict(list), defaultdict(list)
        for a, ao, b, bo, ov, _ in links:
            self.out[(a, ao)].append(((b, bo), ov))
            self.out[(b, flip(bo))].append(((a, flip(ao)), ov))
        for src, lst in list(self.out.items()):
            self.out[src] = sorted(set(lst))
            for dst, ov in self.out[src]:
                self.inn[dst].append((src, ov))

    def oseq(self, node, o):
        s = self.seq[node]
        return s if o == "+" else rc(s)


def canon_link(a, ao, b, bo):
    x = f"{a}{ao}>{b}{bo}"
    y = f"{b}{flip(bo)}>{a}{flip(ao)}"
    return min(x, y)


def canon_trav(a, ao, x, xo, b, bo):
    if xo == "+":
        return f"{a}{ao} {x}+ {b}{bo}"
    return f"{b}{flip(bo)} {x}+ {a}{flip(ao)}"


def write_fasta(path, recs):
    with open(path, "w") as fh:
        for name, s in recs:
            fh.write(f">{name}\n{s}\n")


def minimap2(args, ref, reads, threads, extra):
    cmd = [args.minimap2, "-x", args.preset, "-c", "-t", str(threads)] + extra + [ref] + reads
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    for line in p.stdout:
        f = line.rstrip("\n").split("\t")
        tags = {}
        for t in f[12:]:
            k, typ, v = t.split(":", 2)
            tags[k] = v
        yield {"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4],
               "t": f[5], "tlen": int(f[6]), "ts": int(f[7]), "te": int(f[8]),
               "match": int(f[9]), "alen": int(f[10]), "tp": tags.get("tp", "P"),
               "AS": int(tags.get("AS", 0))}
    if p.wait() != 0:
        raise SystemExit(f"minimap2 failed: {' '.join(cmd)}")


def decoy_reads(args, reads):
    """Reads that lie (almost) entirely on the decoy sequence, e.g. the plastid genome."""
    if not args.decoy:
        return set()
    cover = defaultdict(int)
    qlen = {}
    for h in minimap2(args, args.decoy, reads, args.threads, ["--secondary=no"]):
        if h["match"] / max(h["alen"], 1) >= args.min_identity:
            cover[h["q"]] += h["qe"] - h["qs"]
            qlen[h["q"]] = h["qlen"]
    return {q for q, c in cover.items() if c >= args.decoy_fraction * qlen[q]}


def anchor_len(args, flank):
    """The anchor a read must cover on a flank; a flank segment shorter than the anchor needs
    all but its last end_slack bases (a quarter for very short ones)."""
    return max(1, min(args.anchor, flank - min(args.end_slack, flank // 4)))


def build_probes(g, args):
    """Link probes, traversal probes and their spanning windows."""
    F = args.flank
    probes = {}  # name -> dict(seq, lo, hi, kind, key, group)
    for a, ao, b, bo, ov, tags in g.links:
        left = g.oseq(a, ao)[-F:]
        right = g.oseq(b, bo)[ov:][:F]
        if not right:
            continue
        lf = len(left)
        anc_l, anc_r = anchor_len(args, lf), anchor_len(args, len(right))
        key = canon_link(a, ao, b, bo)
        probes[f"L{len(probes)}"] = dict(seq=left + right, lo=lf - anc_l, hi=lf + anc_r, kind="link",
                                         key=key, group=key, junction=lf, anchors=(anc_l, anc_r))
    trav = {}
    for x in g.seq:
        ins, outs = g.inn[(x, "+")], g.out[(x, "+")]
        if not ins or not outs or (len(ins) < 2 and len(outs) < 2):
            continue
        xs = g.seq[x]
        for (a, ao), ov1 in ins:
            for (b, bo), ov2 in outs:
                mid = xs[ov1:]
                right = g.oseq(b, bo)[ov2:][:F]
                left = g.oseq(a, ao)[-F:]
                if not mid or not right:
                    continue
                key = canon_trav(a, ao, x, "+", b, bo)
                if key in trav:
                    continue
                lf = len(left)
                anc_l, anc_r = anchor_len(args, lf), anchor_len(args, len(right))
                name = f"T{len(trav)}"
                trav[key] = name
                probes[name] = dict(seq=left + mid + right, lo=lf - anc_l, hi=lf + len(mid) + anc_r,
                                    kind="traversal", key=key, group=x, mid=len(mid),
                                    entry=f"{a}{ao}", exit=f"{b}{bo}", anchors=(anc_l, anc_r),
                                    junction=lf)
    return probes


def assign(args, probes, reads, skip, workdir):
    ref = os.path.join(workdir, "probes.fa")
    write_fasta(ref, [(n, p["seq"]) for n, p in probes.items()])
    # every chain is kept (-P) so that a read is compared against all traversals of a segment
    best = defaultdict(list)  # (read, kind, group) -> [(identity, AS, probe)]
    near = defaultdict(set)   # probe -> reads with an anchor reaching the junction
    for h in minimap2(args, ref, reads, args.threads, ["-P"]):
        if h["q"] in skip:
            continue
        ident = h["match"] / max(h["alen"], 1)
        if ident < args.min_identity:
            continue
        p = probes[h["t"]]
        j = p["junction"]
        # a link's opportunity: reads with a full anchor that reach its junction from either side
        if p["kind"] == "link" and ((h["ts"] <= p["lo"] and h["te"] >= j - 20)
                                    or (h["ts"] <= j + 20 and h["te"] >= p["hi"])):
            near[h["t"]].add(h["q"])
        # A read on another path can share the first part of a flank and still pass the anchor;
        # its alignment then stops while both the read and the probe go on. Supporting alignments
        # must reach the probe's end or the read's end on both sides.
        if h["strand"] == "+":
            read_left, read_right = h["qs"], h["qlen"] - h["qe"]
        else:
            read_left, read_right = h["qlen"] - h["qe"], h["qs"]
        stops_early = (min(h["ts"], read_left) > args.end_slack
                       or min(len(p["seq"]) - h["te"], read_right) > args.end_slack)
        if h["ts"] <= p["lo"] and h["te"] >= p["hi"] and not stops_early:
            best[(h["q"], p["kind"], p["group"])].append((ident, h["AS"], h["t"]))
    support = defaultdict(int)
    ambiguous = defaultdict(int)
    for (q, kind, group), cands in best.items():
        per = {}
        for ident, AS, t in cands:
            if t not in per or ident > per[t][0]:
                per[t] = (ident, AS)
        ranked = sorted(per.items(), key=lambda kv: (-kv[1][0], -kv[1][1]))
        if len(ranked) > 1 and ranked[0][1][0] - ranked[1][1][0] < args.tie_identity:
            for t, _ in ranked:
                if ranked[0][1][0] - per[t][0] < args.tie_identity:
                    ambiguous[t] += 1
            continue
        support[ranked[0][0]] += 1
    return support, ambiguous, near


def segment_splits(args, g, reads, skip, workdir):
    """Adjacencies seen in reads split across segment ends that the graph does not connect."""
    ref = os.path.join(workdir, "segments.fa")
    write_fasta(ref, sorted(g.seq.items()))
    hits = defaultdict(list)
    for h in minimap2(args, ref, reads, args.threads, ["--secondary=no"]):
        if h["q"] in skip or h["match"] / max(h["alen"], 1) < args.min_identity:
            continue
        if h["qe"] - h["qs"] < min(args.split_min_aln, 0.8 * h["tlen"]):
            continue
        hits[h["q"]].append(h)
    max_ov = max([l[4] for l in g.links] + [0])
    found = defaultdict(set)
    for q, hs in hits.items():
        hs.sort(key=lambda h: h["qs"])
        for h1, h2 in zip(hs, hs[1:]):
            gap = h2["qs"] - h1["qe"]
            if gap > args.split_slack or gap < -(max_ov + args.split_slack):
                continue
            o1 = "+" if h1["strand"] == "+" else "-"
            o2 = "+" if h2["strand"] == "+" else "-"
            end1 = h1["tlen"] - h1["te"] if o1 == "+" else h1["ts"]
            start2 = h2["ts"] if o2 == "+" else h2["tlen"] - h2["te"]
            if end1 > args.split_slack or start2 > args.split_slack:
                continue
            a, b = (h1["t"], o1), (h2["t"], o2)
            if connected(g, a, b, gap + end1 + start2 + args.split_slack):
                continue
            found[canon_link(a[0], a[1], b[0], b[1])].add(q)
    return found


def connected(g, a, b, budget):
    """Is there a walk a -> ... -> b whose interior spells at most `budget` bases?"""
    stack = [(a, 0, 0)]
    while stack:
        node, used, depth = stack.pop()
        for nxt, ov in g.out[node]:
            if nxt == b:
                return True
            cost = used + len(g.seq[nxt[0]]) - ov
            if depth < 6 and cost <= budget:
                stack.append((nxt, cost, depth + 1))
    return False


def restrict(g, args):
    """Keep only segments mostly covered by --restrict-fasta (e.g. a tool's own final mitogenome)."""
    q = os.path.join(args.out, "all_segments.fa")
    write_fasta(q, sorted(g.seq.items()))
    cover = defaultdict(list)
    cmd = [args.minimap2, "-x", "asm20", "-t", str(args.threads), "--secondary=no", args.restrict_fasta, q]
    for line in subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.splitlines():
        f = line.split("\t")
        cover[f[0]].append((int(f[2]), int(f[3])))
    keep = set()
    for s, iv in cover.items():
        iv.sort()
        tot, end = 0, -1
        for a, b in iv:
            if b > end:
                tot += b - max(a, end)
                end = b
        if tot >= args.restrict_cover * len(g.seq[s]):
            keep.add(s)
    g.seq = {s: v for s, v in g.seq.items() if s in keep}
    g.links = [l for l in g.links if l[0] in keep and l[2] in keep]
    for d in (g.out, g.inn):
        for k in list(d):
            if k[0] not in keep:
                del d[k]
            else:
                d[k] = [x for x in d[k] if x[0][0] in keep]


def length_bin(n):
    for hi, name in ((1000, "<1kb"), (5000, "1-5kb"), (10000, "5-10kb"), (20000, "10-20kb")):
        if n < hi:
            return name
    return ">=20kb"


def missing_share_upper(n):
    return 1.0 if n == 0 else 1.0 - 0.05 ** (1.0 / n)


def call_status(reads, opp, args):
    """opp: independent reads that could have shown it (the junction's or segment's spanners).
    A graph may keep a minor form (a few percent of molecules), so no read among a few dozen
    is no contradiction: that needs the 95% upper bound of its share below --max-absent-share."""
    if reads >= args.min_support:
        return "supported"
    if reads > 0:
        return "weak"
    if opp == 0:
        return "untestable"
    return "contradicted" if missing_share_upper(opp) < args.max_absent_share else "not_seen"


def load_evidence(path, strip, g):
    if not path:
        return None
    ev = json.load(open(path))

    def name(x):
        n, o = x[:-1], x[-1]
        return (strip(n), o)

    # names after compaction: leaving a segment end / entering one (None inside a merged chain)
    left = lambda e: g.map_out.get(e)
    enter = lambda e: g.map_in.get(e)
    links, travs, reps = {}, {}, []
    for l in ev.get("links", []):
        a, b = left(name(l["from"])), enter(name(l["to"]))
        if a and b:
            links[canon_link(*a, *b)] = l.get("reads", 0)
    for t in ev.get("traversals", []):
        p = [name(x) for x in t["path"].split()]
        if len(p) == 3 and left(p[0]) and enter(p[1]) and enter(p[2]):
            travs[canon_trav(*left(p[0]), *enter(p[1]), *enter(p[2]))] = t.get("reads", 0)
    for r in ev.get("repeats", []):
        x = enter(name(r["repeat"]))
        e = [left(name(v)) for v in r["entries"]]
        o = [enter(name(v)) for v in r["exits"]]
        if x and all(e) and all(o):
            reps.append(dict(repeat=x, entries=e, exits=o, pairing_reads=r["pairing_reads"]))
    return dict(links=links, traversals=travs, repeats=reps)


def main():
    ap = ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--reads", required=True, nargs="+", help="independent reads (FASTQ/FASTA, .gz ok)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--evidence", help="OVASM evidence.json for the same graph (optional)")
    ap.add_argument("--decoy", help="FASTA; reads lying on it entirely are ignored (e.g. plastid)")
    ap.add_argument("--decoy-fraction", type=float, default=0.9)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--preset", default="map-hifi")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--flank", type=int, default=2000)
    ap.add_argument("--anchor", type=int, default=500)
    ap.add_argument("--min-identity", type=float, default=0.97)
    ap.add_argument("--tie-identity", type=float, default=0.001)
    ap.add_argument("--end-slack", type=int, default=150,
                    help="an alignment may end this far short of both the probe end and the read end")
    ap.add_argument("--min-support", type=int, default=2)
    ap.add_argument("--max-absent-share", type=float, default=0.05,
                    help="a link or traversal without reads is contradicted when its share's 95%% upper "
                         "bound is below this (>=59 spanning reads at 0.05)")
    ap.add_argument("--max-missing-share", type=float, default=0.05,
                    help="a 2x2 repeat is 'apart' when the absent pairing's 95%% upper share is below this")
    ap.add_argument("--split-min-aln", type=int, default=1000)
    ap.add_argument("--split-slack", type=int, default=300)
    ap.add_argument("--restrict-fasta", help="evaluate only segments covered by this FASTA (asm20)")
    ap.add_argument("--restrict-cover", type=float, default=0.5)
    ap.add_argument("--no-compact", action="store_true", help="do not merge non-branching chains")
    ap.add_argument("--min-segment", type=int, default=0,
                    help="ignore links and traversals touching segments shorter than this")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    g = Graph(args.gfa)

    names = set(g.names)

    def strip(n):  # evidence.json may drop the GFA's '<sample>.' prefix
        if n in names:
            return n
        for s in g.names:
            if s.endswith("." + n):
                return s
        return n

    if args.restrict_fasta:
        restrict(g, args)
    if not args.no_compact:
        g.compact()
    ev = load_evidence(args.evidence, strip, g)
    skip = decoy_reads(args, args.reads)
    probes = build_probes(g, args)
    support, ambiguous, near = assign(args, probes, args.reads, skip, args.out)
    splits = segment_splits(args, g, args.reads, skip, args.out)

    seglen = {s: len(v) for s, v in g.seq.items()}
    small = lambda *ns: any(seglen[n] < args.min_segment for n in ns)

    rows_l = []
    for n, p in probes.items():
        if p["kind"] != "link":
            continue
        a, b = p["key"].split(">")
        if small(a[:-1], b[:-1]):
            continue
        r = support.get(n, 0)
        opp = len(near.get(n, ()))
        rows_l.append(dict(link=p["key"], heldout_reads=r, opportunity=opp, ambiguous=ambiguous.get(n, 0),
                           status=call_status(r, opp, args),
                           ovasm_reads=(ev["links"].get(p["key"]) if ev else None)))
    rows_t = []
    groups = defaultdict(list)
    for n, p in probes.items():
        if p["kind"] != "traversal":
            continue
        groups[p["group"]].append(n)
    for x, names in groups.items():
        total = sum(support.get(n, 0) for n in names)
        amb = sum(ambiguous.get(n, 0) for n in names)
        for n in names:
            p = probes[n]
            r = support.get(n, 0)
            rows_t.append(dict(traversal=p["key"], segment=x, length=p["mid"], bin=length_bin(p["mid"]),
                               heldout_reads=r, group_spanning=total, group_ambiguous=amb,
                               status=call_status(r, total, args),
                               ovasm_reads=(ev["traversals"].get(p["key"], 0) if ev else None)))
    tindex = {r["traversal"]: r for r in rows_t}
    rows_r = []
    if ev:
        for rp in ev["repeats"]:
            (x, xo), e, o = rp["repeat"], rp["entries"], rp["exits"]
            if len(e) != 2 or len(o) != 2:
                continue

            def tr(a, b):
                key = canon_trav(*a, x, xo, *b)
                return tindex.get(key, {}).get("heldout_reads", 0)
            p1 = tr(e[0], o[0]) + tr(e[1], o[1])
            p2 = tr(e[0], o[1]) + tr(e[1], o[0])
            rows_r.append(dict(repeat=f"{x}{xo}", length=seglen[x], bin=length_bin(seglen[x]),
                               ovasm_pairings=rp["pairing_reads"], ovasm_call=pair_call(*rp["pairing_reads"], args),
                               heldout_pairings=[p1, p2], heldout_call=pair_call(p1, p2, args),
                               ovasm_minor=minor(*rp["pairing_reads"]), heldout_minor=minor(p1, p2)))
    rows_m = []
    for key, qs in sorted(splits.items(), key=lambda kv: -len(kv[1])):
        a, b = key.split(">")
        if small(a[:-1], b[:-1]):
            continue
        rows_m.append(dict(link=key, reads=len(qs), status="missed" if len(qs) >= args.min_support else "single"))

    def tsv(path, rows):
        with open(path, "w") as fh:
            if not rows:
                return
            cols = list(rows[0])
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(json.dumps(r[c]) if isinstance(r[c], (list, dict)) else str(r[c]) for c in cols) + "\n")
    tsv(os.path.join(args.out, "links.tsv"), rows_l)
    tsv(os.path.join(args.out, "traversals.tsv"), rows_t)
    tsv(os.path.join(args.out, "repeats.tsv"), rows_r)
    tsv(os.path.join(args.out, "missed_links.tsv"), rows_m)

    def count(rows, key="status"):
        c = defaultdict(int)
        for r in rows:
            c[r[key]] += 1
        return dict(c)

    sl, st = count(rows_l), count(rows_t)
    missed = sum(1 for r in rows_m if r["status"] == "missed")
    sup_l, con_l = sl.get("supported", 0), sl.get("contradicted", 0)
    summary = dict(
        label=args.label, gfa=os.path.abspath(args.gfa), reads=[os.path.abspath(r) for r in args.reads],
        segments=len(g.seq), graph_bases=sum(seglen.values()), links=len(rows_l), link_status=sl,
        link_precision=ratio(sup_l, sup_l + con_l), missed_links=missed,
        link_recall=ratio(sup_l, sup_l + missed),
        branching_segments=len(groups), traversals=len(rows_t), traversal_status=st,
        traversal_precision=ratio(st.get("supported", 0), st.get("supported", 0) + st.get("contradicted", 0)),
        traversals_by_bin=bin_table(rows_t), decoy_reads_ignored=len(skip),
        params={k: v for k, v in vars(args).items() if k not in ("reads", "gfa", "out")})
    if ev:
        summary["ovasm_vs_heldout"] = compare_ovasm(rows_l, rows_t, rows_r)
    json.dump(summary, open(os.path.join(args.out, "summary.json"), "w"), indent=1)
    print(json.dumps({k: summary[k] for k in ("label", "segments", "links", "link_status", "link_precision",
                                               "missed_links", "link_recall", "traversals", "traversal_status")}))


def ratio(a, b):
    return round(a / b, 4) if b else None


def minor(p1, p2):
    return round(min(p1, p2) / (p1 + p2), 4) if p1 + p2 else None


def pair_call(p1, p2, args):
    n = p1 + p2
    if min(p1, p2) > 0:
        return "interconvertible"
    if n == 0:
        return "not_crossed"
    return "apart" if missing_share_upper(n) < args.max_missing_share else "too_few_reads"


def bin_table(rows):
    t = defaultdict(lambda: defaultdict(int))
    for r in rows:
        t[r["bin"]][r["status"]] += 1
    return {k: dict(v) for k, v in t.items()}


def compare_ovasm(rows_l, rows_t, rows_r):
    """OVASM's own read calls (recruited reads) against the independent reads."""
    out = {}
    for name, rows in (("links", rows_l), ("traversals", rows_t)):
        c = defaultdict(int)
        for r in rows:
            o = r["ovasm_reads"]
            if o is None:
                continue
            c[f"ovasm_{'yes' if o > 0 else 'no'}__heldout_{r['status']}"] += 1
        out[name] = dict(c)
    c = defaultdict(int)
    for r in rows_r:
        c[f"{r['ovasm_call']}->{r['heldout_call']}"] += 1
    out["repeat_calls"] = dict(c)
    return out


if __name__ == "__main__":
    main()
