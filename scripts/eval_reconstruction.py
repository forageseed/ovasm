#!/usr/bin/env python3
"""Score organelle assemblies (a short-read reconstruction, a panel genome used as is, a de novo
short-read assembly ...) against a truth assembly, with two scoring schemes kept side by side.

v1 (as used before 2026-11-14): minimap2 -x asm5 --secondary=no, every alignment counts;
  truth_cov = truth bases under any alignment, query_cov = query bases under any alignment,
  identity = matches / alignment columns over all alignments.
v2:
  * circular truth molecules are doubled, so a query that starts elsewhere on the circle still
    aligns in one piece (positions are taken modulo the molecule); alignments are not joined across
    long gaps (--no-long-join), so a query that joins two truth molecules through a repeat is
    scored piece by piece rather than as one low-identity chain;
  * only primary alignments (tp:A:P) count, and only blocks of >= MIN_BLOCK bp at >= MIN_ID
    identity are "good": recall = truth bases under good blocks, precision = query bases under
    good blocks, identity = matches / columns over primary blocks >= MIN_BLOCK (any identity);
  * copy deficit: truth bases reached only by secondary alignments (a query copy that also
    stands for another truth copy: a repeat the query holds fewer times than the truth);
  * truth completeness: a truth with linear molecules >= MIN_BLOCK is flagged incomplete; with
    the sample's own long reads mapped to the query (--support-depth, `samtools depth -a`
    output), query bases outside good blocks that those reads cover at >= SUPPORT_FRACTION of
    their median depth are "read-supported, absent from truth", and precision_supported counts
    them as correct.
v1 and v2 are written into one row so that a table carries both.

Usage: eval_reconstruction.py --truth T.fa --query Q.fa --label NAME [--support-depth D.tsv]
       [--minimap2 PATH] [--threads N] >> table.tsv   (--header prints the header line only)
"""
from argparse import ArgumentParser
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

MIN_BLOCK = 1000
MIN_ID = 0.99
SUPPORT_FRACTION = 0.2

COLUMNS = [
    "label", "truth_molecules", "truth_bp", "truth_linear_bp", "truth_complete", "query_molecules", "query_bp",
    "query_n50", "v1_truth_cov", "v1_query_cov", "v1_identity", "v2_recall", "v2_precision", "v2_identity",
    "v2_longest_good_block", "v2_copy_deficit_bp", "v2_query_unaligned_bp", "v2_read_supported_absent_bp",
    "v2_precision_supported",
]


def read_fasta(path):
    recs, name, buf, header = [], None, [], ""
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if name is not None:
                recs.append((name, header, "".join(buf)))
            header = line[1:]
            name, buf = header.split()[0], []
        elif line:
            buf.append(line.strip())
    if name is not None:
        recs.append((name, header, "".join(buf)))
    return recs


def circular(header):
    return "circular=true" in header.split()


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def union_len(iv):
    return sum(e - s for s, e in merge(iv))


def fold(s, e, length):
    """Interval [s, e) of a doubled molecule back onto [0, length)."""
    if length <= 0:
        return []
    if e - s >= length:
        return [(0, length)]
    s, e = s % length, s % length + (e - s)
    return [(s, e)] if e <= length else [(s, length), (0, e - length)]


def n50(lengths):
    total, acc = sum(lengths), 0
    for x in sorted(lengths, reverse=True):
        acc += x
        if acc * 2 >= total:
            return x
    return 0


def parse_paf(lines):
    rows = []
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        tags = set(f[12:])
        rows.append({"q": f[0], "qs": int(f[2]), "qe": int(f[3]), "t": f[5], "tl": int(f[6]), "ts": int(f[7]),
                     "te": int(f[8]), "match": int(f[9]), "cols": int(f[10]),
                     "primary": "tp:A:P" in tags or not any(t.startswith("tp:A:") for t in tags)})
    return rows


def score_v1(rows, truth_len, query_len):
    tiv, qiv, m, c = {}, {}, 0, 0
    for r in rows:
        tiv.setdefault(r["t"], []).append((r["ts"], r["te"]))
        qiv.setdefault(r["q"], []).append((r["qs"], r["qe"]))
        m += r["match"]
        c += r["cols"]
    tb = sum(truth_len.values())
    qb = sum(query_len.values())
    return (sum(union_len(v) for v in tiv.values()) / tb if tb else 0.0,
            sum(union_len(v) for v in qiv.values()) / qb if qb else 0.0,
            m / c if c else 0.0)


def score_v2(rows, truth_len, query_len, support=None):
    """rows from a PAF against the doubled truth (with secondaries); truth_len: real lengths."""
    good_t, any_t, good_q, prim_t = {}, {}, {}, {}
    m = c = 0
    longest = 0
    for r in rows:
        ok = r["cols"] >= MIN_BLOCK and r["match"] >= MIN_ID * r["cols"]
        L = truth_len.get(r["t"], 0)
        if ok:
            for iv in fold(r["ts"], r["te"], L):
                any_t.setdefault(r["t"], []).append(iv)
        if not r["primary"]:
            continue
        if r["cols"] >= MIN_BLOCK:
            m += r["match"]
            c += r["cols"]
        if ok:
            longest = max(longest, r["te"] - r["ts"])
            for iv in fold(r["ts"], r["te"], L):
                good_t.setdefault(r["t"], []).append(iv)
            good_q.setdefault(r["q"], []).append((r["qs"], r["qe"]))
    tb = sum(truth_len.values())
    qb = sum(query_len.values())
    recall_bp = sum(union_len(v) for v in good_t.values())
    reach_bp = sum(union_len(v) for v in any_t.values())
    prec_bp = sum(union_len(v) for v in good_q.values())
    supported = 0
    if support:
        for q, L in query_len.items():
            depth = support.get(q)
            if not depth:
                continue
            med = statistics.median(depth) or 0
            if med <= 0:
                continue
            covered = [0] * L
            for s, e in merge(good_q.get(q, [])):
                for i in range(s, min(e, L)):
                    covered[i] = 1
            supported += sum(1 for i in range(min(L, len(depth)))
                             if not covered[i] and depth[i] >= SUPPORT_FRACTION * med)
    return {
        "v2_recall": recall_bp / tb if tb else 0.0,
        "v2_precision": prec_bp / qb if qb else 0.0,
        "v2_identity": m / c if c else 0.0,
        "v2_longest_good_block": longest,
        "v2_copy_deficit_bp": max(0, reach_bp - recall_bp),
        "v2_query_unaligned_bp": qb - prec_bp,
        "v2_read_supported_absent_bp": supported if support else "",
        "v2_precision_supported": (prec_bp + supported) / qb if support and qb else "",
    }


def read_depth(path):
    d = {}
    for line in open(path):
        name, pos, depth = line.split("\t")[:3]
        d.setdefault(name, []).append(int(depth))
    return d


def minimap(mm2, threads, target, query, extra):
    return subprocess.run([mm2, "-c", "-x", "asm5", "-t", str(threads), *extra, target, query],
                          check=True, capture_output=True, text=True).stdout.splitlines()


def evaluate(truth, query, label, mm2="minimap2", threads=4, support=None):
    trecs, qrecs = read_fasta(truth), read_fasta(query)
    truth_len = {n: len(s) for n, _, s in trecs}
    query_len = {n: len(s) for n, _, s in qrecs}
    linear_bp = sum(len(s) for _, h, s in trecs if not circular(h) and len(s) >= MIN_BLOCK)
    row = {"label": label, "truth_molecules": len(trecs), "truth_bp": sum(truth_len.values()),
           "truth_linear_bp": linear_bp, "truth_complete": "yes" if linear_bp == 0 else "no",
           "query_molecules": len(qrecs), "query_bp": sum(query_len.values()), "query_n50": n50(list(query_len.values()))}
    if not qrecs:
        row.update({c: 0 for c in COLUMNS if c not in row})
        return row
    v1 = score_v1(parse_paf(minimap(mm2, threads, truth, query, ["--secondary=no"])), truth_len, query_len)
    row.update(dict(zip(["v1_truth_cov", "v1_query_cov", "v1_identity"], v1)))
    with tempfile.TemporaryDirectory() as d:
        t2 = Path(d) / "truth2.fa"
        t2.write_text("".join(f">{n}\n{s + s if circular(h) else s}\n" for n, h, s in trecs))
        # --no-long-join: asm5 otherwise chains across a gap of tens of kb, so a molecule that joins two
        # subcircles of the truth through a repeat (a conformation the truth graph also lists) aligned as one
        # chain at 0.6-0.8 identity and its bases were not counted (rice 710: recall 0.372 -> 0.946)
        rows = parse_paf(minimap(mm2, threads, str(t2), query, ["-N", "50", "-p", "0.5", "--no-long-join"]))
    row.update(score_v2(rows, truth_len, query_len, read_depth(support) if support else None))
    return row


def fmt(v):
    return f"{v:.5f}" if isinstance(v, float) else str(v)


def main():
    ap = ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--truth")
    ap.add_argument("--query")
    ap.add_argument("--label", default="")
    ap.add_argument("--support-depth")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--header", action="store_true")
    a = ap.parse_args()
    if a.header:
        print("\t".join(COLUMNS))
        return
    if not (a.truth and a.query):
        ap.error("--truth and --query are required")
    row = evaluate(a.truth, a.query, a.label, a.minimap2, a.threads, a.support_depth)
    print("\t".join(fmt(row.get(c, "")) for c in COLUMNS))


if __name__ == "__main__":
    sys.exit(main())
