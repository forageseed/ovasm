"""Original-read validation of assembly records (docs/schoenoplectus-read-validation.md).

Input: the records (FASTA) and a minimap2 PAF of the original reads against them, made with
`-c --cs --secondary=no`. Per record: primary-alignment depth, discordant positions (depth >= 5
and >= half of the covering alignments show a mismatch, deletion or an insertion there),
spanning alignments at several flanks (an alignment spans p when it covers [p-f, p+f]), and
reads that continue past a linear end. The same code is applied to every set of records.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np

CS = re.compile(r'(:[0-9]+|\*[a-z][a-z]|[+\-][A-Za-z]+|~[a-z]{2}[0-9]+[a-z]{2})')
END_WINDOW = 2000      # alignments touching the last 2 kb of a linear end
END_CLIP = 1000        # ... that leave at least 1 kb of read beyond it
# An indel of at least this size inside an alignment is a structural difference, not a read
# error (the usual structural-variant boundary; HiFi indel errors are a few bp): the alignment is
# cut there for spanning, and the site is counted separately. Added after the synthetic
# self-test showed a 10-kb misjoin aligned through as one alignment with a 10-kb insertion.
SV_INDEL = 50


def read_fasta(path):
    names, lengths, headers = [], {}, {}
    name = None
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'):
            name = line[1:].split()[0]
            names.append(name)
            lengths[name] = 0
            headers[name] = line[1:]
        elif name is not None:
            lengths[name] += len(line.strip())
    return names, lengths, headers


def intervals(mask):
    """Merged [start, end) runs of True in a boolean array."""
    if not mask.any():
        return []
    d = np.diff(np.concatenate(([0], mask.astype(np.int8), [0])))
    return list(zip(np.flatnonzero(d == 1).tolist(), np.flatnonzero(d == -1).tolist()))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--fasta', type=Path, required=True)
    ap.add_argument('--paf', type=Path, required=True)
    ap.add_argument('--out-json', type=Path, required=True)
    ap.add_argument('--flanks', default='500,1000,5000')
    ap.add_argument('--main-flank', type=int, default=1000)
    ap.add_argument('--min-span', type=int, default=3)
    ap.add_argument('--call-depth', type=int, default=5)
    a = ap.parse_args()
    flanks = [int(x) for x in a.flanks.split(',')]
    names, lengths, headers = read_fasta(a.fasta)
    depth = {n: np.zeros(lengths[n] + 1, np.int64) for n in names}
    disagree = {n: np.zeros(lengths[n] + 1, np.int64) for n in names}
    span = {(n, f): np.zeros(lengths[n] + 1, np.int64) for n in names for f in flanks}
    structural = {n: np.zeros(lengths[n] + 1, np.int64) for n in names}
    ends = {n: {'start': 0, 'end': 0, 'start_continuing': 0, 'end_continuing': 0} for n in names}
    alignments = 0
    reads = set()
    for line in open(a.paf):
        f = line.rstrip('\n').split('\t')
        tags = {t[:2]: t[5:] for t in f[12:]}
        if tags.get('tp', 'P') != 'P':
            continue
        q, ql, qs, qe, strand = f[0], int(f[1]), int(f[2]), int(f[3]), f[4]
        t, tl, ts, te = f[5], int(f[6]), int(f[7]), int(f[8])
        if t not in depth:
            continue
        alignments += 1
        reads.add(q)
        depth[t][ts] += 1
        depth[t][te] -= 1
        # read sequence left beyond the record's start / end (strand-aware)
        before = qs if strand == '+' else ql - qe
        after = ql - qe if strand == '+' else qs
        if ts < END_WINDOW:
            ends[t]['start'] += 1
            ends[t]['start_continuing'] += before >= END_CLIP
        if te > tl - END_WINDOW:
            ends[t]['end'] += 1
            ends[t]['end_continuing'] += after >= END_CLIP
        pos = ts
        dis = disagree[t]
        # One alignment counts at most once per position: a mismatch and an insertion at the same
        # place are one disagreement. (Fixed after the first competitive run: an insertion next to
        # a mismatch landed on the same position in one strand's coordinates and on the next one
        # in the other's, so the same reads were counted twice on one record only.)
        hit = set()
        segments = []          # continuous pieces: no indel of SV_INDEL or more inside
        seg_start = ts
        for op in CS.findall(tags.get('cs', '')):
            c = op[0]
            if c == ':':
                pos += int(op[1:])
            elif c == '*':
                hit.add(pos)
                pos += 1
            elif c == '+':
                hit.add(min(pos, tl - 1))
                if len(op) - 1 >= SV_INDEL:
                    structural[t][min(pos, tl - 1)] += 1
                    segments.append((seg_start, pos))
                    seg_start = pos
            elif c == '-':
                n_del = len(op) - 1
                hit.update(range(pos, pos + n_del))
                if n_del >= SV_INDEL:
                    structural[t][pos] += 1
                    segments.append((seg_start, pos))
                    seg_start = pos + n_del
                pos += n_del
        if pos != te:
            raise SystemExit(f'cs walk ended at {pos}, alignment end {te}, read {q}')
        if hit:
            dis[np.fromiter(hit, np.int64, len(hit))] += 1
        segments.append((seg_start, te))
        for s0, s1 in segments:
            for fl in flanks:
                if s1 - s0 >= 2 * fl:
                    span[(t, fl)][s0 + fl] += 1
                    span[(t, fl)][s1 - fl + 1] -= 1
    records = []
    for n in names:
        L = lengths[n]
        d = np.cumsum(depth[n])[:L]
        dis = disagree[n][:L]
        called = d >= a.call_depth
        discordant = called & (dis * 2 >= np.maximum(d, 1))
        sv = structural[n][:L]
        sv_sites = called & (sv * 2 >= np.maximum(d, 1)) & (sv > 0)
        rec = {
            'record': n, 'length': L, 'header': headers[n],
            'median_depth': float(np.median(d)) if L else 0.0,
            'min_depth': int(d.min()) if L else 0,
            'bp_depth_below_3': int((d < 3).sum()),
            'discordant_positions': int(discordant.sum()),
            'discordant_examples': [
                {'pos': int(p), 'depth': int(d[p]), 'disagree': int(dis[p])}
                for p in np.flatnonzero(discordant)[:20]
            ],
            'structural_discordant_sites': int(sv_sites.sum()),
            'structural_examples': [
                {'pos': int(p), 'depth': int(d[p]), 'large_indels': int(sv[p])}
                for p in np.flatnonzero(sv_sites)[:20]
            ],
            'positions_with_any_large_indel': int((sv > 0).sum()),
            'ends': ends[n],
            'spanning': {},
        }
        for fl in flanks:
            s = np.cumsum(span[(n, fl)])[:L]
            interior = slice(fl, max(fl, L - fl))
            si = s[interior]
            if si.size == 0:
                rec['spanning'][fl] = {'evaluable_bp': 0}
                continue
            zero = intervals(si == 0)
            i = int(si.argmin())
            rec['spanning'][fl] = {
                'evaluable_bp': int(si.size),
                'min': int(si[i]), 'min_at': int(i + fl),
                'median': float(np.median(si)),
                'bp_below_min_span': int((si < a.min_span).sum()),
                'breakpoints': [[int(x + fl), int(y + fl)] for x, y in zero][:50],
                'breakpoint_count': len(zero),
            }
        m = rec['spanning'].get(a.main_flank, {})
        rec['read_validated'] = bool(
            m.get('evaluable_bp', 0) > 0 and m['breakpoint_count'] == 0
            and m['min'] >= a.min_span and rec['discordant_positions'] == 0
            and rec['structural_discordant_sites'] == 0)
        records.append(rec)
    out = {
        'fasta': str(a.fasta), 'paf': str(a.paf), 'alignments_used': alignments,
        'reads_aligned': len(reads), 'flanks': flanks, 'main_flank': a.main_flank,
        'min_span': a.min_span, 'call_depth': a.call_depth,
        'records_validated': sum(r['read_validated'] for r in records),
        'records': records,
    }
    a.out_json.write_text(json.dumps(out, indent=1))
    for r in records:
        m = r['spanning'].get(a.main_flank, {})
        print(f"{r['record'][:28]:28s} len={r['length']:>8d} depth~{r['median_depth']:>6.0f} "
              f"discordant={r['discordant_positions']:>4d} sv={r['structural_discordant_sites']:>3d} "
              f"span{a.main_flank}_min={m.get('min')} "
              f"breaks={m.get('breakpoint_count')} ends(start {r['ends']['start_continuing']}/{r['ends']['start']}, "
              f"end {r['ends']['end_continuing']}/{r['ends']['end']}) validated={r['read_validated']}")


if __name__ == '__main__':
    main()
