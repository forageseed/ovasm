"""Collect a batch's results into tables.

  samples.tsv      one row per sample x organelle: run status, time, memory, molecules,
                   circularity, decisiveness and the assessment's headline numbers
  repeats.tsv      one row per repeat call of chromosomes.json (recombining repeats with their
                   pairing counts and minor share)
  genes.tsv        one row per molecule: length, depth in single-copy units, gene call

usage: aggregate.py BATCH_DIR [MANIFEST.tsv]
"""
import csv
import json
import sys
from pathlib import Path


def load(p):
    try:
        return json.loads(Path(p).read_text())
    except (OSError, ValueError):
        return None


def main():
    batch = Path(sys.argv[1])
    meta = {}
    if len(sys.argv) > 2:
        meta = {r['sample']: r for r in csv.DictReader(open(sys.argv[2]), delimiter='\t')}
    samples, repeats, genes = [], [], []
    for d in sorted(p for p in batch.iterdir() if p.is_dir() and not p.name.startswith('.') and not p.name.endswith('.tmp')):
        st = load(d / 'status.json')
        if st is None:
            continue
        m = meta.get(d.name, {})
        for org in ('mitochondrion', 'plastid'):
            s = load(d / 'ovasm' / org / 'summary.json') or {}
            a = load(d / 'ovasm' / org / 'assessment.json') or {}
            ch = s.get('chromosomes') or {}
            cx = a.get('complexity') or {}
            cn = a.get('copy_number') or {}
            jn = a.get('junctions') or {}
            ps = a.get('plastid_structure') or {}
            lengths = s.get('lengths') or []
            samples.append({
                'sample': d.name,
                'species': st.get('species'),
                'source': st.get('source'),
                'clade': m.get('clade', ''),
                'organelle': org,
                'exit_code': st.get('exit_code'),
                'wall_seconds': st.get('wall_seconds'),
                'peak_rss_gb': st.get('peak_rss_gb'),
                'assembled': bool(s),
                'k': s.get('k'),
                'molecules': s.get('molecules'),
                'total_length': sum(x for x in lengths if isinstance(x, int)),
                'longest': max(lengths) if lengths else None,
                'circular': sum(1 for c in (s.get('circular') or []) if c),
                'decisive': s.get('decisive'),
                'chromosomes_min': (ch.get('range') or [None, None])[0],
                'chromosomes_max': (ch.get('range') or [None, None])[-1],
                'single_copy_depth': (a.get('depth') or {}).get('single_copy_depth'),
                'copy_number_agreement': cn.get('agreement'),
                'shallow_bases': (a.get('shallow') or {}).get('bases'),
                'unwalked_bases': (a.get('coverage_of_graph') or {}).get('unwalked_bases'),
                'junctions_used': jn.get('used'),
                'junctions_without_reads': jn.get('without_reads'),
                'junction_min_reads': jn.get('min_reads'),
                'complexity': cx.get('class'),
                'segments': cx.get('segments'),
                'components': cx.get('components'),
                'dead_ends': cx.get('dead_ends'),
                'double_forks': cx.get('double_forks'),
                'recombining_repeats': len(cx.get('recombining_repeats') or []),
                'plastid_structure': ps.get('status'),
                'lsc': ps.get('lsc'),
                'ssc': ps.get('ssc'),
                'ir': ps.get('ir'),
            })
            for c in ch.get('chromosomes') or []:
                for r in c.get('repeats') or []:
                    p = r.get('pairings') or []
                    repeats.append({
                        'sample': d.name, 'organelle': org, 'chromosome': c.get('index'),
                        'segment': r.get('segment'), 'length': r.get('length'), 'call': r.get('call'),
                        'pairing_a': p[0] if len(p) == 2 else None,
                        'pairing_b': p[1] if len(p) == 2 else None,
                        'minor_share': (min(p) / sum(p)) if len(p) == 2 and sum(p) else None,
                    })
            for g in a.get('genes') or []:
                genes.append({
                    'sample': d.name, 'organelle': org, 'molecule': g.get('molecule'),
                    'length': g.get('length'), 'depth_copies': g.get('depth_copies'),
                    'target_genes': len(g.get('target_genes') or []),
                    'other_genes': len(g.get('other_genes') or []), 'call': g.get('call'),
                })
    for name, rows in (('samples.tsv', samples), ('repeats.tsv', repeats), ('genes.tsv', genes)):
        if not rows:
            continue
        with open(batch / name, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
            w.writeheader()
            w.writerows(rows)
        print(f'{name}: {len(rows)} rows', file=sys.stderr)


if __name__ == '__main__':
    main()
