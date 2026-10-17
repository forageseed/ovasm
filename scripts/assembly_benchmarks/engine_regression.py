"""Paired control/candidate regression of the assembly engine on existing recruited inputs.

Each case runs the same reads through both binaries: the default product route
(automatic k and min_count, then ``evidence`` and ``linearize``) and the homopolymer-run
encoding (k 1001 runs, the control's automatic min_count, then ``primary``). A published
reference, when one exists, is used only after assembly, to score each output with the
unchanged per-record protocol of ``score_independent_multichrom``. These inputs were seen
during development: the results are development regressions, not blind accuracy.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from score_independent_multichrom import read_fasta, score  # noqa: E402

FIELDS = ['k', 's', 'reads', 'bases', 'min_count', 'min_count_auto', 'node_depth', 'nodes_kept', 'edges_kept',
          'unitigs', 'circular_unitigs', 'links', 'total_length', 'n50', 'longest']


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(binary: Path, argv: list, out: Path, name: str, threads: int) -> dict:
    started = time.monotonic()
    timing = out / f'{name}.time'
    with (out / f'{name}.log').open('w') as log:
        rc = subprocess.call(['/usr/bin/time', '-v', '-o', str(timing), str(binary), *map(str, argv)],
                             stdout=log, stderr=subprocess.STDOUT,
                             env={**os.environ, 'RAYON_NUM_THREADS': str(threads)})
    text = timing.read_text() if timing.exists() else ''
    m = re.search(r'Maximum resident set size \(kbytes\):\s+(\d+)', text)
    step = {'name': name, 'exit_code': rc, 'wall_seconds': time.monotonic() - started,
            'gnu_max_rss_kib': int(m[1]) if m else None}
    if rc:
        raise RuntimeError(f'{name} failed ({rc}): {(out / f"{name}.log").read_text()[-800:]}')
    return step


def fasta_lengths(path: Path) -> list[int]:
    if not path.exists():
        return []
    lengths, n = [], None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            if n is not None:
                lengths.append(n)
            n = 0
        elif n is not None:
            n += len(line.strip())
    if n is not None:
        lengths.append(n)
    return sorted(lengths, reverse=True)


def route(binary: Path, reads: Path, out: Path, encoding: str, min_count: int | None,
          threads: int, auto_min_count: bool = False) -> dict:
    out.mkdir(parents=True)
    gfa, report = out / 'graph.gfa', out / 'graph.json'
    steps = []
    if encoding == 'bases':
        steps.append(run(binary, ['assemble', '--reads', reads, '--threads', threads,
                                  '--out', gfa, '--out-json', report], out, 'assemble', threads))
        steps.append(run(binary, ['evidence', '--gfa', gfa, '--reads', reads, '--out', out / 'evidence.json',
                                  '--k', 21, '--diag-tol', 30, '--min-side-hits', 10, '--min-anchor', 500,
                                  '--bootstrap', 0, '--threads', threads], out, 'evidence', threads))
        steps.append(run(binary, ['linearize', '--gfa', gfa, '--evidence', out / 'evidence.json',
                                  '--out-fasta', out / 'sequences.fasta', '--out-json', out / 'linear.json',
                                  '--min-reads', 2], out, 'linearize', threads))
    else:
        cutoff = [] if auto_min_count else ['--min-count', min_count]
        steps.append(run(binary, ['assemble', '--rust-hpc', '--reads', reads, '--k', 1001, '--s', 31,
                                  *cutoff, '--threads', threads,
                                  '--out', gfa, '--out-json', report], out, 'assemble', threads))
        steps.append(run(binary, ['primary', '--gfa', gfa, '--out', out / 'sequences.fasta',
                                  '--out-json', out / 'linear.json'], out, 'primary', threads))
    r = json.loads(report.read_text())
    return {'steps': steps, 'report': {k: r.get(k) for k in FIELDS},
            'gfa_sha256': sha(gfa), 'sequence_lengths': fasta_lengths(out / 'sequences.fasta')}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--control', type=Path, required=True)
    ap.add_argument('--candidate', type=Path, required=True)
    ap.add_argument('--cases', type=Path, required=True,
                    help='JSON list of {name, organelle, reads, truth|null}')
    ap.add_argument('--minimap2', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--threads', type=int, default=8)
    ap.add_argument('--candidate-auto-min-count', action='store_true',
                    help='candidate homopolymer route omits --min-count (the engine decides); '
                         'the control keeps the base route\'s value')
    a = ap.parse_args()
    root = Path(__file__).resolve().parents[2]
    a.out.mkdir(parents=True, exist_ok=False)
    cases = json.loads(a.cases.read_text())
    receipt = {'formal_source': str(root), 'test_source': str(Path(__file__).resolve()),
               'control': str(a.control.resolve()), 'control_sha256': sha(a.control),
               'candidate': str(a.candidate.resolve()), 'candidate_sha256': sha(a.candidate),
               'candidate_auto_min_count': a.candidate_auto_min_count,
               'scope': __doc__.split('\n')[0], 'reference_use': 'post-assembly scoring only',
               'physical_molecule_count': None, 'physical_topology_identified': False, 'cases': []}

    def save():
        (a.out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')

    for case in cases:
        reads = Path(case['reads'])
        st = reads.stat()
        entry = {**case, 'input_size': st.st_size, 'input_mtime_ns': st.st_mtime_ns, 'routes': {}}
        receipt['cases'].append(entry)
        for encoding in ['bases', 'homopolymer_runs']:
            pair = {}
            for label, binary in [('control', a.control), ('candidate', a.candidate)]:
                out = a.out / case['name'] / encoding / label
                try:
                    cutoff = None
                    if encoding == 'homopolymer_runs':
                        cutoff = entry['routes']['bases']['control']['report']['min_count']
                    result = route(binary, reads, out, encoding, cutoff, a.threads,
                                   auto_min_count=a.candidate_auto_min_count and label == 'candidate')
                    if case.get('truth') and result['sequence_lengths']:
                        # hifiasm-style names: a trailing 'c' marks a circular record, 'l' linear
                        circular = {n for n in read_fasta(Path(case['truth'])) if n.endswith('c')}
                        s = score(a.minimap2, Path(case['truth']), out / 'sequences.fasta', circular,
                                  out / 'score')
                        result['score'] = {k: s[k] for k in (
                            'relaxed_accepted', 'strict_accepted', 'matched_chromosomes', 'truth_records',
                            'candidate_records', 'truth_bp', 'candidate_bp', 'extra_sequence_bp')}
                        result['score']['aggregate_coverage'] = {
                            n: round(v['aggregate_coverage'], 6) for n, v in s['per_chromosome'].items()}
                    result['status'] = 'completed'
                except Exception as error:  # first failures are retained, never retried
                    result = {'status': 'failed', 'error': str(error)}
                pair[label] = result
                entry['routes'][encoding] = pair
                save()
            c, d = pair['control'], pair['candidate']
            pair['gfa_identical'] = c.get('gfa_sha256') == d.get('gfa_sha256') and 'gfa_sha256' in c
            if 'score' in c and 'score' in d:
                pair['no_record_loss'] = d['score']['matched_chromosomes'] >= c['score']['matched_chromosomes']
                pair['no_coverage_loss'] = all(
                    d['score']['aggregate_coverage'][n] >= c['score']['aggregate_coverage'][n] - 1e-4
                    for n in c['score']['aggregate_coverage'])
            save()
        entry['input_unchanged'] = (reads.stat().st_size, reads.stat().st_mtime_ns) == (st.st_size, st.st_mtime_ns)
    save()
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
