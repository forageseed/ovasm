"""Finite development comparison on existing recruited inputs; never treats them as truth."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def fasta(path):
    records = []
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            records.append('')
        elif line.strip():
            records[-1] += line.strip().upper()
    return records


def equivalent(a, b):
    if len(a) != len(b):
        return False
    # Sequence representation comparison only; not physical circularity evidence.
    rc = b.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]
    return a in b + b or a in rc + rc


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--binary', type=Path, required=True)
    ap.add_argument('--baseline-directory', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[2]
    binary = args.binary.resolve()
    assert binary.is_relative_to(root / 'target'), 'Use a build from the formal source repository'
    args.out.mkdir(parents=True, exist_ok=False)
    sources = [root / 'Cargo.toml', root / 'Cargo.lock', root / 'build.rs', *sorted((root / 'src').glob('*.rs')), Path(__file__).resolve()]
    hashes = {str(p): sha(p) for p in sources}
    receipt = {'run_id': 'ovasm-independent-rust-hpc-existing-ZM4-v1',
               'formal_source': str(root), 'test_source': str(Path(__file__).resolve()),
               'source_hashes': hashes, 'binary': str(binary), 'binary_sha256': sha(binary),
               'baseline_directory': str(args.baseline_directory.resolve()),
               'scope': 'previously seen recruited HiFi input, fixed existing k/s/cutoff; development comparison, not truth',
               'independent_accuracy_score': None, 'physical_molecule_count': None,
               'physical_topology_identified': False, 'cases': [], 'status': 'running'}

    def save():
        (args.out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    save()

    def step(case, name, command):
        log = args.out / case['organelle'] / f'{name}.log'
        timing = log.with_suffix('.time')
        start = time.monotonic()
        with log.open('w') as stream:
            proc = subprocess.run(['/usr/bin/time', '-v', '-o', str(timing), *map(str, command)],
                                  stdout=stream, stderr=subprocess.STDOUT, timeout=1800)
        match = re.search(r'Maximum resident set size \(kbytes\):\s+(\d+)', timing.read_text())
        result = {'name': name, 'command': list(map(str, command)), 'exit_code': proc.returncode,
                  'wall_seconds': time.monotonic()-start, 'gnu_max_rss_kib': int(match[1]) if match else None,
                  'log': str(log), 'gnu_time': str(timing)}
        case['steps'].append(result)
        save()
        if proc.returncode:
            raise RuntimeError(log.read_text()[-2000:])

    for name in ['mitochondrion', 'plastid']:
        case = {'organelle': name, 'steps': [], 'status': 'running'}
        receipt['cases'].append(case)
        out = args.out / name
        out.mkdir()
        try:
            baseline = args.baseline_directory / name
            metadata = json.loads((baseline / 'assembly.json').read_text())
            input_path = args.baseline_directory / 'recruit' / f'{name}.fastq'
            st = input_path.stat()
            case['input'] = {'path': str(input_path.resolve()), 'size': st.st_size, 'mtime_ns': st.st_mtime_ns,
                             'expected_records': metadata['reads'], 'expected_bases': metadata['bases']}
            case['fixed_existing_parameters'] = {k: metadata[k] for k in ['k', 's', 'min_count']}
            step(case, 'assemble', [binary, 'assemble', '--rust-hpc', '--reads', input_path,
                '--k', metadata['k'], '--s', metadata['s'], '--min-count', metadata['min_count'],
                '--threads', 8, '--out', out / 'graph.gfa', '--out-json', out / 'graph.json'])
            graph = json.loads((out / 'graph.json').read_text())
            assert graph['uses_syncasm'] is False
            assert (graph['reads'], graph['bases']) == (metadata['reads'], metadata['bases'])
            assert input_path.stat().st_mtime_ns == st.st_mtime_ns
            case['graph'] = {k: graph[k] for k in ['reads', 'bases', 'compressed_bases', 'unitigs', 'links', 'total_length', 'consensus']}
            case['segment_length_includes_link_overlaps'] = True
            step(case, 'evidence', [binary, 'evidence', '--gfa', out / 'graph.gfa', '--reads', input_path,
                 '--bootstrap', 0, '--threads', 8, '--out', out / 'evidence.json'])
            ev = json.loads((out / 'evidence.json').read_text())
            assert ev['reads_scanned'] == metadata['reads']
            case['raw_read_evidence'] = {k: v for k, v in ev.items() if not isinstance(v, (list, dict))}
            step(case, 'primary', [binary, 'primary', '--gfa', out / 'graph.gfa',
                 '--out', out / 'primary.fasta', '--out-json', out / 'primary.json'])
            old, new = fasta(baseline / 'molecules.fasta'), fasta(out / 'primary.fasta')
            case['primary_record_lengths'] = list(map(len, new))
            case['baseline_record_lengths'] = list(map(len, old))
            case['baseline_sequence_matches_allowing_orientation_rotation'] = [
                [j for j, b in enumerate(new) if equivalent(a, b)] for a in old]
            case['status'] = 'completed'
        except Exception as exc:
            case['status'] = 'failed'
            case['error'] = str(exc)
        save()
    receipt['source_unchanged'] = all(sha(p) == hashes[str(p)] for p in sources)
    receipt['binary_unchanged'] = sha(binary) == receipt['binary_sha256']
    receipt['status'] = 'completed' if all(c['status'] == 'completed' for c in receipt['cases']) else 'failed'
    save()
    print(json.dumps({k: v for k, v in receipt.items() if k not in ['source_hashes']}, indent=2))


if __name__ == '__main__':
    main()
