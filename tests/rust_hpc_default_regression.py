"""Compare the unmodified default assembler with a published binary on server inputs."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def main():
    ap = argparse.ArgumentParser(__doc__)
    ap.add_argument('--control', type=Path, required=True)
    ap.add_argument('--candidate', type=Path, required=True)
    ap.add_argument('--inputs', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[1]
    assert args.candidate.resolve().is_relative_to(root / 'target')
    args.out.mkdir(parents=True, exist_ok=False)
    receipt = {'formal_source': str(root), 'test_source': str(Path(__file__).resolve()),
               'control': str(args.control.resolve()), 'candidate': str(args.candidate.resolve()),
               'scope': 'default assembler output preservation on previously seen recruited MT/PT inputs',
               'cases': []}
    for name in ['mitochondrion', 'plastid']:
        metadata = json.loads((args.inputs / name / 'assembly.json').read_text())
        paths = []
        for label, binary in [('control', args.control), ('candidate', args.candidate)]:
            out = args.out / name / label
            out.mkdir(parents=True)
            cmd = [str(binary.resolve()), 'assemble', '--reads', str(args.inputs / 'recruit' / f'{name}.fastq'),
                   '--k', str(metadata['k']), '--s', str(metadata['s']), '--min-count', str(metadata['min_count']),
                   '--threads', '8', '--out', str(out / 'graph.gfa'), '--out-json', str(out / 'graph.json')]
            with (out / 'command.log').open('w') as stream:
                p = subprocess.run(cmd, stdout=stream, stderr=subprocess.STDOUT, timeout=600)
            assert p.returncode == 0, (name, label, p.returncode)
            paths.append(out)
        reports = [json.loads((p / 'graph.json').read_text()) for p in paths]
        # Exact output GFA and the existing graph/science fields must agree.
        fields = ['k', 's', 'reads', 'bases', 'sampled_kmers', 'distinct_nodes', 'min_count', 'node_depth',
                  'nodes_kept', 'edges_kept', 'cleaning', 'rescue', 'overlap_joins', 'gaps_filled',
                  'boundary_bp_moved', 'links_blunted', 'links_with_overlap', 'unitigs', 'circular_unitigs',
                  'links', 'total_length', 'n50', 'longest']
        case = {'organelle': name,
                'gfa_byte_identical': (paths[0] / 'graph.gfa').read_bytes() == (paths[1] / 'graph.gfa').read_bytes(),
                'scientific_field_comparison': {k: reports[0][k] == reports[1][k] for k in fields},
                'gfa_sha256': hashlib.sha256((paths[1] / 'graph.gfa').read_bytes()).hexdigest()}
        receipt['cases'].append(case)
    receipt['passed'] = all(c['gfa_byte_identical'] and all(c['scientific_field_comparison'].values()) for c in receipt['cases'])
    (args.out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))
    assert receipt['passed']


if __name__ == '__main__':
    main()
