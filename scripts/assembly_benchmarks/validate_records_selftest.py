"""Synthetic self-test of validate_records.py: a correct record must pass, a record with a
misjoin (a 10-kb deletion joining 20,000 to 30,000) and one substitution (at 5,000) must fail at
those places, and a record joined to unrelated sequence must show a breakpoint at the join.
Records sit inside a longer molecule (10 kb continues on each side), as real records do.

usage: python validate_records_selftest.py --minimap2 /path/to/minimap2 [--work DIR]
Needs numpy. Writes only under --work (default: a new temporary directory)."""
import argparse
import json
import random
import subprocess
import sys
import tempfile
from pathlib import Path

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--minimap2', required=True)
ap.add_argument('--work', type=Path)
a = ap.parse_args()
work = a.work or Path(tempfile.mkdtemp(prefix='ovasm-validate-selftest-'))
work.mkdir(parents=True, exist_ok=True)
checker = Path(__file__).with_name('validate_records.py')

r = random.Random(7)
rand = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
g = rand(40000)
mol = rand(10000) + g + rand(10000)
bad = list(g[:20000] + g[30000:])
bad[5000] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[bad[5000]]
records = {'good': g, 'bad': ''.join(bad), 'chimera': g[:15000] + rand(15000)}
comp = str.maketrans('ACGT', 'TGCA')
with open(work / 'reads.fq', 'w') as f:
    for i in range(300):
        s = r.randrange(0, len(mol) - 15000)
        seq = mol[s:s + 15000]
        if i % 2:
            seq = seq.translate(comp)[::-1]
        f.write(f'@r{i}\n{seq}\n+\n{"I" * len(seq)}\n')
result = {}
for name, seq in records.items():
    (work / f'{name}.fa').write_text(f'>{name}\n{seq}\n')
    with open(work / f'{name}.paf', 'w') as out:
        subprocess.run([a.minimap2, '-x', 'map-hifi', '-c', '--cs', '--secondary=no', '-t', '4',
                        str(work / f'{name}.fa'), str(work / 'reads.fq')],
                       stdout=out, stderr=subprocess.DEVNULL, check=True)
    subprocess.run([sys.executable, str(checker), '--fasta', str(work / f'{name}.fa'),
                    '--paf', str(work / f'{name}.paf'), '--out-json', str(work / f'{name}.json')],
                   stdout=subprocess.DEVNULL, check=True)
    result[name] = json.loads((work / f'{name}.json').read_text())['records'][0]
good, bad, chim = result['good'], result['bad'], result['chimera']
assert good['read_validated'] and good['discordant_positions'] == 0, good
assert good['spanning']['1000']['breakpoint_count'] == 0
assert good['ends']['start_continuing'] > 0 and good['ends']['end_continuing'] > 0, good['ends']
assert not bad['read_validated']
assert any(abs(x['pos'] - 5000) <= 1 for x in bad['discordant_examples']), bad['discordant_examples']
assert any(x <= 20000 <= y for x, y in bad['spanning']['1000']['breakpoints'])
assert not chim['read_validated']
assert any(x <= 15000 <= y for x, y in chim['spanning']['1000']['breakpoints'])
print(f'validate_records self-test passed ({work})')
