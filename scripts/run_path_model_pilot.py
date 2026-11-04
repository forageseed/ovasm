#!/usr/bin/env python3
"""Server development pilot; runs formal Rust evidence/solver and formal Python audit."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '.agent-work/path-model-pilot-20261104'
OUT.mkdir(exist_ok=True)
os.chdir(ROOT)
spec=importlib.util.spec_from_file_location('audit',ROOT/'scripts/audit_path_models.py')
audit=importlib.util.module_from_spec(spec); spec.loader.exec_module(audit)


def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def source_manifest():
    return {str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/'src').glob('*')) if p.is_file()}


before=source_manifest()
(OUT/'source-before.json').write_text(json.dumps(before,indent=2))
panel=ROOT/'.agent-work/zm4-prior-ranking-20261104'
matched=[]
for manifest in panel.glob('v*-source.sha256'):
    version=manifest.name.split('-')[0]
    if not version[1:].isdigit() or not (panel/(version+'-verified')).is_file(): continue
    old=dict((n.strip(),h) for h,n in (line.split(None,1) for line in manifest.read_text().splitlines()))
    if all((ROOT/n).is_file() and digest(ROOT/n)==h for n,h in old.items()):
        matched.append((int(version[1:]),manifest,panel/(version+'-verified')))
if not matched: raise RuntimeError('no verified formal binary matches all current Rust source hashes')
_,manifest,source_binary=max(matched)
differences=[]
binary=OUT/'ovasm-formal-baseline'
shutil.copy2(source_binary,binary)
expected=digest(source_binary)
assert digest(binary)==expected
(OUT/'provenance.json').write_text(json.dumps(dict(cwd=str(ROOT),binary=str(binary),binary_sha256=expected,
    source_manifest=str(manifest),source_differences=differences,
    source_binary=str(source_binary),
    scope='All recorded Rust source hashes match the formal build. Only evidence and linearize controls are run; full WGS assembly is not validated here.',
    audit_source=str(ROOT/'scripts/audit_path_models.py'),audit_sha256=digest(ROOT/'scripts/audit_path_models.py')),indent=2))


def run(command, name):
    with (OUT/(name+'.log')).open('w') as f:
        subprocess.run(['/usr/bin/time','-v','timeout','90','taskset','-c','40-43','nice','-n','19',*map(str,command)],stdout=f,stderr=subprocess.STDOUT,check=True)


def read_json(p): return json.loads(p.read_text())


def native_case(name, nodes, edges, reads):
    d=OUT/name; d.mkdir(exist_ok=True)
    g=d/'graph.gfa'; f=d/'reads.fa'
    g.write_text(''.join(f'S\t{n}\t{s}\tdp:f:40\n' for n,s in nodes.items())+
                 ''.join(f'L\t{a}\t+\t{b}\t+\t0M\n' for a,b in edges))
    f.write_text(''.join(f'>{name}.{i}\n{s}\n' for i,s in enumerate(reads)))
    e=d/'evidence.json'
    run([binary,'check','evidence','--gfa',g,'--reads',f,'--out',e,'--threads','2','--bootstrap','0'],name+'-evidence')
    return d,g,e


def candidates(g, d, covers):
    c=d/'candidates.json'
    c.write_text(json.dumps(dict(graph=str(g),components=[dict(best=covers[0],alternatives=covers[1:],
                 candidates=len(covers),search_truncated=False,cover_size_limited=False)])))
    return c


def cover(paths):
    return dict(molecules=[dict(path=p,circular=c) for p,c in paths],weak_bridges=0,topological_bridges=0,copy_deviation=0)


rng=random.Random(20261104)
nodes={n:''.join(rng.choices('ACGT',k=1600)) for n in ['A','B','C','D','R','M','S']}
edges=[('A','R'),('C','R'),('R','M'),('M','S'),('S','B'),('S','D'),('B','C'),('D','A')]
real=[('A','R','M','S','B'),('C','R','M','S','D')]
reads=[''.join(nodes[n] for n in p) for p in real for _ in range(10)]
d,g,e=native_case('joint-linkage',nodes,edges,reads)
true=cover([('A+ R+ M+ S+ B+',False),('C+ R+ M+ S+ D+',False)])
false=cover([('A+ R+ M+ S+ D+',False),('C+ R+ M+ S+ B+',False)])
c=candidates(g,d,[false,true]); result=audit.audit(g,e,c)
(d/'audit.json').write_text(json.dumps(result,indent=2))
assert result['components'][0]['top_tied_candidates']==[1],result
scores=result['components'][0]['ordered_candidates']
assert scores[0]['matched_record_walk_occurrences']['local']==scores[1]['matched_record_walk_occurrences']['local']

# Same graph, both pair combinations present equally: never invent a unique winner.
mixed=reads+[''.join(nodes[n] for n in p) for p in [('A','R','M','S','D'),('C','R','M','S','B')] for _ in range(10)]
d,g,e=native_case('coexisting-linkage',nodes,edges,mixed)
c=candidates(g,d,[false,true]); result=audit.audit(g,e,c)
(d/'audit.json').write_text(json.dumps(result,indent=2))
assert result['components'][0]['top_tied_candidates']==[0,1]

# Closure control: native evidence from an open chain, then actual closure-spanning reads.
ring_nodes={n:nodes[n] for n in ['A','B','C']}
ring_edges=[('A','B'),('B','C'),('C','A')]
for close in (False,True):
    name='ring-observed' if close else 'ring-unobserved'
    rs=[nodes['A']+nodes['B']+nodes['C']]*10
    if close: rs += [nodes['B']+nodes['C']+nodes['A']]*10
    d,g,e=native_case(name,ring_nodes,ring_edges,rs)
    c=candidates(g,d,[cover([('A+ B+ C+',False)]),cover([('A+ B+ C+',True)])])
    result=audit.audit(g,e,c); (d/'audit.json').write_text(json.dumps(result,indent=2))
    assert result['components'][0]['top_tied_candidates']==([1] if close else [0,1])

# Let the actual native solver enumerate a repeat-mediated closure with no spanning support.
d,g,e=native_case('native-solver-topology',{n:nodes[n] for n in ['A','B','C','D','R']},
    [('A','B'),('B','R'),('R','A'),('C','R'),('R','D')],
    [nodes['A']+nodes['B']]*20+[nodes['C']+nodes['R']+nodes['D']]*20)
c=d/'linearization.json'
run([binary,'solve','linearize','--gfa',g,'--evidence',e,'--out-fasta',d/'molecules.fa','--out-json',c,
     '--max-molecules','5','--max-alternatives','1000','--max-steps','100000'],'native-solver-topology')
result=audit.audit(g,e,c); (d/'audit.json').write_text(json.dumps(result,indent=2))

# Archived real outputs are development inputs, not fresh product acceptance.
panel=ROOT/'.agent-work/zm4-prior-ranking-20261104'
g=panel/'v16-wgs/mitochondrion/de-novo/assembly.gfa'
c=OUT/'real-native-candidates.json'
run([binary,'solve','linearize','--gfa',g,'--evidence',g.parent/'evidence.json',
     '--out-fasta',OUT/'real-native-molecules.fa','--out-json',c,
     '--max-molecules','12','--max-alternatives','1000','--max-steps','5000000'],'real-native-solver')
native=read_json(c)
assert not native['search_truncated'] and not native['cover_size_limited']
summary=[]
for e in sorted(panel.glob('SRR*.target-graph-evidence.json')):
    try:
        r=audit.audit(g,e,c,20_000_000)
    except ValueError as err:
        summary.append(dict(input=str(e),error=str(err))); continue
    (OUT/(e.stem+'.audit.json')).write_text(json.dumps(r,indent=2))
    summary.append(dict(input=str(e),scope='archived real-data diagnostic; not independent truth',
        observed_walk_types=r['observed_walk_types'],joint_walk_types=r['joint_walk_types'],
        components=[dict(top_ties=len(x['top_tied_candidates']),top_candidates=x['top_tied_candidates'],
                         original_excluded=x['original_excluded_from_top'],search_truncated=x['search_truncated'],
                         cover_size_limited=x['cover_size_limited'],candidate_list_limited=x['candidate_list_limited']) for x in r['components']]))
(OUT/'real-summary.json').write_text(json.dumps(summary,indent=2))
after=source_manifest(); (OUT/'source-after.json').write_text(json.dumps(after,indent=2))
assert before==after,'Concurrent source changes: do not treat this run as a stable formal snapshot'
print(json.dumps(dict(native_controls=5,real_native_candidates=native['candidates'],real=summary),indent=2))
