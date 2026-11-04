#!/usr/bin/env python3
"""Formal-source native graph-solver validation with within-library held-out reads.

Synthetic truth is used only for evaluation, never supplied to native assembly/solver.
Real held-out linkage is not independent biological topology truth.
"""
import argparse
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import random
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('path_audit',ROOT/'scripts/audit_path_models.py')
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def read_fasta(p):
    records=[];name=None;seq=[]
    for line in p.read_text().splitlines():
        if line.startswith('>'):
            if name is not None:records.append((name,''.join(seq)))
            name=line[1:];seq=[]
        else:seq.append(line.strip())
    if name is not None:records.append((name,''.join(seq)))
    return records


def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]


def kmers(sequences,k=31):
    return {min(s[i:i+k],rc(s[i:i+k])) for s in sequences for i in range(len(s)-k+1)}


def gfa_data(p):
    nodes={};overlaps={}
    for line in p.read_text().splitlines():
        f=line.split('\t')
        if f[0]=='S':nodes[f[1]]=f[2]
        elif f[0]=='L':
            a,b=f[1]+f[2],f[3]+f[4]
            overlaps[(a,b)]=int(f[5][:-1]);overlaps[(audit.flip(b),audit.flip(a))]=int(f[5][:-1])
    return nodes,overlaps


def check_fasta_spelling(gfa,fasta):
    nodes,edges=gfa_data(gfa)
    def seq(label):return nodes[label[:-1]] if label[-1]=='+' else rc(nodes[label[:-1]])
    checked=0
    for name,dna in read_fasta(fasta):
        tags=dict(t.split('=',1) for t in name.split() if '=' in t)
        path=tags['path'].split(',');expected=seq(path[0])
        for a,b in zip(path,path[1:]):expected+=seq(b)[edges[(a,b)]:]
        if tags.get('circular')=='true':
            overlap=edges[(path[-1],path[0])]
            if overlap:expected=expected[:-overlap]
        assert dna==expected,'FASTA differs from graph walk: '+name
        checked+=1
    return checked


def summarize_fit(report,heldout):
    paths=[]
    for c in report['components']:
        for m in (c.get('best') or {}).get('molecules',[]):paths.append((m['path'].split(),m['circular']))
    evidence=heldout['thread_evidence'];total=evidence['accepted_records'];explained=0
    for w in evidence['paths']:
        if any(audit.contains(p,c,w['path'].split()) for p,c in paths):explained+=w['reads']
    family=[]
    for c in report['components']:
        covers=c.get('witness_configurations') or ([c['best']] if c.get('best') else [])
        for cover in covers:
            family.extend((m['path'].split(),m['circular']) for m in cover['molecules'])
    family_explained=sum(w['reads'] for w in evidence['paths'] if any(audit.contains(p,c,w['path'].split()) for p,c in family))
    envelope=[]
    families=[c['read_path_cover'] for c in report['components'] if c.get('read_path_cover') is not None]
    for c in report['components']:
        covers=c.get('ambiguity_configurations') or c.get('witness_configurations') or ([c['best']] if c.get('best') else [])
        for cover in covers:envelope.extend((m['path'].split(),m['circular']) for m in cover['molecules'])
    envelope_explained=sum(w['reads'] for w in evidence['paths'] if any(audit.contains(p,c,w['path'].split()) for p,c in envelope))
    return dict(eligible_heldout_records=total,explained_heldout_records=explained,
                heldout_linkage_fraction=explained/total if total else None,
                explanatory_family_heldout_records=family_explained,
                explanatory_family_heldout_fraction=family_explained/total if total else None,
                family_size=sum(len(c.get('witness_configurations',[])) for c in report['components']),
                family_minimum_proven=all(c['minimum_family_size_proven'] for c in families) if families else None,
                ambiguity_envelope_heldout_records=envelope_explained,
                ambiguity_envelope_heldout_fraction=envelope_explained/total if total else None,
                ambiguity_envelope_complete=all(c.get('equivalent_families_complete',False) for c in families) if families else None,
                molecules=report['molecules'],decisive=report['decisive'],
                search_truncated=report['search_truncated'],cover_size_limited=report['cover_size_limited'],
                read_model_budget_exhausted=any(c.get('read_model_budget_exhausted',False) for c in report['components']),
                candidates=report['candidates'],partial_paths=len(report['partial_paths']))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary',type=Path,required=True)
    p.add_argument('--source-manifest',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--sample',action='append',default=[],help='name=/absolute/recruited.fastq[.gz]')
    p.add_argument('--synthetic-case',action='append',choices=['none','synthetic-single','synthetic-shared-repeat','synthetic-long-repeat','synthetic-mixed-noisy','synthetic-mixed'],
                   help='Explicit subset; omitted includes all controls, including the noisy upstream-assembly stress test')
    a=p.parse_args();a.binary=a.binary.resolve();a.out=a.out.resolve();a.out.mkdir(parents=True,exist_ok=True)
    os.chdir(ROOT)
    analysis_sources={str(x):sha(x) for x in [Path(__file__).resolve(),ROOT/'scripts/audit_path_models.py']}
    manifest=json.loads(a.source_manifest.read_text())
    for name,h in manifest.items():assert sha(ROOT/name)==h,'source mismatch before validation: '+name
    runs=[]
    def call(name,args,timeout=300):
        t=time.monotonic()
        with (a.out/(name+'.log')).open('w') as f:
            proc=subprocess.run(['/usr/bin/time','-v','timeout',str(timeout),str(a.binary),*map(str,args)],stdout=f,stderr=subprocess.STDOUT)
        runs.append(dict(name=name,args=list(map(str,args)),exit=proc.returncode,wall_seconds=time.monotonic()-t))
        (a.out/'commands.json').write_text(json.dumps(runs,indent=2))
        if proc.returncode:raise RuntimeError(f'{name}: exit {proc.returncode}')
    rng=random.Random(20261104)
    parts={n:''.join(rng.choices('ACGT',k=3000 if n!='R' else 2500)) for n in 'ABCDR'}
    circles=[parts['A']+parts['R']+parts['B'],parts['C']+parts['R']+parts['D']]
    cases=[]
    long_repeat=''.join(rng.choices('ACGT',k=12000))
    recombinant=parts['A']+parts['R']+parts['D']+parts['C']+parts['R']+parts['B']
    synthetic=[('synthetic-single',[circles[0]],False,False),
               ('synthetic-shared-repeat',circles,False,False),
               ('synthetic-long-repeat',[parts['A']+long_repeat+parts['B'],parts['C']+long_repeat+parts['D']],False,True),
               ('synthetic-mixed-noisy',circles+[recombinant],True,True),
               ('synthetic-mixed',circles+[recombinant],False,True)]
    for name,genomes,noisy,ambiguous in synthetic:
        d=a.out/name;d.mkdir(exist_ok=True)
        train=d/'train.fastq';held=d/'heldout.fastq'
        for file,count in [(train,80),(held,40)]:
            with file.open('w') as f:
                for gi,g in enumerate(genomes):
                    for i in range(count//2 if noisy and gi==2 else count):
                        start=rng.randrange(len(g));s=(g+g)[start:start+6500]
                        if rng.randrange(2):s=rc(s)
                        if noisy:
                            s=''.join(rng.choice([b for b in 'ACGT' if b!=ch]) if rng.random()<0.001 else ch for ch in s)
                        f.write(f'@{name}.{file.stem}.{gi}.{i}\n{s}\n+\n'+('I'*len(s))+'\n')
        if a.synthetic_case is None or name in a.synthetic_case:cases.append((name,train,held,genomes,ambiguous))
    for item in a.sample:
        name,path=item.split('=',1)
        if not name or Path(name).name!=name:raise ValueError('sample name must be a single path component')
        source=Path(path).resolve();d=a.out/name;d.mkdir(exist_ok=True)
        train=d/'train.fastq';held=d/'heldout.fastq';counts=[0,0]
        opener=gzip.open if source.suffix=='.gz' else open
        with opener(source,'rt') as src,train.open('w') as tr,held.open('w') as ho:
            while True:
                header=src.readline()
                if not header:break
                seq,plus,qual=src.readline(),src.readline(),src.readline()
                if not header.startswith('@') or not plus.startswith('+') or len(seq.strip())!=len(qual.strip()):raise ValueError('invalid FASTQ record')
                key=header[1:].split()[0].encode();is_held=int.from_bytes(hashlib.sha256(key).digest()[:8],'big')%5==0
                (ho if is_held else tr).write(header+seq+plus+qual);counts[int(is_held)]+=1
        (d/'input.json').write_text(json.dumps(dict(path=str(source),sha256=sha(source),train_records=counts[0],heldout_records=counts[1],scope='within-library ID-grouped holdout; graph assembled from training only'),indent=2))
        cases.append((name,train,held,None,False))
    summary=[]
    for name,train,held,truth,expect_ambiguity in cases:
        d=train.parent
        call(name+'-assemble',['run','hifi','--reads',train,'--read-set','target','--organelle','mitochondrion','--threads','4','--out',d/'native'])
        result=d/'native'
        if not (result/'assembly.gfa').exists():result=result/'mitochondrion'
        g=result/'assembly.gfa';e=result/'evidence.json'
        new=json.loads((result/'linearization.json').read_text());assert new['ranking_model']=='single_record_linkage_v1'
        legacy=json.loads(e.read_text());del legacy['thread_evidence'];(d/'legacy-evidence.json').write_text(json.dumps(legacy))
        call(name+'-legacy',['solve','linearize','--gfa',g,'--evidence',d/'legacy-evidence.json','--out-json',d/'legacy.json','--out-fasta',d/'legacy.fa'])
        call(name+'-heldout',['check','evidence','--gfa',g,'--reads',held,'--out',d/'heldout-evidence.json','--threads','4','--bootstrap','0'])
        ho=json.loads((d/'heldout-evidence.json').read_text());old=json.loads((d/'legacy.json').read_text())
        row=dict(sample=name,scope='synthetic known sequence' if truth else 'within-library held-out linkage, not physical truth',
                 model=summarize_fit(new,ho),legacy=summarize_fit(old,ho),
                 graph_sha256=sha(g),training_sha256=sha(train),heldout_sha256=sha(held),
                 fasta_records_checked=check_fasta_spelling(g,result/'molecules.fasta'))
        if new.get('configurations_fasta'):
            row['configuration_fasta_records_checked']=check_fasta_spelling(g,Path(new['configurations_fasta']))
            assert json.loads((result/'summary.json').read_text())['files']['configurations_fasta']==new['configurations_fasta'],'pipeline summary omitted configuration family'
        p_lines=[line.split('\t') for line in (result/'organelle.gfa').read_text().splitlines() if line.startswith('P\t')]
        assert len({f[1] for f in p_lines})==len(p_lines),'GFA path names collide across configurations'
        for f in p_lines:
            if int(f[1].split('#')[1])>1:assert 'st:Z:unresolved_configuration' in f,'alternative path lost uncertainty provenance'
        row['exported_gfa_paths']=len(p_lines)
        if truth:
            nodes,links=gfa_data(g);k=kmers([x+x[:30] for x in truth]);found=kmers(nodes.values())
            def sequence(label):return nodes[label[:-1]] if label[-1]=='+' else rc(nodes[label[:-1]])
            found |= kmers([sequence(a)[-30:]+sequence(b)[ov:ov+30] for (a,b),ov in links.items()])
            row['truth_31mer_recall']=len(k&found)/len(k)
            assert row['truth_31mer_recall']==1.0,row
            if name in ['synthetic-mixed-noisy','synthetic-mixed']:
                assert row['model']['ambiguity_envelope_heldout_fraction']==1.0,row
            elif not expect_ambiguity:
                assert row['model']['heldout_linkage_fraction'] in [None,1.0],row
            if expect_ambiguity:assert not row['model']['decisive'],row
            if name=='synthetic-long-repeat':
                assert any(c['best']['read_support']['unphased_branch_visits']>0 for c in new['components'] if c.get('best')),row
                assert all(c['read_path_cover']['represented_local_adjacencies']==c['read_path_cover']['supported_local_adjacencies']
                           for c in new['components'] if c.get('read_path_cover')),row
            assert not row['model']['search_truncated'] and not row['model']['read_model_budget_exhausted'],row
        summary.append(row);(a.out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(row),flush=True)
    for name,h in manifest.items():assert sha(ROOT/name)==h,'source changed during validation: '+name
    for name,h in analysis_sources.items():assert sha(Path(name))==h,'analysis source changed during validation: '+name
    (a.out/'provenance.json').write_text(json.dumps(dict(root=str(ROOT),binary=str(a.binary),binary_sha256=sha(a.binary),source_manifest=str(a.source_manifest.resolve()),source_stable=True,synthetic_cases=a.synthetic_case or 'all including noisy assembly stress',analysis_sources=analysis_sources,script=str(Path(__file__).resolve()),script_sha256=sha(__file__)),indent=2))


if __name__=='__main__':main()
