#!/usr/bin/env python3
"""Read-only sequence concordance of existing HiFi and WGS+pan assemblies.

Exact canonical k-mer agreement measures sequence content, not base QV or complete
physical topology. Circular records include words crossing the chosen origin.
"""
import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import subprocess


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def fasta(p):
    records=[];header=None;seq=[]
    with p.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if header is not None:records.append(record(header,''.join(seq)))
                header=line[1:].strip();seq=[]
            else:seq.append(line.strip().upper())
    if header is not None:records.append(record(header,''.join(seq)))
    return records


def record(header,seq):
    tags=dict(t.split('=',1) for t in header.split() if '=' in t)
    assert set(seq)<=set('ACGT'),'non-ACGT assembly needs explicit ambiguous-word treatment'
    if 'length' in tags:assert int(tags['length'])==len(seq)
    return dict(header=header,seq=seq,circular=tags.get('circular')=='true',configuration=tags.get('configuration'))


def rc(s):
    return s.translate(str.maketrans('ACGT','TGCA'))[::-1]


def words(records,k):
    result=set()
    for rec in records:
        s=rec['seq'];s+=s[:k-1] if rec['circular'] else '';reverse=rc(s);n=len(s)
        result.update(min(s[i:i+k],reverse[n-i-k:n-i]) for i in range(n-k+1))
    return result


def agreement(reference,query):
    shared=len(reference&query)
    return dict(hifi_distinct_words=len(reference),wgs_distinct_words=len(query),shared_words=shared,
                hifi_words_recovered=shared/len(reference) if reference else None,
                wgs_words_found_in_hifi=shared/len(query) if query else None,
                hifi_only_words=len(reference-query),wgs_only_words=len(query-reference))


def exact_molecule_match(records,hifi):
    if len(records)!=len(hifi):return False
    used=set()
    for rec in records:
        s=rec['seq'];match=None
        for i,h in enumerate(hifi):
            if i in used or len(h['seq'])!=len(s):continue
            if h['circular']:
                doubled=h['seq']*2
                same=s in doubled or rc(s) in doubled
            else:same=s==h['seq'] or rc(s)==h['seq']
            if same:match=i;break
        if match is None:return False
        used.add(match)
    return True


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--hifi-root',type=Path,required=True)
    p.add_argument('--wgs-run',action='append',required=True,help='label=/existing/pan/output')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--k',type=int,default=31)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    result={'scope':'Same-plant identity confirmed by the user on 2026-11-04; existing assembly concordance, not independent topology truth or QV',
            'k':a.k,'main_record_min_bp':1000,'script_sha256':sha(Path(__file__)),'inputs':{},'comparisons':[]}
    def load(path):
        path=path.resolve();result['inputs'][str(path)]=sha(path);return fasta(path)
    for org in ['mitochondrion','plastid']:
        h=load(a.hifi_root/org/'molecules.fasta');reference=words(h,a.k)
        for entry in a.wgs_run:
            label,path=entry.split('=',1);root=Path(path);rep=load(root/(org+'.fasta'))
            candidates=load(root/(org+'.configurations.fasta'));groups=defaultdict(list)
            for rec in candidates:groups[int(rec['configuration'])].append(rec)
            rows=[];exact=[];exact_closed=[];content_groups=[]
            for config,recs in sorted(groups.items()):
                main=[rec for rec in recs if len(rec['seq'])>=1000]
                fit=agreement(reference,words(main,a.k));fit.update(configuration=config,lengths=[len(r['seq']) for r in main],circular=[r['circular'] for r in main])
                rows.append(fit)
                if fit['hifi_only_words']==0 and fit['wgs_only_words']==0:content_groups.append(config)
                if exact_molecule_match(main,h):
                    exact.append(config)
                    if sorted(r['circular'] for r in main)==sorted(r['circular'] for r in h):exact_closed.append(config)
            best=max(rows,key=lambda x:(x['hifi_words_recovered'],x['wgs_words_found_in_hifi']))
            row=dict(run=label,organelle=org,hifi_lengths=[len(r['seq']) for r in h],
                     representative_lengths=[len(r['seq']) for r in rep],
                     representative_all=agreement(reference,words(rep,a.k)),
                     representative_main=agreement(reference,words([rec for rec in rep if len(rec['seq'])>=1000],a.k)),
                     candidates=len(rows),best_sequence_content_candidate=best,
                     identical_word_content_candidates=content_groups,
                     exact_main_molecule_sequence_candidates_up_to_rotation_rc=exact,
                     exact_main_sequence_and_circle_flag_candidates=exact_closed,
                     candidate_rows=rows)
            result['comparisons'].append(row)
            (a.out/'comparison.json').write_text(json.dumps(result,indent=2))
            print(json.dumps({key:value for key,value in row.items() if key!='candidate_rows'}),flush=True)
    for path,h in result['inputs'].items():assert sha(Path(path))==h,'input changed: '+path
    result['input_hashes_stable']=True
    (a.out/'comparison.json').write_text(json.dumps(result,indent=2))


if __name__=='__main__':main()
