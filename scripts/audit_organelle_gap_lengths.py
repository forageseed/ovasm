#!/usr/bin/env python3
"""Reproduce exact assembly gaps and original-read length support from existing PAFs.

This reads existing files only and writes derived audit JSON. It does not patch an
assembly or turn an assembly/read self-consistency test into independent truth.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def fasta(p):
    rows=[];name=None;parts=[]
    for line in p.read_text().splitlines():
        if line.startswith('>'):
            if name is not None:rows.append((name,''.join(parts)))
            name=line[1:].split()[0];parts=[]
        else:parts.append(line.strip().upper())
    if name is not None:rows.append((name,''.join(parts)))
    return rows


def rc(s):
    return s.translate(str.maketrans('ACGT','TGCA'))[::-1]


def paf(line):
    f=line.rstrip().split('\t')
    return f,{tag[:2]:tag[5:] for tag in f[12:]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit-dir',type=Path,required=True)
    args=parser.parse_args();d=args.audit_dir
    names=['plastid-hifi.fasta','plastid-wgs.fasta','pt.paf','pt-hifi-reads.paf']
    hashes={n:sha(d/n) for n in names}
    reference=fasta(d/'plastid-hifi.fasta')[0][1]
    query_name,query=fasta(d/'plastid-wgs.fasta')[0]
    events=[];ranges=[]
    for line in (d/'pt.paf').open():
        f,tags=paf(line)
        if f[0]!=query_name:continue
        qstart,qend=int(f[2]),int(f[3]);ranges.append((qstart,qend))
        segment=query[qstart:qend]
        if f[4]=='-':segment=rc(segment)
        target=int(f[7]);q=0
        for ns,op in re.findall(r'(\d+)([MIDNSHP=X])',tags['cg']):
            n=int(ns)
            assert op in ['=','D'],'assembly comparison also contains non-deletion differences'
            if op=='=':
                assert reference[target:target+n]==segment[q:q+n]
                target+=n;q+=n
            else:
                events.append(dict(hifi_start0=target,hifi_end0=target+n,length=n,
                                   sequence=reference[target:target+n],
                                   query_boundary0=qstart+q if f[4]=='+' else qend-q))
                target+=n
        assert q==len(segment) and target==int(f[8])
    cursor=0
    for lo,hi in sorted(ranges):
        assert lo==cursor,'query contains unexamined or overlapping alignment spans'
        cursor=hi
    assert cursor==len(query)
    missing={i for e in events for i in range(e['hifi_start0'],e['hifi_end0'])}
    shortened=''.join(base for i,base in enumerate(reference) if i not in missing)
    assert len(shortened)==len(query) and (query in shortened*2 or rc(query) in shortened*2)
    bounds=sorted({x for e in events for x in [e['hifi_start0']-300,e['hifi_end0']+300]})
    fits=[{} for e in events];alignment_records=0
    for line in (d/'pt-hifi-reads.paf').open():
        f,tags=paf(line);alignment_records+=1
        if int(f[11])<20 or int(f[9])/int(f[10])<.98:continue
        target=int(f[7]);q=0;maps={}
        for ns,op in re.findall(r'(\d+)([MIDNSHP=X])',tags['cg']):
            n=int(ns)
            if op in ['M','=','X']:
                if op=='=':
                    for b in bounds:
                        if target+20<=b<=target+n-20:maps[b]=q+b-target
                target+=n;q+=n
            elif op in ['D','N']:target+=n
            elif op in ['I','S']:q+=n
        for i,e in enumerate(events):
            lo,hi=e['hifi_start0']-300,e['hifi_end0']+300
            if lo not in maps or hi not in maps:continue
            delta=maps[hi]-maps[lo]-(hi-lo);score=int(f[9]);identifier=f[0]
            if identifier not in fits[i] or score>fits[i][identifier][0]:fits[i][identifier]=(score,delta)
    for e,by_id in zip(events,fits):
        hist=Counter(value[1] for value in by_id.values())
        e.update(eligible_unique_records=len(by_id),
                 hifi_length_within2=sum(n for delta,n in hist.items() if abs(delta)<=2),
                 wgs_short_length_within2=sum(n for delta,n in hist.items() if abs(delta+e['length'])<=2),
                 delta_histogram=hist.most_common(25))
    for name,h in hashes.items():assert sha(d/name)==h,'input changed: '+name
    result=dict(scope='Exact local sequence gaps and same-source original-read span consistency; not independent biological truth',
                hifi_length=len(reference),wgs_main_length=len(query),sum_deleted_bp=len(missing),
                exact_circular_equality_after_deletions=True,alignment_records=alignment_records,
                read_filters=dict(mapq_min=20,identity_min=.98,flank_bp=300,exact_boundary_anchor_bp=20,length_tolerance_bp=2),
                events=events,inputs_sha256=hashes,inputs_stable=True,script_sha256=sha(Path(__file__)))
    (d/'plastid-gap-audit-reproduced.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='events'}))
    for e in events:print(json.dumps({k:v for k,v in e.items() if k not in ['sequence','delta_histogram']}))


if __name__=='__main__':main()
