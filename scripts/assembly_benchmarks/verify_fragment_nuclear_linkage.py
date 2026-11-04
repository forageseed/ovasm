#!/usr/bin/env python3
"""Verify exact fragment + two physical 500bp flanks in original read records.

Consumes existing context alignments. A positive result establishes a nuclear-
connected copy; it is not an exhaustive absence test for organellar copies.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def fasta(p):
    rows={};name=None;parts=[]
    for line in p.read_text().splitlines():
        if line.startswith('>'):
            if name is not None:rows[name]=''.join(parts)
            name=line[1:].split()[0];parts=[]
        else:parts.append(line.strip().upper())
    if name is not None:rows[name]=''.join(parts)
    return rows


def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def words(s):return {min(s[i:i+31],rc(s[i:i+31])) for i in range(len(s)-30)}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['raw-fastq','contexts','paf','nuclear','mitochondrial','fragment','out']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();inputs={n:getattr(a,n) for n in ['raw_fastq','contexts','paf','nuclear','mitochondrial','fragment']}
    hashes={n:sha(path) for n,path in inputs.items()}
    raw={}
    with a.raw_fastq.open() as handle:
        while True:
            header=handle.readline()
            if not header:break
            s=handle.readline().strip();plus=handle.readline();quality=handle.readline().strip()
            assert header.startswith('@') and plus.startswith('+') and len(s)==len(quality)
            key=header[1:].split()[0]
            assert key not in raw or raw[key]==s
            raw[key]=s
    contexts=fasta(a.contexts);nuclear=fasta(a.nuclear);fragment=next(iter(fasta(a.fragment).values()))
    mt_words=set()
    for s in fasta(a.mitochondrial).values():mt_words|=words(s+s[:30])
    best={}
    for line in a.paf.open():
        f=line.rstrip().split('\t');tags={x[:2]:x[5:] for x in f[12:]}
        if f[5] not in nuclear or int(f[2])!=0 or int(f[3])!=1453:continue
        identifier,position,strand=f[0].rsplit('|',2);position=int(position)
        if identifier not in raw:continue
        source=raw[identifier];context=contexts[f[0]]
        assert position>=500 and position+len(fragment)+500<=len(source)
        assert context==source[position-500:position+len(fragment)+500]
        assert source[position:position+len(fragment)]==(fragment if strand=='+' else rc(fragment))
        q=context if f[4]=='+' else rc(context);reference=nuclear[f[5]];t=int(f[7]);qi=0;edits=0
        for ns,op in re.findall(r'(\d+)([MIDNSHP=X])',tags['cg']):
            n=int(ns)
            if op=='=':
                assert q[qi:qi+n]==reference[t:t+n];qi+=n;t+=n
            elif op=='X':qi+=n;t+=n;edits+=n
            elif op=='I':qi+=n;edits+=n
            elif op=='D':t+=n;edits+=n
            else:raise AssertionError('unexpected context CIGAR operation')
        assert qi==len(q) and t==int(f[8]) and edits==int(tags['NM'])
        if int(f[9])/int(f[10])<.98:continue
        center=q[500:953];start=reference.find(center)
        if start<550:continue
        left_specific=words(reference[start-550:start])-mt_words
        right_specific=words(reference[start+453:start+453+550])-mt_words
        left=len(words(q[:500])&left_specific);right=len(words(q[953:])&right_specific)
        if min(left,right)<100:continue
        row=dict(id=identifier,movie=identifier.split('/')[0],original_read_bp=len(source),
                 exact_fragment_bp=453,physical_flank_bp_each=500,context_bp=len(q),
                 reference=f[5],target_span0=[int(f[7]),int(f[8])],identity=int(f[9])/int(f[10]),
                 edits=edits,left_nuclear_specific_words=left,right_nuclear_specific_words=right,
                 source_sequence_sha256=hashlib.sha256(source.encode()).hexdigest())
        if identifier not in best or row['identity']>best[identifier]['identity']:best[identifier]=row
    assert len(best)>=3,'multiple independent original molecule criterion not met'
    for name,path in inputs.items():assert sha(path)==hashes[name],'input changed: '+name
    rows=list(best.values());result=dict(verified=True,conclusion='nuclear-connected fragment containing mitochondrial homology; NUMT classification supported',
        independent_original_records=len(rows),distinct_movies=len({r['movie'] for r in rows}),
        excludes_all_possible_organellar_copies=False,rows=rows,inputs_sha256=hashes,
        script_sha256=sha(Path(__file__)))
    a.out.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':main()
