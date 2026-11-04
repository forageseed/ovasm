#!/usr/bin/env python3
"""Independently enumerate minimum explanatory families from a full candidate report.

No native solver scoring or set-cover implementation is imported. Proof applies only
to the supplied finite candidates, never to physical topology or all graph walks.
"""
import argparse
from collections import defaultdict
import itertools
import json
from pathlib import Path


def flip(s):
    return s[:-1] + ('-' if s[-1] == '+' else '+')


def edge(a,b):
    return min((a,b),(flip(b),flip(a)))


def contains(path,circular,query):
    if not path or not query or (not circular and len(query)>len(path)):
        return False
    for q in (query,[flip(x) for x in reversed(query)]):
        for i in range(len(path) if circular else len(path)-len(q)+1):
            if all(path[(i+j)%len(path)] == x for j,x in enumerate(q)):
                return True
    return False


def verify(report,evidence,max_combinations=1000000):
    strict=evidence['thread_evidence']
    repeats=set(evidence.get('repeat_segments',[]))
    results=[]
    for component in report['components']:
        family=component.get('read_path_cover')
        if not family:
            continue
        covers=[component['best']]+component['alternatives']
        assert len(covers)==component['candidates'],'full candidate report required'
        anchors=set(component['anchors']);nodes=set(anchors)
        adjacency=defaultdict(set)
        for link in strict['links']:
            a,b=link['from'][:-1],link['to'][:-1]
            adjacency[a].add(b);adjacency[b].add(a)
        queue=list(nodes)
        while queue:
            for b in adjacency[queue.pop()]:
                if b not in nodes and (b in repeats or b in anchors):
                    nodes.add(b);queue.append(b)
        queries=[p['path'].split() for p in strict['paths'] if p['reads']>=2
                 and any(x[:-1] in anchors and x[:-1] not in repeats for x in p['path'].split())]
        links=[edge(p['from'],p['to']) for p in strict['links'] if p['records']>=2
               and p['from'][:-1] in nodes and p['to'][:-1] in nodes]
        masks=[]
        for cover in covers:
            walks=[(m['path'].split(),m['circular']) for m in cover['molecules']]
            used=set()
            for path,circ in walks:
                used.update(edge(a,b) for a,b in zip(path,path[1:]))
                if circ and path:used.add(edge(path[-1],path[0]))
            mask=0
            for i,q in enumerate(queries):
                if any(contains(p,c,q) for p,c in walks):mask|=1<<i
            for i,e in enumerate(links,len(queries)):
                if e in used:mask|=1<<i
            masks.append(mask)
        target=0
        for m in masks:target|=m
        tested=0;minima=[];exhausted=False
        for size in range(len(family['candidate_indices'])+1):
            for indices in itertools.combinations(range(len(masks)),size):
                tested+=1
                if tested>max_combinations:
                    exhausted=True;break
                union=0
                for i in indices:union|=masks[i]
                if union==target:minima.append(indices)
            if exhausted or minima:break
        if not exhausted:
            assert minima,'native family did not cover representable evidence'
            assert len(minima[0])==len(family['candidate_indices']),'minimum family size differs'
            members=sorted(set(itertools.chain.from_iterable(minima)))
            if family['equivalent_families_complete']:
                assert len(minima)==family['equivalent_minimum_families'],'equivalent family count differs'
                assert members==family['equivalent_candidate_indices'],'ambiguity envelope differs'
        results.append(dict(anchors=sorted(anchors),candidates=len(covers),combinations_tested=tested,
                            independently_complete=not exhausted,minimum_size=len(minima[0]) if minima else None,
                            equivalent_families=len(minima),supported_paths=len(queries),supported_local_links=len(links)))
    return results


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report',type=Path,required=True)
    p.add_argument('--evidence',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    result=verify(json.loads(a.report.read_text()),json.loads(a.evidence.read_text()))
    a.out.write_text(json.dumps(result,indent=2))
    print(json.dumps(result))


if __name__=='__main__':
    main()
