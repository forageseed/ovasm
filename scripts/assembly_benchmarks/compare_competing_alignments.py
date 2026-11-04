#!/usr/bin/env python3
"""Compare alignments against explicit organelle and nuclear references.

Input must be a single PAF generated with the same parameters against a combined
reference, with secondary alignments and AS:i scores retained. An alignment score
comparison is not a posterior and cannot rule out unrepresented loci or copies.
"""
import argparse
from collections import Counter
import json
import math
from pathlib import Path

def compare(paf, organelle_ids, nuclear_ids, min_fraction=.9, min_identity=.95, margin=10):
    org, nuclear = set(organelle_ids), set(nuclear_ids)
    if not org or not nuclear or org & nuclear:
        raise ValueError("specify nonempty, disjoint organelle and nuclear target IDs")
    if not (math.isfinite(min_fraction) and 0 < min_fraction <= 1 and math.isfinite(min_identity) and 0 < min_identity <= 1 and margin >= 0):
        raise ValueError("invalid thresholds")
    reads = {}
    with paf.open() as file:
        for line in file:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                raise ValueError("malformed PAF")
            name, target = f[0], f[5]
            if target not in org and target not in nuclear:
                continue
            qlen, start, end, matches, block = map(int, (f[1], f[2], f[3], f[9], f[10]))
            if qlen <= 0 or not (0 <= start <= end <= qlen and 0 <= matches <= block):
                raise ValueError("invalid PAF coordinates or counts")
            record = reads.setdefault(name, {"length": qlen, "organelle": None, "nuclear": None})
            if qlen != record["length"]:
                raise ValueError("inconsistent query lengths")
            tags = [x[5:] for x in f[12:] if x.startswith("AS:i:")]
            if len(tags) != 1:
                raise ValueError("exactly one AS:i score required")
            score = int(tags[0])
            if end-start < min_fraction*qlen or not block or matches/block < min_identity:
                continue
            group = "organelle" if target in org else "nuclear"
            current = record[group]
            if current is None or score > current["score"]:
                record[group] = {"score": score, "target": target, "query_start": start, "query_end": end, "matches": matches, "block": block}
    rows = []
    for name, row in sorted(reads.items()):
        o, n = row["organelle"], row["nuclear"]
        if o is None or n is None:
            state = "missing_qualifying_competitor"
            delta = None
        elif abs((n["query_end"] - n["query_start"]) - (o["query_end"] - o["query_start"])) > 0.05 * row["length"]:
            state = "alignment_span_incomparable"
            delta = None
        else:
            delta = n["score"]-o["score"]
            state = "nuclear_score_higher" if delta > margin else "organelle_score_higher" if delta < -margin else "score_ambiguous"
        rows.append({"read": name, "status": state, "nuclear_minus_organelle_score": delta, **row})
    return {"scope": "Comparison of scores against specified references only. Neither taxonomic origin nor molecule type is certified; a missing alignment is not evidence of absence and unrepresented loci are not assessed.",
            "parameters": {"min_query_fraction":min_fraction,"min_identity":min_identity,"score_margin":margin,"max_query_span_difference_fraction":0.05},
            "summary":dict(Counter(x["status"] for x in rows)),"reads":rows}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--paf",type=Path,required=True)
    p.add_argument("--organelle-id",action="append",required=True)
    p.add_argument("--nuclear-id",action="append",required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--min-query-fraction",type=float,default=.9)
    p.add_argument("--min-identity",type=float,default=.95)
    p.add_argument("--score-margin",type=int,default=10)
    args=p.parse_args()
    result=compare(args.paf,args.organelle_id,args.nuclear_id,args.min_query_fraction,args.min_identity,args.score_margin)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result["summary"],indent=2))

if __name__=="__main__":
    main()
