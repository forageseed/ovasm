"""Compare completed validation outputs without requiring atlas status.json.
Usage: python3 -I compare_validation.py RUN_ROOT ZENODO_DATA MINIMAP2 OUTPUT_JSON
Pending samples remain explicit. Coordinates are zero-based half-open.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def fasta(path):
    records = {}
    name = None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in records:
                raise ValueError(f"duplicate FASTA id: {name}")
            records[name] = ""
        elif line.strip():
            if name is None:
                raise ValueError(f"sequence before header: {path}")
            records[name] += line.strip()
    return records


def compare(query, target, mm, work):
    windows = work / "windows.fa"
    with windows.open("w") as out:
        for name, seq in fasta(query).items():
            starts = list(range(0, max(len(seq)-10000, 0)+1, 5000))
            if starts[-1]+10000 < len(seq):
                starts.append(len(seq)-10000)
            for start in starts:
                out.write(f">{name}:{start}\n{seq[start:start+10000]}\n")
    proc = subprocess.run([str(mm), "-x", "asm5", "-t", "4", "-N", "50", "-p", "0.5", str(target), str(windows)], capture_output=True, text=True, check=True)
    lengths = {k: len(v) for k, v in fasta(target).items()}
    spans = {k: [] for k in lengths}
    match = aligned = 0
    for line in proc.stdout.splitlines():
        c = line.split("\t")
        if int(c[10]) < 1000:
            continue
        spans[c[5]].append((int(c[7]), int(c[8])))
        if "tp:A:P" in c[12:]:
            match += int(c[9]); aligned += int(c[10])
    result = {}
    for name, size in lengths.items():
        merged = []
        for a,b in sorted(spans[name]):
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(b, merged[-1][1])
            else:
                merged.append([a,b])
        missing = []
        cursor = 0
        for a,b in merged:
            if a > cursor: missing.append([cursor,a])
            cursor = b
        if cursor < size: missing.append([cursor,size])
        result[name] = dict(length=size, covered_bp=sum(b-a for a,b in merged), uncovered_intervals=missing)
    total = sum(lengths.values())
    return dict(coverage=sum(x["covered_bp"] for x in result.values())/total if total else None,
                identity=match/aligned if aligned else None, references=result)


def main():
    root, refs, mm, output = map(Path, sys.argv[1:])
    if output.exists(): raise FileExistsError(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for sample in sorted(root.glob("*/*")):
        if not sample.is_dir(): continue
        summary = sample / "ovasm/summary.json"
        row = dict(variant=sample.parent.name, sample=sample.name, summary=str(summary), status="pending")
        timing = sample / "time.txt"
        timing_text = timing.read_text() if timing.exists() else None
        row.update(command_record=str(timing), command_record_text=timing_text)
        rows.append(row)
        if not summary.exists():
            timing = sample / "time.txt"
            if timing.exists() and "Exit status:" in timing.read_text():
                row["status"] = "finished_without_summary"
                row["command_record"] = str(timing)
            continue
        species = sample.name.split("_")[:2]
        key = species[0][:3].capitalize()+species[1][:4].capitalize()
        hits = [p for p in refs.iterdir() if key in p.name]
        if len(hits) != 1: raise ValueError(f"ambiguous reference: {sample.name}: {hits}")
        failed = json.loads(summary.read_text()).get("failed", [])
        row.update(status="partial" if failed else "complete", failed=failed)
        for org, ref in [("mitochondrion", "MITO.fasta"), ("plastid", "PLTD.fasta")]:
            query = sample / "ovasm" / org / "molecules.fasta"
            target = hits[0] / ref
            if not query.exists() or not target.exists():
                row[org] = dict(status="missing_assembly_or_reference")
                continue
            with tempfile.TemporaryDirectory(dir=output.parent) as tmp:
                row[org] = dict(query=str(query), reference=str(target), query_sha256=hashlib.sha256(query.read_bytes()).hexdigest(), reference_sha256=hashlib.sha256(target.read_bytes()).hexdigest(), oatk=compare(query,target,mm,Path(tmp)), ovasm=compare(target,query,mm,Path(tmp)))
    if not rows: raise ValueError("no sample directories")
    output.write_text(json.dumps(dict(coordinates="zero-based half-open; windowed asm5 alignment, minimum alignment 1 kb", samples=rows),indent=2)+"\n")
    for row in rows:
        print(row["variant"], row["sample"], row["status"], {k:row.get(k,{}).get("oatk",{}).get("coverage") for k in ["mitochondrion","plastid"]})

if __name__ == "__main__": main()
