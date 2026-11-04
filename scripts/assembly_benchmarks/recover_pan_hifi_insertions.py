#!/usr/bin/env python3
"""Experimentally propose HiFi-supported insertions into one selected pan FASTA record.

The source reads may be reused by pan reconstruction. This is not independent validation,
organelle-origin classification, complete consensus, or a physical topology caller.
"""
import argparse
import json
import importlib.util
from pathlib import Path
import subprocess

spec = importlib.util.spec_from_file_location("recover_read_supported_insertions", Path(__file__).with_name("recover_read_supported_insertions.py"))
recovery = importlib.util.module_from_spec(spec)
spec.loader.exec_module(recovery)


def select_record(path, record_id):
    found = []
    name, parts = None, []
    with path.open() as source:
        for line in source:
            if line.startswith(">"):
                if name == record_id:
                    found.append("".join(parts).upper())
                fields = line[1:].split()
                if not fields:
                    raise ValueError("empty FASTA header")
                name, parts = fields[0], []
            elif name is None:
                if line.strip():
                    raise ValueError("sequence before FASTA header")
            elif name == record_id:
                parts.append(line.strip())
        if name == record_id:
            found.append("".join(parts).upper())
    if len(found) != 1 or not found[0] or set(found[0]) - set("ACGTN"):
        raise ValueError("record must occur exactly once and contain nonempty ACGTN sequence")
    return found[0]


def run(args):
    if args.threads < 1 or args.min_reads < 1:
        raise ValueError("threads and min-reads must be positive")
    source_hash = recovery.digest(args.pan_fasta)
    sequence = select_record(args.pan_fasta, args.record)
    read_stat = args.hifi.stat()
    args.out_dir.mkdir(parents=True, exist_ok=False)
    ref = args.out_dir / "selected-reference.fasta"
    ref.write_text(">selected_pan_representative\n" + sequence + "\n")
    paf = args.out_dir / "target-hifi.paf"
    cmd = [args.minimap2, "-x", "map-hifi", "--secondary=yes", "-N", "20", "-c", "-t", str(args.threads), str(ref), str(args.hifi)]
    with paf.open("x") as output, (args.out_dir / "minimap2.log").open("x") as log:
        subprocess.run(cmd, stdout=output, stderr=log, check=True)
    if recovery.digest(args.pan_fasta) != source_hash:
        raise ValueError("pan FASTA changed during mapping")
    end_stat = args.hifi.stat()
    if (read_stat.st_size, read_stat.st_mtime_ns) != (end_stat.st_size, end_stat.st_mtime_ns):
        raise ValueError("HiFi input changed during mapping")
    settings = argparse.Namespace(reference=ref, paf=paf, reads=args.hifi, out=args.out_dir / "candidate", min_reads=args.min_reads, min_mapq=20, min_identity=0.95, min_flank=200, min_length=30, max_length=500, min_modal_fraction=0.7)
    report = recovery.run(settings)
    result = {"scope": "Experimental read-supported insertion candidate for the explicitly selected pan representative; no physical topology or organelle origin is established, and reused target reads are not independent validation.", "pan_fasta": str(args.pan_fasta), "pan_fasta_sha256": source_hash, "selected_record": args.record, "hifi_reads": str(args.hifi), "hifi_size": read_stat.st_size, "hifi_mtime_ns": read_stat.st_mtime_ns, "minimap2_command": cmd, "event_count": len(report["events"]), "candidate_length": report["candidate_length"], "candidate": str(args.out_dir / "candidate.fasta"), "details": str(args.out_dir / "candidate.json")}
    with (args.out_dir / "workflow.json").open("x") as f:
        json.dump(result, f, indent=2)
        f.write("\n")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pan-fasta", type=Path, required=True)
    p.add_argument("--record", required=True, help="exact FASTA record ID, explicitly selected")
    p.add_argument("--hifi", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True, help="new directory; never overwritten")
    p.add_argument("--minimap2", default="minimap2")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--min-reads", type=int, default=1000)
    print(json.dumps(run(p.parse_args()), indent=2))


if __name__ == "__main__":
    main()
