#!/usr/bin/env python3
"""Produce an explicitly experimental read-supported insertion candidate from PAF and FASTQ.

The reference must have one linearized record. PAF must be generated against this exact
record, for example minimap2 -x map-hifi --secondary=yes -N 20 -c reference.fa reads.fq.
This is not a topology caller, a probabilistic posterior, or an independent validation.
"""
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import re

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def reference(path):
    records = []
    for part in path.read_text().split(">"):
        lines = part.splitlines()
        if lines:
            records.append((lines[0].split()[0], "".join(lines[1:]).upper()))
    if len(records) != 1 or not records[0][1] or set(records[0][1]) - set("ACGTN"):
        raise ValueError("reference must contain one nonempty ACGTN FASTA record")
    return records[0]


def run(args):
    if args.min_reads < 1 or not 0 <= args.min_mapq <= 255 or args.min_flank < 0 or args.min_length < 1 or args.max_length < args.min_length or not 0 <= args.min_identity <= 1 or not 0 <= args.min_modal_fraction <= 1:
        raise ValueError("invalid support, length, flank, identity, or modal fraction threshold")
    fasta = Path(str(args.out) + ".fasta")
    report = Path(str(args.out) + ".json")
    if fasta.exists() or report.exists():
        raise FileExistsError("refusing to overwrite output")
    reference_hash, paf_hash = digest(args.reference), digest(args.paf)
    target, ref = reference(args.reference)
    input_stat = args.reads.stat()
    events = defaultdict(dict)
    paf_records = 0
    with args.paf.open() as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                raise ValueError("malformed PAF")
            paf_records += 1
            q, qlen, qs, qe, strand = f[0], int(f[1]), int(f[2]), int(f[3]), f[4]
            if f[5] != target or int(f[6]) != len(ref):
                raise ValueError("PAF target name or length does not match reference")
            if strand not in ("+", "-"):
                raise ValueError("invalid PAF strand")
            if int(f[11]) < args.min_mapq or int(f[9]) / max(1, int(f[10])) < args.min_identity or "tp:A:P" not in f[12:]:
                continue
            ts, te = int(f[7]), int(f[8])
            tpos, qpos = ts, qs if strand == "+" else qlen - qe
            cigar = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
            if cigar is None:
                raise ValueError("PAF must include cg:Z: CIGAR (-c)")
            operations = CIGAR.findall(cigar)
            if "".join(n + op for n, op in operations) != cigar:
                raise ValueError("invalid CIGAR")
            for num, op in operations:
                n = int(num)
                if op == "I" and args.min_length <= n <= args.max_length and tpos - ts >= args.min_flank and te - tpos >= args.min_flank:
                    entries = events[(tpos, n)]
                    new = (qpos, strand, qlen)
                    if q not in entries:
                        entries[q] = new
                    elif entries[q] != new:
                        entries[q] = None  # exclude conflicting mappings for this event
                if op in "MI=X":
                    qpos += n
                if op in "MDN=X":
                    tpos += n
                if op in "SHP":
                    raise ValueError("clipped/padded PAF CIGAR is unsupported")
    positions = defaultdict(list)
    ambiguous_event_read_pairs = 0
    for (pos, length), names in events.items():
        ambiguous_event_read_pairs += sum(v is None for v in names.values())
        names = {name: v for name, v in names.items() if v is not None}
        if len(names) >= args.min_reads:
            positions[pos].append((len(names), length, names))
    chosen = {p: max(v, key=lambda x: (x[0], -x[1])) for p, v in positions.items()}
    needed = defaultdict(list)
    for pos, (_, length, names) in chosen.items():
        for name, (qpos, strand, qlen) in names.items():
            needed[name].append((pos, length, qpos, strand, qlen))
    seqcounts = defaultdict(Counter)
    seen = set()
    scanned = 0
    opener = gzip.open if args.reads.suffix == ".gz" else open
    with opener(args.reads, "rt") as fh:
        while True:
            head = fh.readline()
            if not head:
                break
            seq, plus, qual = fh.readline().strip(), fh.readline(), fh.readline().strip()
            if not head.startswith("@") or not plus.startswith("+") or len(seq) != len(qual):
                raise ValueError("invalid or truncated FASTQ")
            scanned += 1
            name = head[1:].split()[0]
            if name not in needed:
                continue
            if name in seen:
                raise ValueError("duplicate selected read name: " + name)
            seen.add(name)
            reverse = None
            for pos, length, qpos, strand, qlen in needed[name]:
                if len(seq) != qlen:
                    raise ValueError("PAF query length differs from FASTQ")
                if strand == "-":
                    if reverse is None:
                        reverse = seq.translate(COMPLEMENT)[::-1]
                    oriented = reverse
                else:
                    oriented = seq
                insert = oriented[qpos:qpos + length].upper()
                if len(insert) != length:
                    raise ValueError("PAF query coordinates exceed FASTQ")
                seqcounts[(pos, length)][insert] += 1
    if not needed.keys() <= seen:
        raise ValueError("some selected PAF query names are missing from FASTQ")
    after = args.reads.stat()
    if (input_stat.st_size, input_stat.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError("FASTQ changed during scan")
    if digest(args.reference) != reference_hash or digest(args.paf) != paf_hash:
        raise ValueError("reference or PAF changed during scan")
    last = 0
    pieces, details, rejected = [], [], []
    for pos in sorted(chosen):
        support, length, _ = chosen[pos]
        sequence, modal_count = seqcounts[(pos, length)].most_common(1)[0]
        if modal_count / support < args.min_modal_fraction:
            rejected.append({"position": pos, "length": length, "reads": support, "modal_reads": modal_count})
            continue
        pieces.extend((ref[last:pos], sequence))
        last = pos
        details.append({"position": pos, "length": length, "read_support": support, "modal_sequence_reads": modal_count, "inserted_sequence": sequence, "alternate_eligible_lengths": [{"length": n, "read_support": c} for c, n, _ in positions[pos] if n != length]})
    pieces.append(ref[last:])
    sequence = "".join(pieces)
    result = {"scope": "Experimental read-supported insertion candidate; no physical topology, organelle origin, unbiased abundance, or calibrated posterior is established. Repeats and PAF mapping can move or collapse events; same-read validation is not independent. Substitutions and deletions are not corrected.", "reference": str(args.reference), "reference_sha256": reference_hash, "paf": str(args.paf), "paf_sha256": paf_hash, "reads": str(args.reads), "reads_size": input_stat.st_size, "reads_mtime_ns": input_stat.st_mtime_ns, "paf_records": paf_records, "read_records": scanned, "ambiguous_event_read_pairs_excluded": ambiguous_event_read_pairs, "thresholds": {"min_reads": args.min_reads, "min_mapq": args.min_mapq, "min_identity": args.min_identity, "min_flank": args.min_flank, "min_length": args.min_length, "max_length": args.max_length, "min_modal_fraction": args.min_modal_fraction}, "reference_length": len(ref), "candidate_length": len(sequence), "candidate_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(), "events": details, "rejected_low_modal_fraction": rejected}
    with fasta.open("x") as f:
        f.write(">experimental_read_supported_candidate topology=unresolved\n")
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + "\n")
    with report.open("x") as f:
        json.dump(result, f, indent=2)
        f.write("\n")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--paf", type=Path, required=True)
    p.add_argument("--reads", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True, help="new output prefix; FASTA and JSON must not exist")
    p.add_argument("--min-reads", type=int, default=1000)
    p.add_argument("--min-mapq", type=int, default=20)
    p.add_argument("--min-identity", type=float, default=0.95)
    p.add_argument("--min-flank", type=int, default=200)
    p.add_argument("--min-length", type=int, default=30)
    p.add_argument("--max-length", type=int, default=500)
    p.add_argument("--min-modal-fraction", type=float, default=0.7)
    a = p.parse_args()
    result = run(a)
    print(json.dumps({"events": len(result["events"]), "candidate_length": result["candidate_length"], "scope": result["scope"]}))


if __name__ == "__main__":
    main()
