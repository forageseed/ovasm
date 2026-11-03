#!/usr/bin/env python3
"""One-pass nested FASTQ subsampling. Fractions are sampling probabilities, not nuclear depth."""
import argparse
from collections import Counter
from contextlib import ExitStack
from decimal import Decimal
import gzip
import hashlib
import json
from pathlib import Path
import sys


def records(handle):
    """Yield complete raw FASTQ records, supporting wrapped sequence/quality."""
    index = 0
    while True:
        header = handle.readline()
        if not header:
            return
        index += 1
        if not header.startswith(b"@") or not header.rstrip(b"\r\n")[1:]:
            raise ValueError(f"record {index}: invalid FASTQ header")
        raw = [header]
        length = 0
        while True:
            line = handle.readline()
            if not line:
                raise ValueError(f"record {index}: missing + line")
            raw.append(line)
            if line.startswith(b"+"):
                break
            length += len(line.rstrip(b"\r\n"))
        if not length:
            raise ValueError(f"record {index}: empty sequence")
        quality_length = 0
        while quality_length < length:
            line = handle.readline()
            if not line:
                raise ValueError(f"record {index}: truncated qualities")
            raw.append(line)
            quality_length += len(line.rstrip(b"\r\n"))
        if quality_length != length:
            raise ValueError(f"record {index}: quality length differs from sequence")
        yield index, header.rstrip(b"\r\n")[1:].split()[0], length, b"".join(raw)


def fraction(text):
    value = Decimal(text)
    if not value.is_finite() or not 0 < value <= 1:
        raise argparse.ArgumentTypeError("fraction must be finite and in (0, 1]")
    return value


def distribution(counts):
    reads = sum(counts.values())
    bases = sum(k * n for k, n in counts.items())
    return {
        "reads": reads, "bases": bases,
        "mean_length": bases / reads if reads else None,
        "min_length": min(counts) if counts else None,
        "max_length": max(counts) if counts else None,
        "length_histogram": [[k, counts[k]] for k in sorted(counts)],
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reads", action="append", required=True, type=Path)
    p.add_argument("--source-id", action="append",
                   help="Stable source IDs, one per --reads; default input ordinal:basename")
    p.add_argument("--fraction", action="append", required=True, type=fraction)
    p.add_argument("--seed", action="append", required=True, type=int)
    p.add_argument("--out", required=True, type=Path, help="New directory; existing paths refused")
    args = p.parse_args(argv)
    paths = [x.resolve(strict=True) for x in args.reads]
    if len(set(paths)) != len(paths):
        p.error("duplicate input paths are not allowed")
    ids = args.source_id or [f"{i}:{x.name}" for i, x in enumerate(paths)]
    if len(ids) != len(paths) or len(set(ids)) != len(ids):
        p.error("source IDs must be unique and match the input count")
    fractions = sorted(set(args.fraction))
    seeds = sorted(set(args.seed))
    snapshots = [x.stat() for x in paths]
    # No side effects until arguments and inputs have been validated.
    args.out.mkdir(parents=True, exist_ok=False)
    marker = args.out / "INCOMPLETE"
    marker.write_text("Subsampling has not completed successfully. Do not benchmark these outputs.\n")
    manifest = {
        "status": "incomplete", "schema_version": 1,
        "sampling": "nested Bernoulli thresholds on uniform 64-bit BLAKE2b",
        "key": "JSON [seed, source_id, 1-based record_index, read_id_hex], UTF-8",
        "source_id_policy": "Explicit IDs preferred; defaults ordinal:basename. Preserve IDs and record order when relocating inputs.",
        "record_note": "Sequence and quality bytes preserved; missing final newline is added to keep concatenated FASTQ valid.",
        "coverage_note": "Fractions are not nuclear coverage. Actual bases may differ from expected; reads are never truncated.",
        "seeds": seeds, "fractions": [str(f) for f in fractions],
        "sources": [{"path": str(x), "source_id": sid, "size": st.st_size,
                     "mtime_ns": st.st_mtime_ns} for x, sid, st in zip(paths, ids, snapshots)],
        "subsets": [],
    }
    manifest_path = args.out / "manifest.json"

    def save():
        tmp = args.out / "manifest.json.tmp"
        tmp.write_text(json.dumps(manifest, indent=2) + "\n")
        tmp.replace(manifest_path)

    save()
    try:
        with ExitStack() as stack:
            outputs = []
            for seed in seeds:
                for f in fractions:
                    name = f"seed-{seed}.fraction-{format(f, 'f')}.fastq.gz"
                    raw = stack.enter_context((args.out / name).open("xb"))
                    writer = stack.enter_context(gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0))
                    entry = {"seed": seed, "fraction": str(f), "path": name}
                    manifest["subsets"].append(entry)
                    outputs.append((seed, int(f * (1 << 64)), writer, Counter(), hashlib.sha256(), entry))
            for path, sid, before, source in zip(paths, ids, snapshots, manifest["sources"]):
                counts = Counter()
                digest = hashlib.sha256()
                opener = gzip.open if path.suffix.lower() == ".gz" else open
                with opener(path, "rb") as handle:
                    for idx, read_id, length, record in records(handle):
                        counts[length] += 1
                        digest.update(record)
                        hashes = {}
                        for seed in seeds:
                            key = json.dumps([seed, sid, idx, read_id.hex()], separators=(",", ":")).encode()
                            hashes[seed] = int.from_bytes(hashlib.blake2b(key, digest_size=8).digest(), "big")
                        for seed, threshold, writer, hist, out_digest, _ in outputs:
                            if hashes[seed] < threshold:
                                emitted = record if record.endswith(b"\n") else record + b"\n"
                                writer.write(emitted)
                                hist[length] += 1
                                out_digest.update(emitted)
                after = path.stat()
                if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
                    raise RuntimeError(f"input changed while reading: {path}")
                source.update(distribution(counts))
                source["uncompressed_fastq_sha256"] = digest.hexdigest()
            total_bases = sum(x["bases"] for x in manifest["sources"])
            for _, _, _, hist, digest, entry in outputs:
                entry.update(distribution(hist))
                entry["actual_base_fraction"] = entry["bases"] / total_bases if total_bases else None
                entry["uncompressed_fastq_sha256"] = digest.hexdigest()
        for entry in manifest["subsets"]:
            entry["compressed_bytes"] = (args.out / entry["path"]).stat().st_size
        manifest["status"] = "complete"
        save()
        marker.unlink()
    except BaseException as exc:
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        save()
        raise
    return 0


if __name__ == "__main__":
    sys.exit(main())
