#!/usr/bin/env python3
"""Recover exact original FASTQ records for coordinate-bearing seed crops."""

from argparse import ArgumentParser
from collections import Counter, defaultdict
from dataclasses import dataclass
import gzip
import hashlib
import io
import json
import os
from pathlib import Path

from audit_seed_source import SPAN
from discover_graph_frontier import fastx
from prepare_read_admission import SHA256, LockedFile, digest, file_state, open_input, read_json, require
from validate_source_exclusion import new_output, publish_outputs


PLAN_SCHEMA = "ovasm.seed-crop-originals-plan.v1"
REPORT_SCHEMA = "ovasm.seed-crop-originals.v1"


class DigestReader(io.RawIOBase):
    """Hash the file bytes delivered to the FASTQ reader, including compression."""

    def __init__(self, stream):
        self.stream = stream
        self.hash = hashlib.sha256()
        self.bytes_read = 0

    def readable(self):
        return True

    def readinto(self, buffer):
        size = self.stream.readinto(buffer)
        self.hash.update(memoryview(buffer)[:size])
        self.bytes_read += size
        return size


@dataclass
class Original:
    input_index: int
    record_index: int
    header: bytes
    sequence: bytes
    quality: bytes

    def location(self):
        return self.input_index, self.record_index

    def binding(self):
        return {
            "input_index0": self.input_index, "record_index0": self.record_index,
            "read_id_hex": self.header.hex(), "length": len(self.sequence),
            "sequence_sha256": digest(self.sequence), "quality_sha256": digest(self.quality),
        }

    def fastq(self):
        def line(value):
            return value + (b"\r\n" if value.endswith(b"\r") else b"\n")

        return line(b"@" + self.header) + line(self.sequence) + b"+\n" + line(self.quality)


def read_crops(entries, root):
    require(isinstance(entries, list) and entries, "provide at least one candidate set")
    labels = set()
    proposals = defaultdict(set)
    rows = []
    sources = []
    for entry in entries:
        require(isinstance(entry, dict) and set(entry) == {"id", "path", "sha256"}, "invalid candidate set fields")
        label = entry["id"]
        require(isinstance(label, str) and label.strip() and label not in labels, "invalid candidate set ID")
        labels.add(label)
        source = LockedFile.lock(root, {k: entry[k] for k in ("path", "sha256")})
        sources.append(source)
        with open_input(source) as stream:
            require(stream.peek(1)[:1] in (b"", b">"), "candidate crops must be FASTA")
            for ordinal, (header, sequence, _) in enumerate(fastx(stream, preserve_original=True)):
                name, *fields = header.split()
                spans = [SPAN.fullmatch(field) for field in fields if field.startswith(b"span=")]
                require(len(spans) == 1 and spans[0] is not None, "candidate header requires one span=START-END")
                lo, hi = map(int, spans[0].groups())
                require(0 <= lo < hi and hi - lo == len(sequence), "candidate span disagrees with crop length")
                crop = (name, lo, hi, sequence)
                proposals[name].add(crop)
                rows.append((label, ordinal, crop))
    return proposals, rows, sources


def source_path(entry, root):
    require(isinstance(entry, dict) and {"path", "sha256"} <= entry.keys()
            and entry.keys() <= {"path", "sha256", "id_field"}, "invalid original input fields")
    path, expected, field = entry["path"], entry["sha256"], entry.get("id_field", 1)
    require(isinstance(path, str) and path.strip() and isinstance(expected, str) and SHA256.fullmatch(expected),
            "original input requires a path and lowercase SHA-256")
    require(type(field) is int and field > 0, "original id_field must be a positive integer")
    path = (root / path).resolve(strict=True)
    require(path.is_file(), f"expected a regular original file: {path}")
    return path, expected, field


def collect(stream, proposals, matches, counts, names, in_bounds, index, field):
    require(stream.peek(1)[:1] in (b"", b"@"), "original records must be four-line FASTQ")
    total = 0
    for ordinal, (header, sequence, quality) in enumerate(fastx(stream, preserve_original=True)):
        total += 1
        fields = header.split()
        require(len(fields) >= field, "original header lacks requested identifier field")
        name = fields[field - 1]
        if name not in proposals:
            continue
        names.add(name)
        original = Original(index, ordinal, header, sequence, quality)
        for crop in proposals[name]:
            _, lo, hi, wanted = crop
            if hi > len(sequence):
                continue
            in_bounds.add(crop)
            if sequence[lo:hi] != wanted:
                continue
            counts[crop] += 1
            if counts[crop] == 1:
                matches[crop] = original
            else:
                matches.pop(crop, None)
    return total


def scan_originals(entries, root, proposals):
    require(isinstance(entries, list) and entries, "provide at least one original FASTQ")
    matches, counts, names, in_bounds = {}, Counter(), set(), set()
    files, states, bindings = set(), [], []
    for index, entry in enumerate(entries):
        path, expected, field = source_path(entry, root)
        with path.open("rb", buffering=0) as stream:
            state = file_state(os.fstat(stream.fileno()))
            require(state[:2] not in files, "repeated original input file")
            files.add(state[:2])
            tracker = DigestReader(stream)
            with io.BufferedReader(tracker) as raw:
                if raw.peek(2).startswith(b"\x1f\x8b"):
                    with gzip.GzipFile(fileobj=raw, mode="rb") as parsed:
                        total = collect(parsed, proposals, matches, counts, names, in_bounds, index, field)
                else:
                    total = collect(raw, proposals, matches, counts, names, in_bounds, index, field)
            require(tracker.bytes_read == state[2] and tracker.hash.hexdigest() == expected,
                    f"original SHA-256 or byte count mismatch: {path}")
            require(file_state(os.fstat(stream.fileno())) == state, f"original changed while scanning: {path}")
        states.append((path, state))
        bindings.append({"path": str(path), "sha256": expected, "bytes": state[2],
                         "id_field": field, "records": total})
    for path, state in states:
        require(file_state(path.stat()) == state, f"original changed while scanning: {path}")
    return matches, counts, names, in_bounds, bindings


def prepare(plan_path):
    root = plan_path.resolve().parent
    plan_file = LockedFile.lock(root, {"path": str(plan_path.resolve()), "sha256": digest(plan_path.read_bytes())})
    plan = read_json(plan_file)
    require(isinstance(plan, dict) and set(plan) == {"schema", "scope", "candidate_sets", "raw_reads"}
            and plan["schema"] == PLAN_SCHEMA, "invalid seed crop original plan schema or fields")
    require(isinstance(plan["scope"], str) and plan["scope"].strip(), "provide a scope")
    proposals, rows, sources = read_crops(plan["candidate_sets"], root)
    matches, counts, names, in_bounds, raw = scan_originals(plan["raw_reads"], root, proposals)
    for source in (plan_file, *sources):
        source.check()
    originals = {match.location(): match for match in matches.values()}
    ordered = sorted(originals)
    original_rows = {location: originals[location].binding() | {"output_record_index0": i}
                     for i, location in enumerate(ordered)}
    output = b"".join(originals[location].fastq() for location in ordered)
    decisions = []
    for label, ordinal, crop in rows:
        name, lo, hi, sequence = crop
        original = matches.get(crop)
        if original is not None:
            state = "bound"
        elif counts[crop] > 1:
            state = "multiple_exact_originals"
        elif crop in in_bounds:
            state = "crop_bytes_mismatch"
        elif name in names:
            state = "span_outside_originals"
        else:
            state = "identifier_not_observed"
        decisions.append({
            "candidate_set": label, "candidate_record_index0": ordinal, "candidate_id_hex": name.hex(),
            "span": [lo, hi], "crop_sha256": digest(sequence), "exact_original_matches": counts[crop],
            "state": state, "original_output_record_index0": (original_rows[original.location()]["output_record_index0"]
                                                              if original else None),
        })
    report = {
        "schema": REPORT_SCHEMA, "scope": plan["scope"], "plan": plan_file.binding(),
        "coordinate_system": "zero-based half-open spans and zero-based record indices",
        "raw_reads": raw, "candidate_sets": [
            source.binding() | {"id": entry["id"]} for source, entry in zip(sources, plan["candidate_sets"])
        ],
        "candidate_counts": Counter(row["state"] for row in decisions), "candidates": decisions,
        "originals": list(original_rows.values()), "output_records": len(ordered),
    }
    return report, output


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path)
    parser.add_argument("--out", type=Path, required=True, help="new original identity report")
    parser.add_argument("--out-reads", type=Path, required=True, help="new original FASTQ when reads bind")
    args = parser.parse_args()
    try:
        out, reads = new_output(args.out), new_output(args.out_reads)
        require(out != reads, "report and reads need different paths")
        report, data = prepare(args.plan)
        report["output_fastq"] = {"path": str(reads), "sha256": digest(data)} if data else None
        outputs = [(out, (json.dumps(report, indent=2, sort_keys=True) + "\n").encode())]
        if data:
            outputs.insert(0, (reads, data))
        publish_outputs(outputs)
    except (OSError, EOFError, ValueError, TypeError, KeyError) as error:
        parser.exit(1, f"seed crop original binding: {error}\n")
    print(f'{len(report["candidates"])} candidate crops, {report["output_records"]} bound originals')


if __name__ == "__main__":
    main()
