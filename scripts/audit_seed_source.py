#!/usr/bin/env python3
"""Audit native protein seed crops against locked local source contexts."""
from argparse import ArgumentParser
from collections import defaultdict
from contextlib import contextmanager
from hashlib import file_digest, sha256
from pathlib import Path
import gzip
import json
import os
import re
import shutil
import struct
import subprocess
import tempfile
from threading import Thread


MANIFEST = "ovasm.seed-source-audit.v1"
COMPLEMENT = bytes.maketrans(b"ACGT", b"TGCA")
SPAN = re.compile(rb"span=(\d+)-(\d+)")
ROLES = {"target", "competing"}
REQUIRED_SOURCE_ROLES = {"mitochondrion", "plastid", "nuclear"}
REFERENCE_SOURCE_ROLES = REQUIRED_SOURCE_ROLES | {"unplaced"}
EXCLUDE_COMPETING = "exclude-competing"
TARGET_ONLY = "target-supported-only"
PLACEMENT_MAGIC = b"OVXPOS1"
PLACEMENT_END = 0xffffffff
PLACEMENT_HIT = re.compile(rb"H\t([0-9]+)\t([0-9]+)\t([0-9]+)\n")


def require(condition, message):
    if not condition:
        raise ValueError(message)


@contextmanager
def locked_file(reference, base):
    path = (base / reference["path"]).resolve(strict=True)
    expected = reference["sha256"]
    with path.open("rb") as source:
        stat = os.fstat(source.fileno())
        require(file_digest(source, "sha256").hexdigest() == expected, f"SHA-256 mismatch: {path}")
        source.seek(0)
        yield path, source
        final = os.fstat(source.fileno())
        require((stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns) ==
                (final.st_size, final.st_mtime_ns, final.st_ctime_ns), f"input changed: {path}")


def fasta_records(source):
    header = None
    pieces = []
    for line in source:
        if line.startswith(b">"):
            if header is not None:
                yield header, pieces
            header, pieces = line, []
        else:
            require(header is not None, "FASTA has content before its first header")
            pieces.append(line)
    if header is not None:
        yield header, pieces


def unpack_fasta(header, lines):
    fields = header[1:].split()
    require(fields, "empty FASTA identifier")
    return fields, b"".join(line.strip() for line in lines).upper()


def raw_records(source, id_index):
    prefix = source.read(2)
    source.seek(0)
    if prefix == b"\x1f\x8b":
        with gzip.GzipFile(fileobj=source) as stream:
            yield from fastx_records(stream, id_index)
    else:
        yield from fastx_records(source, id_index)


def fastx_records(source, id_index):
    prefix = source.peek(1)[:1]
    if prefix == b">":
        for header, lines in fasta_records(source):
            fields, sequence = unpack_fasta(header, lines)
            require(len(fields) > id_index, "original header lacks requested identifier field")
            yield fields[id_index], sequence
    elif prefix == b"@":
        for header in source:
            sequence, plus, quality = (source.readline().rstrip(b"\r\n") for _ in range(3))
            fields = header[1:].split()
            require(header.startswith(b"@") and len(fields) > id_index and sequence and plus.startswith(b"+")
                    and len(quality) == len(sequence), "invalid four-line FASTQ record")
            require(len(fields) > id_index, "original header lacks requested identifier field")
            yield fields[id_index], sequence.upper()
    else:
        require(not prefix, "expected FASTA or four-line FASTQ")


def placements(sequence, window):
    """Return up to two positions, enough to establish ambiguity."""
    found = []
    for strand, query in (("+", window), ("-", window.translate(COMPLEMENT)[::-1])):
        start = sequence.find(query)
        while start >= 0:
            found.append((start, start + len(query), strand))
            if len(found) == 2:
                return found
            start = sequence.find(query, start + 1)
    return found


def write_placement_records(stream, sequences, errors):
    try:
        with stream:
            stream.write(PLACEMENT_MAGIC + b"\n")
            for sequence in sequences:
                require(len(sequence) < PLACEMENT_END, "exact placement record exceeds the frame limit")
                stream.write(struct.pack("<I", len(sequence)))
                stream.write(sequence)
            stream.write(struct.pack("<I", PLACEMENT_END))
    except (OSError, ValueError) as error:
        errors.append(error)


def read_placement_hits(stream, sequences, patterns, scope):
    header = PLACEMENT_MAGIC + b"\t" + scope.encode("ascii") + b"\n"
    require(stream.readline(128) == header, "invalid exact placement response header")
    terminal = f"E\t{len(sequences)}\t{len(patterns)}\n".encode("ascii")
    hits = defaultdict(list)
    previous_record = -1
    while line := stream.readline(128):
        if line == terminal:
            require(not stream.read(1), "trailing exact placement response")
            return hits
        match = PLACEMENT_HIT.fullmatch(line)
        require(match is not None, "invalid exact placement hit or terminal summary")
        record, pattern, start = map(int, match.groups())
        require(previous_record <= record < len(sequences) and pattern < len(patterns),
                "invalid exact placement record or pattern index")
        sequence, literal = sequences[record], patterns[pattern]
        require(start + len(literal) <= len(sequence) and sequence.startswith(literal, start),
                "invalid exact placement coordinates or literal")
        previous_record = record
        key = (None if scope == "panel" else record, pattern)
        positions = hits[key]
        require(len(positions) < 2 and (not positions or positions[-1] < (record, start)),
                "duplicate, unordered or excess exact placement hits")
        positions.append((record, start))
    raise ValueError("missing exact placement terminal summary")


def batch_placement_lookup(binary, scope, sequences, words):
    def oriented(word):
        return (("+", word), ("-", word.translate(COMPLEMENT)[::-1]))

    oriented_words = {word: oriented(word) for word in words}
    patterns = tuple(dict.fromkeys(literal for pair in oriented_words.values() for _, literal in pair))
    pattern_ids = {literal: index for index, literal in enumerate(patterns)}
    with tempfile.TemporaryDirectory(prefix="ovasm-exact-placements-") as tmp:
        path = Path(tmp) / "patterns.txt"
        with path.open("wb") as output:
            for literal in patterns:
                output.write(literal + b"\n")
        try:
            with subprocess.Popen(
                [binary, "find", "exact-placements", "--patterns", path, "--scope", scope],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            ) as process:
                errors = []
                writer = Thread(target=write_placement_records, args=(process.stdin, sequences, errors))
                writer.start()
                try:
                    hits = read_placement_hits(process.stdout, sequences, patterns, scope)
                    status = process.wait()
                    require(status == 0, f"exact placement backend exit status {status}")
                finally:
                    if process.poll() is None:
                        process.kill()
                    writer.join()
                require(not errors, f"cannot deliver exact placement input: {errors}")
        except (OSError, ValueError) as error:
            raise ValueError(f"OVASM exact placement backend failed ({scope}, {binary}): {error}") from error

    def lookup(word, record=None):
        matches = [(record_index, start, start + len(word), strand)
                   for strand, literal in oriented_words[word]
                   for record_index, start in hits.get((record, pattern_ids[literal]), ())]
        # Reference order takes precedence across records; the legacy search
        # emits all forward positions before reverse positions in each record.
        return sorted(matches, key=lambda hit: (hit[0], hit[3] == "-", hit[1]))[:2]

    return lookup


def legacy_panel_placements(panel, word):
    matches = []
    for ref, (_, records) in panel.items():
        for contig, sequence in records.items():
            matches.extend((ref, contig, *position) for position in placements(sequence, word))
            if len(matches) >= 2:
                return matches[:2]
    return matches


def declared_reference_scope(manifest, target):
    required = manifest.get("required_reference_roles")
    require(isinstance(required, list) and
            all(isinstance(role, str) and role in REQUIRED_SOURCE_ROLES for role in required)
            and len(set(required)) == len(required),
            "required_reference_roles must list distinct mitochondrion, plastid or nuclear roles")
    require({target, "nuclear"} <= set(required),
            "required_reference_roles must include the target and nuclear roles")
    require(target != "plastid" or "mitochondrion" in required,
            "a plastid target also requires the mitochondrial reference role")
    references = manifest.get("references")
    require(isinstance(references, list), "references must be a list")
    declared = []
    reference_ids = set()
    for entry in references:
        require(isinstance(entry, dict), "invalid reference entry")
        reference_id = entry.get("id")
        require(isinstance(reference_id, str) and reference_id and reference_id not in reference_ids,
                "invalid or repeated reference identifier")
        reference_ids.add(reference_id)
        source_role = entry.get("source_role")
        require(isinstance(source_role, str) and source_role in REFERENCE_SOURCE_ROLES,
                "each reference requires a declared mitochondrion, plastid, nuclear or unplaced source_role")
        legacy_role = "target" if source_role == target else "competing"
        require(entry.get("role") == legacy_role, "reference source_role disagrees with target/competing role")
        declared.append({"id": reference_id, "source_role": source_role})
    provided = {entry["source_role"] for entry in declared}
    missing = sorted(set(required) - provided)
    return {
        "policy": "require_declared_reference_roles",
        "reference_role_authority": "user_declared",
        "scope": "supplied_named_reference_records",
        "required_roles": sorted(required),
        "provided_roles": sorted(provided),
        "missing_required_roles": missing,
        "declared_references": declared,
        "source_windows_evaluated": not missing,
        "raw_records_evaluated": not missing,
        "physical_origin": "unresolved",
    }


def load_panel(manifest, base):
    panel = {}
    for entry in manifest["references"]:
        reference_id, role = entry["id"], entry["role"]
        require(isinstance(reference_id, str) and reference_id and reference_id not in panel
                and role in ROLES, "invalid or repeated reference identifier or role")
        with locked_file(entry, base) as (_, source):
            records = {}
            for header, lines in fasta_records(source):
                fields, sequence = unpack_fasta(header, lines)
                record = fields[0].decode("ascii")
                require(sequence and record not in records, "empty or repeated reference record")
                records[record] = sequence
        require(records, "empty reference FASTA")
        panel[reference_id] = (role, records)
    require({role for role, _ in panel.values()} == ROLES,
            "provide both target and competing references")
    return panel


def load_windows(manifest, panel, base, placement_binary=None):
    windows = {}
    for window in manifest["windows"]:
        name, reference, record = window["id"], window["reference"], window["record"]
        require(isinstance(name, str) and name and name not in windows,
                "invalid or repeated window identifier")
        role, records = panel[reference]
        start, end = window["interval"]
        boundary = window["boundary"]
        require(all(type(value) is int for value in (start, end, boundary))
                and 0 <= start < boundary < end <= len(records[record]), f"invalid window: {name}")
        word = records[record][start:end]
        require(not set(word).difference(b"ACGT"), f"ambiguous bases in source window: {name}")
        with locked_file(window["provenance"], base):
            pass
        windows[name] = (role, word, {
            "reference": reference, "record": record, "interval": [start, end],
            "boundary": boundary, "role": role, "sequence_sha256": sha256(word).hexdigest(),
            "provenance_sha256": window["provenance"]["sha256"],
        })
    require(windows, "provide at least one source window")
    if placement_binary is not None:
        entries = [(ref, record, sequence) for ref, (_, records) in panel.items() for record, sequence in records.items()]
        locate = batch_placement_lookup(placement_binary, "panel", [entry[2] for entry in entries],
                                        (word for _, word, _ in windows.values()))
    for _, word, report in windows.values():
        if placement_binary is None:
            matches = legacy_panel_placements(panel, word)
        else:
            matches = [(entries[i][0], entries[i][1], start, end, strand)
                       for i, start, end, strand in locate(word)]
        start, end = report["interval"]
        report["unique_in_named_panel"] = matches == [(report["reference"], report["record"], start, end, "+")]
        report["observed_placements"] = matches
    return windows


def load_seeds(reference, base):
    seeds = []
    with locked_file(reference, base) as (_, source):
        for header, lines in fasta_records(source):
            fields, sequence = unpack_fasta(header, lines)
            spans = [SPAN.fullmatch(field) for field in fields[1:] if field.startswith(b"span=")]
            span = tuple(map(int, spans[0].groups())) if len(spans) == 1 and spans[0] else None
            seeds.append((fields[0], span, sequence, header + b"".join(lines)))
    require(seeds, "empty seed FASTA")
    return seeds


def load_originals(references, base, wanted):
    require(isinstance(references, list) and references, "provide at least one original read file")
    originals = {}
    for index, reference in enumerate(references):
        with locked_file(reference, base) as (_, source):
            field = reference.get("id_field", 1)
            require(type(field) is int and field > 0, "original id_field must be a positive integer")
            for name, sequence in raw_records(source, field - 1):
                if name in wanted:
                    # A repeated record identifier prevents unique original-coordinate binding.
                    originals[name] = None if name in originals else (index, sequence)
    return originals


def load_strict_originals(references, base, seeds):
    """Join exact proposals without storing every original with a repeated name."""
    require(isinstance(references, list) and references, "provide at least one original read file")
    wanted = {name for name, _, _, _ in seeds}
    proposed = defaultdict(set)
    for name, span, seed, _ in seeds:
        if span is not None and span[0] < span[1] and span[1] - span[0] == len(seed):
            proposed[name].add(span)
    keys = {(name, span, seed) for name, span, seed, _ in seeds}
    originals, observed_names, files = {}, set(), set()
    for index, reference in enumerate(references):
        with locked_file(reference, base) as (_, source):
            stat = os.fstat(source.fileno())
            identity = stat.st_dev, stat.st_ino
            require(identity not in files, "strict source audit repeats an original input file")
            files.add(identity)
            field = reference.get("id_field", 1)
            require(type(field) is int and field > 0, "original id_field must be a positive integer")
            for ordinal, (name, sequence) in enumerate(raw_records(source, field - 1)):
                if name not in wanted:
                    continue
                observed_names.add(name)
                for span in proposed[name]:
                    if span[1] > len(sequence):
                        continue
                    key = (name, span, sequence[slice(*span)])
                    if key in keys:
                        originals[key] = None if key in originals else (index, ordinal, sequence)
    return originals, observed_names


def original_placement_lookup(binary, originals, windows):
    sequences = {id(binding[-1]): binding[-1] for binding in originals.values() if binding is not None}
    indices = {identity: index for index, identity in enumerate(sequences)}
    locate = batch_placement_lookup(binary, "records", list(sequences.values()),
                                    (word for _, word, _ in windows.values()))

    def lookup(sequence, word):
        return [(start, end, strand) for _, start, end, strand in locate(word, indices[id(sequence)])]

    return lookup


def unresolved_seed(name, span, seed, reason):
    return {
        "read_id": name.decode("ascii", errors="backslashreplace"),
        "read_id_hex": name.hex(), "crop": span, "seed_sha256": sha256(seed).hexdigest(),
        "state": "unresolved", "reason": reason, "observations": [],
    }


def assess_seed(name, span, seed, originals, windows, strict_names=None, locate=placements):
    row = unresolved_seed(name, span, seed, "original_missing")
    if strict_names is None:
        if name not in originals:
            return row
        binding = originals[name]
        if binding is None:
            row["reason"] = "original_identifier_repeated"
            return row
        index, original = binding
    else:
        if name not in strict_names:
            return row
        key = name, span, seed
        if key not in originals:
            row["reason"] = "crop_does_not_bind_to_original"
            return row
        binding = originals[key]
        if binding is None:
            row["reason"] = "crop_matches_multiple_original_records"
            return row
        index, ordinal, original = binding
        row["original_record_index"] = ordinal
    row["original_input_index"] = index
    row["original_sha256"] = sha256(original).hexdigest()
    if span is None or not 0 <= span[0] < span[1] <= len(original) or original[slice(*span)] != seed:
        row["reason"] = "crop_does_not_bind_to_original"
        return row

    observed_roles, crop_roles = set(), set()
    for window, (role, word, reference) in windows.items():
        positions = locate(original, word)
        if not positions:
            continue
        observed_roles.add(role)
        in_crop = [position for position in positions if span[0] <= position[0] and position[1] <= span[1]]
        reliable = reference["unique_in_named_panel"] and len(positions) == 1 and bool(in_crop)
        if reliable:
            crop_roles.add(role)
        row["observations"].append({
            "window": window, "role": role, "original_positions": positions,
            "inside_crop": in_crop, "qualifies": reliable,
        })
    if len(observed_roles) > 1:
        row["reason"] = "mixed_source_contexts"
    elif crop_roles:
        row["state"] = f"{next(iter(crop_roles))}_supported"
        row["reason"] = "unique_reference_window_in_original_and_crop"
    elif observed_roles:
        row["reason"] = "source_window_ambiguous_or_outside_crop"
    else:
        row["reason"] = "no_source_window_observed"
    return row


def audit(path, retention_policy=EXCLUDE_COMPETING, placement_binary=None, *, require_reference_roles=False):
    data = path.read_bytes()
    manifest = json.loads(data)
    scope = manifest.get("scope")
    require(manifest.get("schema") == MANIFEST and isinstance(scope, str) and scope.strip(),
            "unknown manifest schema or empty named scope")
    target = manifest.get("target")
    require(target is None or target in ("mitochondrion", "plastid"), "invalid source audit target")
    require(retention_policy in (EXCLUDE_COMPETING, TARGET_ONLY), "unknown retention policy")
    strict = retention_policy == TARGET_ONLY
    require(not strict or target is not None, "target-supported-only requires a named organelle target")
    require(not require_reference_roles or strict,
            "required reference roles require the target-supported-only retention policy")
    reference_scope = declared_reference_scope(manifest, target) if require_reference_roles else None
    missing_roles = reference_scope is not None and bool(reference_scope["missing_required_roles"])
    base = path.parent
    windows = {} if missing_roles else load_windows(manifest, load_panel(manifest, base), base, placement_binary)
    seeds = load_seeds(manifest["seeds"], base)
    strict_names = None
    if not missing_roles:
        if strict:
            originals, strict_names = load_strict_originals(manifest["raw_reads"], base, seeds)
        else:
            originals = load_originals(manifest["raw_reads"], base, {seed[0] for seed in seeds})
        locate = placements if placement_binary is None else original_placement_lookup(placement_binary, originals, windows)
    rows = []
    kept = []
    counts = {"target_supported": 0, "competing_supported": 0, "unresolved": 0}
    retained_sha = sha256()
    for name, span, sequence, fasta in seeds:
        if missing_roles:
            row = unresolved_seed(name, span, sequence, "missing_required_reference_roles")
        else:
            row = assess_seed(name, span, sequence, originals, windows, strict_names, locate)
        rows.append(row)
        counts[row["state"]] += 1
        allowed = row["state"] == "target_supported" if strict else row["state"] != "competing_supported"
        if allowed:
            kept.append(fasta)
            retained_sha.update(fasta)
    report = {
        "schema": MANIFEST, "scope": scope, "manifest_sha256": sha256(data).hexdigest(),
        "source_limit": "support is restricted to declared roles and unique exact windows in the named panel",
        "coordinate_system": "0-based half-open original record and reference coordinates",
        "placement_limit": "at most two exact placements per sequence or panel",
        "inputs": {key: manifest[key] for key in ("seeds", "raw_reads", "references")},
        "windows": {name: value[2] for name, value in windows.items()},
        "counts": counts, "retained_records": len(kept), "records": rows,
        "retained_fasta_sha256": retained_sha.hexdigest(),
    }
    if target is not None:
        report["target"] = target
    if strict:
        report["retention_policy"] = "target_supported_only"
    if reference_scope is not None:
        report["reference_panel_scope"] = reference_scope
        report["source_limit"] = (
            "support is restricted to supplied named reference records and unique exact windows. "
            "Reference roles are user declared. Role presence does not establish panel "
            "completeness or physical origin."
        )
    return report, kept


def publish(outputs):
    staged = []
    linked = []
    try:
        for path, parts in outputs:
            with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
                temp = Path(stream.name)
                staged.append((path, temp))
                for part in parts:
                    stream.write(part)
        for path, temp in staged:
            os.link(temp, path)
            linked.append(path)
    except BaseException:
        for path in linked:
            path.unlink()
        raise
    finally:
        for _, temp in staged:
            temp.unlink()


def resolve_placement_binary(explicit):
    if explicit is None:
        found = shutil.which("ovasm")
        require(found is not None, "OVASM executable not found on PATH. Install OVASM or pass "
                "--exact-placement-bin PATH. Use --python-placements for the Python comparison backend.")
        explicit = Path(found)
    try:
        binary = explicit.resolve(strict=True)
    except (OSError, RuntimeError) as error:
        raise ValueError(f"cannot resolve OVASM executable {explicit}: {error}") from error
    require(binary.is_file() and os.access(binary, os.X_OK), f"OVASM executable is not an executable file: {binary}")
    return binary


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--out", type=Path, required=True, help="new JSON audit path")
    parser.add_argument("--out-fasta", type=Path, help="new FASTA containing crops admitted by the retention policy")
    parser.add_argument("--retention-policy", choices=(EXCLUDE_COMPETING, TARGET_ONLY), default=EXCLUDE_COMPETING,
                        help="retain target-supported and unresolved crops, or require positive target evidence")
    parser.add_argument("--require-reference-roles", action="store_true",
                        help="require explicitly declared target and nuclear reference roles before retaining crops")
    backend = parser.add_mutually_exclusive_group()
    backend.add_argument("--exact-placement-bin", type=Path,
                         help="OVASM executable for exact batch placement (default: ovasm on PATH)")
    backend.add_argument("--python-placements", action="store_true",
                         help="explicitly use the Python comparison backend")
    args = parser.parse_args()
    try:
        paths = [args.out] + ([args.out_fasta] if args.out_fasta else [])
        require(len({path.resolve() for path in paths}) == len(paths), "output paths must differ")
        for path in paths:
            require(not (path.exists() or path.is_symlink()), f"output already exists: {path}")
        binary = None if args.python_placements else resolve_placement_binary(args.exact_placement_bin)
        report, retained = audit(args.manifest, args.retention_policy, binary,
                                 require_reference_roles=args.require_reference_roles)
        outputs = [(args.out, [json.dumps(report, indent=2, sort_keys=True).encode() + b"\n"])]
        if args.out_fasta:
            outputs.append((args.out_fasta, retained))
        publish(outputs)
    except (OSError, ValueError, TypeError, KeyError, AttributeError) as error:
        parser.exit(1, f"{parser.prog}: {error}\n")


if __name__ == "__main__":
    main()
