#!/usr/bin/env python3
"""Prepare exact original-read exclusions from registered exterior source windows."""

from argparse import ArgumentParser
from collections import Counter
import json
from pathlib import Path

from discover_graph_frontier import fastx
from prepare_read_admission import (
    DNA, ROLES, LockedFile, digest, fasta, open_input, read_json, require,
    reverse_complement, sequence_occurrences,
)
from validate_source_exclusion import new_output, publish_outputs


PLAN_SCHEMA = "ovasm.read-source-exclusion-plan.v1"
REPORT_SCHEMA = "ovasm.read-source-exclusion-report.v1"
MANIFEST_SCHEMA = "ovasm.read-source-exclusions.v1"


def read_panel(entries, lock, target):
    require(isinstance(entries, list), "references must be an array")
    panel = {}
    identities = set()
    for entry in entries:
        require(isinstance(entry, dict) and set(entry) == {"id", "role", "path", "sha256"},
                "invalid reference fields")
        name, role = entry["id"], entry["role"]
        require(isinstance(name, str) and name and name not in panel and role in ROLES,
                "invalid reference identifier or role")
        source = lock({key: entry[key] for key in ("path", "sha256")})
        require(source.state[:2] not in identities, "repeated reference input file")
        identities.add(source.state[:2])
        panel[name] = role, dict(fasta(source)), source.binding()
    roles = {role for role, _, _ in panel.values()}
    require(target in roles and len(roles) >= 2, "provide target and competing reference roles")
    return panel


def reference_placements(panel, word):
    result = []
    for reference, (role, records, _) in panel.items():
        for name, sequence in records.items():
            result.extend(dict(hit, panel_reference=reference, role=role)
                          for hit in sequence_occurrences(sequence, name, word, 2 - len(result)))
            if len(result) == 2:
                return result
    return result


def read_windows(entries, panel, target, lock):
    require(isinstance(entries, list) and entries, "provide registered source windows")
    windows = {}
    fields = {"id", "reference", "record", "interval", "boundary", "provenance"}
    for entry in entries:
        require(isinstance(entry, dict) and isinstance(entry.get("reference"), str)
                and entry["reference"] in panel, "invalid source window reference")
        reference = entry["reference"]
        role, records, _ = panel[reference]
        competing = role != target
        require(set(entry) == fields | ({"exterior_side"} if competing else set()),
                "competing windows require an exterior side; target windows are veto contexts")
        name, record, interval, boundary = (entry[k] for k in ("id", "record", "interval", "boundary"))
        require(isinstance(name, str) and name and name not in windows and record in records,
                "invalid source window identifier or record")
        require(isinstance(interval, list) and len(interval) == 2
                and all(type(i) is int for i in [*interval, boundary]), "invalid source window coordinates")
        lo, hi = interval
        require(0 <= lo < boundary < hi <= len(records[record]), "source boundary must lie inside its window")
        sequence = records[record]
        word = bytes(sequence[lo:hi])
        require(DNA.fullmatch(word), "source window must be unambiguous ACGT")
        evidence = lock(entry["provenance"])
        placements = reference_placements(panel, word)
        intended = {"panel_reference": reference, "reference": record, "start0": lo,
                    "end0": hi, "strand": "+", "role": role}
        exterior = None
        exterior_digest = None
        exterior_other_roles = []
        if competing:
            side = entry["exterior_side"]
            require(side in ("left", "right"), "invalid competing exterior side")
            exterior = [lo, boundary] if side == "left" else [boundary, hi]
            flank = bytes(sequence[slice(*exterior)])
            exterior_digest = digest(flank)
            rev = reverse_complement(flank)
            exterior_other_roles = sorted({
                r for r, sequences, _ in panel.values()
                if r != role and any(flank in s or rev in s for s in sequences.values())
            })
        qualified = competing and word != reverse_complement(word) and placements == [intended] and not exterior_other_roles
        windows[name] = word, {
            "reference": reference, "record": record, "role": role, "interval": interval,
            "boundary": boundary, "word_sha256": digest(word), "provenance": evidence.binding(),
            "reference_placements": placements, "exterior_interval": exterior,
            "exterior_sha256": exterior_digest, "exterior_other_roles": exterior_other_roles,
            "qualifies_as_competing": qualified,
        }
    require(any(info["role"] != target for _, info in windows.values()), "provide a competing source window")
    return windows


def original_exterior(hit, window):
    if window["exterior_interval"] is None:
        return None
    origin = window["interval"][0]
    lo, hi = window["exterior_interval"]
    if hit["strand"] == "+":
        return [hit["start0"] + lo - origin, hit["start0"] + hi - origin]
    return [hit["end0"] - (hi - origin), hit["end0"] - (lo - origin)]


def observe(sequence, windows, target):
    observed = {}
    for name, (word, window) in windows.items():
        hits = list(sequence_occurrences(sequence, name, word, 2))
        if hits:
            observed[name] = {
                "role": window["role"],
                "original_positions": hits,
                "original_exterior_intervals": [original_exterior(hit, window) for hit in hits],
            }
    roles = {o["role"] for o in observed.values()}
    if not observed:
        return "no_registered_window", observed
    if len(roles) > 1:
        return "mixed_registered_roles", observed
    if roles == {target}:
        return "target_context_only", observed
    if any(len(o["original_positions"]) != 1 or not windows[n][1]["qualifies_as_competing"]
           for n, o in observed.items()):
        return "ambiguous_or_shared_competing_context", observed
    return "competing_context", observed


def scan_originals(inputs, windows, target):
    counts = Counter()
    records = []
    candidates = Counter()
    for file_index, locked in enumerate(inputs):
        with open_input(locked) as source:
            for ordinal, (header, sequence, quality) in enumerate(fastx(source, preserve_original=True)):
                reason, observed = observe(sequence.upper(), windows, target)
                counts[reason] += 1
                if not observed:
                    continue
                if reason == "competing_context":
                    candidates[header] = 0
                records.append({
                    "read_id_hex": header.hex(), "sequence_sha256": digest(sequence),
                    "quality_sha256": None if quality is None else digest(quality),
                    "original_input_index0": file_index, "original_record_ordinal0": ordinal,
                    "original_length": len(sequence), "state": reason, "observations": observed,
                })
    # An earlier record with the same full header may have no registered
    # window. Count selected headers in the full input before exporting.
    if candidates:
        for locked in inputs:
            with open_input(locked) as source:
                for header, _, _ in fastx(source, preserve_original=True):
                    if header in candidates:
                        candidates[header] += 1
    for record in records:
        if record["state"] == "competing_context" and candidates[bytes.fromhex(record["read_id_hex"])] != 1:
            record["state"] = "duplicate_original_header"
            counts["competing_context"] -= 1
            counts["duplicate_original_header"] += 1
    return counts, records


def audit(plan_path):
    locked_files = []

    def lock(spec):
        source = LockedFile.lock(plan_path.resolve().parent, spec)
        locked_files.append(source)
        return source

    plan_file = lock({"path": str(plan_path.resolve()), "sha256": digest(plan_path.read_bytes())})
    plan = read_json(plan_file)
    fields = {"schema", "scope", "target", "references", "raw_reads", "windows"}
    require(isinstance(plan, dict) and set(plan) == fields and plan["schema"] == PLAN_SCHEMA,
            "invalid source exclusion plan schema or fields")
    scope, target = plan["scope"], plan["target"]
    require(isinstance(scope, str) and scope.strip() and target in ROLES[:2], "invalid scope or organelle target")
    require(isinstance(plan["raw_reads"], list) and plan["raw_reads"], "provide original FASTX inputs")
    inputs = [lock(spec) for spec in plan["raw_reads"]]
    require(len({item.state[:2] for item in inputs}) == len(inputs), "repeated original input file")
    panel = read_panel(plan["references"], lock, target)
    windows = read_windows(plan["windows"], panel, target, lock)
    counts, records = scan_originals(inputs, windows, target)
    for source in locked_files:
        source.check()
    return {
        "schema": REPORT_SCHEMA, "plan": plan_file.binding(), "scope": scope, "target": target,
        "source_limit": "Exact registered source contexts in the supplied linear reference records and originals.",
        "coordinates": "zero-based half-open sequence intervals and zero-based original ordinals",
        "raw_reads": [item.binding() for item in inputs],
        "references": [dict(binding, id=name, role=role) for name, (role, _, binding) in panel.items()],
        "windows": {name: detail for name, (_, detail) in windows.items()},
        "counts": counts, "records": records, "eligible_records": counts["competing_context"],
    }


def json_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2) + "\n").encode()


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path)
    parser.add_argument("--out", type=Path, required=True, help="new source-evidence report")
    parser.add_argument("--out-exclusions", type=Path, required=True, help="new manifest if any reads qualify")
    args = parser.parse_args()
    try:
        out, manifest = new_output(args.out), new_output(args.out_exclusions)
        require(out != manifest, "output paths must differ")
        report = audit(args.plan)
        report["manifest_created"] = bool(report["eligible_records"])
        data = json_bytes(report)
        outputs = [(out, data)]
        if report["manifest_created"]:
            proof = {"path": str(out), "sha256": digest(data)}
            fields = ("read_id_hex", "sequence_sha256", "quality_sha256")
            records = [{key: row[key] for key in fields} | {"target": report["target"], "provenance": proof}
                       for row in report["records"] if row["state"] == "competing_context"]
            outputs.append((manifest, json_bytes({
                "schema": MANIFEST_SCHEMA, "scope": report["scope"], "records": records,
            })))
        publish_outputs(outputs)
    except (OSError, ValueError, TypeError, KeyError) as error:
        parser.exit(1, f"source exclusion audit: {error}\n")
    print(f'{sum(report["counts"].values())} originals, {report["eligible_records"]} exclusions')


if __name__ == "__main__":
    main()
