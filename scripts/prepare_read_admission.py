#!/usr/bin/env python3
"""Screen exact original read ends at registered plant organelle junctions."""
import argparse
from collections import Counter
from contextlib import contextmanager
from dataclasses import dataclass
import gzip
import hashlib
import json
import os
from pathlib import Path
import re

ROLES = ("mitochondrion", "plastid", "nuclear")
SCREEN_SCHEMA = "ovasm.read-source-screen.v1"
QUALIFICATION_SCHEMA = "ovasm.read-source-qualification.v1"
SHA256 = re.compile(r"[0-9a-f]{64}")
READ_ID = re.compile(r"(?:[0-9a-f]{2})+")
DNA = re.compile(rb"[ACGT]+")
COMPLEMENT = bytes.maketrans(b"ACGT", b"TGCA")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def file_state(stat):
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


@dataclass(frozen=True)
class LockedFile:
    path: Path
    sha256: str
    state: tuple

    @classmethod
    def lock(cls, parent, binding):
        require(isinstance(binding, dict) and set(binding) == {"path", "sha256"},
                "file binding requires path and sha256")
        require(isinstance(binding["path"], str) and binding["path"].strip(), "empty file path")
        expected = binding["sha256"]
        require(isinstance(expected, str) and SHA256.fullmatch(expected),
                "expected lowercase SHA-256")
        path = (parent / binding["path"]).resolve(strict=True)
        require(path.is_file(), f"expected a regular input file: {path}")
        with path.open("rb") as stream:
            before = file_state(os.fstat(stream.fileno()))
            actual = hashlib.file_digest(stream, "sha256").hexdigest()
        require(actual == expected, f"input SHA-256 mismatch: {path}")
        locked = cls(path, actual, before)
        locked.check()
        return locked

    def check(self):
        require(file_state(self.path.stat()) == self.state, f"input changed during screening: {self.path}")

    def binding(self):
        return {"path": str(self.path), "sha256": self.sha256}


@contextmanager
def open_input(locked):
    locked.check()
    with locked.path.open("rb") as stream:
        require(file_state(os.fstat(stream.fileno())) == locked.state,
                f"input changed during screening: {locked.path}")
        compressed = stream.read(2) == b"\x1f\x8b"
        stream.seek(0)
        if compressed:
            with gzip.GzipFile(fileobj=stream) as raw:
                yield raw
        else:
            yield stream
    locked.check()


def unique_keys(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, f"duplicate JSON field: {key}")
        result[key] = value
    return result


def read_json(locked):
    with open_input(locked) as stream:
        return json.load(stream, object_pairs_hook=unique_keys)


def fasta(locked):
    with open_input(locked) as stream:
        name = None
        seen = set()
        seq = bytearray()
        for line in stream:
            if line.startswith(b">"):
                if name is not None:
                    require(seq, f"empty FASTA sequence: {name}")
                    yield name, seq
                parts = line[1:].split(None, 1)
                require(parts, f"empty FASTA identifier: {locked.path}")
                name = parts[0].decode()
                require(name not in seen, f"duplicate FASTA identifier: {name}")
                seen.add(name)
                seq = bytearray()
            elif line.strip():
                require(name is not None, f"sequence before FASTA identifier: {locked.path}")
                seq.extend(b"".join(line.split()).upper())
        require(name is not None and seq, f"missing FASTA sequence: {locked.path}")
        yield name, seq


def fastq(locked):
    with open_input(locked) as stream:
        while head := stream.readline():
            seq = stream.readline().rstrip(b"\r\n")
            plus = stream.readline()
            qual = stream.readline().rstrip(b"\r\n")
            rid = head.rstrip(b"\r\n")[1:]
            require(head.startswith(b"@") and rid and plus.startswith(b"+")
                    and len(seq) == len(qual) > 0, f"malformed FASTQ record in {locked.path}")
            yield rid, seq, qual


@dataclass(frozen=True)
class Window:
    id: str
    reference: str
    start0: int
    boundary0: int
    end0: int

    @classmethod
    def parse(cls, value):
        require(isinstance(value, dict) and set(value) == {
            "id", "reference", "start0", "boundary0", "end0"}, "invalid junction window fields")
        require(all(isinstance(value[k], str) and value[k].strip() for k in ("id", "reference")),
                "junction window requires id and reference")
        require(all(type(value[k]) is int for k in ("start0", "boundary0", "end0"))
                and 0 <= value["start0"] < value["boundary0"] < value["end0"],
                "junction window must contain bases on both sides of boundary0")
        return cls(**value)

    def binding(self):
        return vars(self)


def reverse_complement(word):
    return word.translate(COMPLEMENT)[::-1]


def validated_scope(obj, schema):
    require(isinstance(obj, dict) and obj.get("schema") == schema, f"expected {schema}")
    require(isinstance(obj.get("scope"), str) and obj["scope"].strip(), "scope is required")


@dataclass
class Qualification:
    evidence: LockedFile
    provenance: LockedFile
    scope: str
    records: dict

    @classmethod
    def load(cls, parent, spec, target, panel, windows):
        evidence = LockedFile.lock(parent, spec)
        raw = read_json(evidence)
        validated_scope(raw, QUALIFICATION_SCHEMA)
        require(set(raw) == {"schema", "scope", "target", "panel_sha256", "provenance", "records"},
                "invalid qualification fields")
        require(raw["target"] == target and raw["panel_sha256"] == {
            role: lock.sha256 for role, lock in panel.items()}, "qualification panel or target differs")
        provenance = LockedFile.lock(evidence.path.parent, raw["provenance"])
        require(isinstance(raw["records"], list), "qualification records must be an array")
        records = {}
        for row in raw["records"]:
            require(isinstance(row, dict) and set(row) == {
                "read_id_hex", "sequence_sha256", "quality_sha256", "window", "verdict"},
                "invalid qualification record fields")
            require(isinstance(row["read_id_hex"], str) and READ_ID.fullmatch(row["read_id_hex"]),
                    "qualification requires a lowercase hexadecimal full FASTQ header")
            require(all(isinstance(row[key], str) and SHA256.fullmatch(row[key])
                        for key in ("sequence_sha256", "quality_sha256")),
                    "qualification requires sequence and quality SHA-256")
            window = Window.parse(row["window"])
            require(windows.get(window.id) == window, "qualification junction window differs")
            require(row["verdict"] in ("target_supported", "competing", "unresolved"),
                    "invalid qualification verdict")
            key = (bytes.fromhex(row["read_id_hex"]), window.id)
            require(key not in records, "duplicate qualification record and junction")
            records[key] = row
        return cls(evidence, provenance, raw["scope"], records)

    def binding(self):
        return {"declarations": self.evidence.binding(), "provenance": self.provenance.binding(),
                "scope": self.scope}


@dataclass
class ScreenPlan:
    source: LockedFile
    scope: str
    target: str
    reads: list
    panel: dict
    windows: dict
    qualification: Qualification | None

    @classmethod
    def load(cls, path):
        path = path.resolve(strict=True)
        before = file_state(path.stat())
        raw = path.read_bytes()
        source = LockedFile(path, digest(raw), before)
        source.check()
        obj = json.loads(raw, object_pairs_hook=unique_keys)
        validated_scope(obj, SCREEN_SCHEMA)
        expected = {"schema", "scope", "target", "reads", "panel", "windows"}
        require(expected <= set(obj) <= expected | {"qualification"}, "invalid screening plan fields")
        target = obj["target"]
        require(target in ROLES[:2], "target must be mitochondrion or plastid")
        require(isinstance(obj["reads"], list) and obj["reads"], "at least one original FASTQ is required")
        reads = [LockedFile.lock(path.parent, spec) for spec in obj["reads"]]
        require(isinstance(obj["panel"], dict) and set(obj["panel"]) == set(ROLES),
                "screening requires mitochondrial, plastid and nuclear panel references")
        panel = {role: LockedFile.lock(path.parent, obj["panel"][role]) for role in ROLES}
        require(isinstance(obj["windows"], list) and obj["windows"], "at least one junction window is required")
        windows = {}
        for value in obj["windows"]:
            window = Window.parse(value)
            require(window.id not in windows, f"duplicate junction window ID: {window.id}")
            windows[window.id] = window
        qualification = (Qualification.load(path.parent, obj["qualification"], target, panel, windows)
                         if "qualification" in obj else None)
        return cls(source, obj["scope"], target, reads, panel, windows, qualification)

    def check(self):
        for locked in (self.source, *self.reads, *self.panel.values()):
            locked.check()
        if self.qualification:
            self.qualification.evidence.check()
            self.qualification.provenance.check()

    def binding(self):
        return {"plan": self.source.binding(), "scope": self.scope, "target": self.target,
                "reads": [read.binding() for read in self.reads],
                "panel": {role: ref.binding() for role, ref in self.panel.items()},
                "windows": [window.binding() for window in self.windows.values()],
                "qualification": self.qualification.binding() if self.qualification else None}


def window_words(plan):
    words = {}
    by_ref = {}
    for window in plan.windows.values():
        by_ref.setdefault(window.reference, []).append(window)
    for name, seq in fasta(plan.panel[plan.target]):
        for window in by_ref.get(name, ()):
            word = bytes(seq[window.start0:window.end0])
            require(len(word) == window.end0 - window.start0 and DNA.fullmatch(word),
                    f"junction window is out of bounds or ambiguous: {window.id}")
            words[window.id] = word
    require(set(words) == set(plan.windows), "junction window reference missing from target panel")
    return words


@dataclass
class Candidate:
    id: bytes
    word: bytes
    sequence_sha256: str
    quality_sha256: str
    windows: list

    def binding(self):
        return {"read_id_hex": self.id.hex(), "sequence_sha256": self.sequence_sha256,
                "quality_sha256": self.quality_sha256}


def find_candidates(reads, words, cap):
    patterns = [(key, word, reverse_complement(word)) for key, word in words.items()]
    candidates = {}
    scan = []
    for source in reads:
        count = 0
        for rid, seq, qual in fastq(source):
            count += 1
            normalized = seq.upper()
            hits = [key for key, word, rc in patterns if word in normalized or rc in normalized]
            if not hits:
                continue
            require(rid not in candidates, f"duplicate candidate FASTQ header: {rid.hex()}")
            if len(candidates) == cap:
                scan.append({"file": source.binding(), "records": count})
                return None, scan
            candidates[rid] = Candidate(rid, normalized, digest(seq), digest(qual), hits)
        scan.append({"file": source.binding(), "records": count})
    # A duplicate header outside a window still makes the original identity ambiguous.
    if candidates:
        seen = set()
        for source in reads:
            for rid, _, _ in fastq(source):
                if rid in candidates:
                    require(rid not in seen, f"duplicate candidate FASTQ header: {rid.hex()}")
                    seen.add(rid)
        require(len(seen) == len(candidates), "candidate FASTQ header disappeared during screening")
    return candidates, scan


def sequence_occurrences(seq, name, word, remaining):
    reverse = reverse_complement(word)
    strands = (("+", word),) if reverse == word else (("+", word), ("-", reverse))
    for strand, pattern in strands:
        pos = seq.find(pattern)
        while pos != -1:
            yield {"reference": name, "start0": pos, "end0": pos + len(pattern), "strand": strand}
            remaining -= 1
            if not remaining:
                return
            pos = seq.find(pattern, pos + 1)


def search_panel(panel, words):
    """Keep at most two positions per role. Two suffices to refute uniqueness."""
    placements = {word: {role: [] for role in ROLES} for word in words}
    for role, ref in panel.items():
        for name, seq in fasta(ref):
            for word, by_role in placements.items():
                matches = by_role[role]
                remaining = 2 - len(matches)
                if remaining:
                    matches.extend(sequence_occurrences(seq, name, word, remaining))
    return placements


def unique_target(placements, target):
    if len(placements[target]) != 1 or any(placements[role] for role in ROLES if role != target):
        return None
    return placements[target][0]


def qualification_verdict(qualification, candidate, window_id):
    if not qualification:
        return "missing"
    row = qualification.records.get((candidate.id, window_id))
    if not row:
        return "missing"
    if (row["sequence_sha256"], row["quality_sha256"]) != (
            candidate.sequence_sha256, candidate.quality_sha256):
        return "record_mismatch"
    return row["verdict"]


def decide_candidate(plan, candidate, words, placements):
    end_hits = placements.get(candidate.word)
    end_target = unique_target(end_hits, plan.target) if end_hits else None
    windows = {}
    for key in candidate.windows:
        window = plan.windows[key]
        hits = placements[words[key]]
        window_target = unique_target(hits, plan.target)
        direct = bool(end_target and window_target
                      and window_target["reference"] == end_target["reference"] == window.reference
                      and window_target["start0"] == window.start0
                      and window_target["end0"] == window.end0
                      and end_target["start0"] <= window.start0 < window.end0 <= end_target["end0"])
        verdict = qualification_verdict(plan.qualification, candidate, key)
        windows[key] = {"exact_target_span": direct, "qualification": verdict}
    if any(w["qualification"] == "competing" for w in windows.values()):
        status = "competing_qualification"
    elif any(w["qualification"] == "unresolved" for w in windows.values()):
        status = "unresolved_qualification"
    elif end_hits and any(end_hits[r] for r in ROLES if r != plan.target):
        status = "competing_exact_record"
    elif any(w["exact_target_span"] and w["qualification"] == "target_supported" for w in windows.values()):
        status = "target_supported"
    elif any(w["exact_target_span"] for w in windows.values()):
        status = "unresolved_qualification"
    else:
        status = "unresolved_exact_placement"
    return candidate.binding() | {"full_record_placements": end_hits, "windows": windows, "status": status}


def screen(plan, cap):
    words = window_words(plan)
    candidates, scan = find_candidates(plan.reads, words, cap)
    report = plan.binding() | {"schema": "ovasm.read-source-screen-report.v1",
                              "candidate_cap": cap, "original_records_scanned": scan}
    if candidates is None:
        return report | {"state": "truncated", "reason": "candidate_cap"}, []
    patterns = set(words.values()) | {c.word for c in candidates.values() if DNA.fullmatch(c.word)}
    placements = search_panel(plan.panel, sorted(patterns))
    decisions = [decide_candidate(plan, candidate, words, placements)
                 for _, candidate in sorted(candidates.items())]
    selected = [{key: row[key] for key in ("read_id_hex", "sequence_sha256", "quality_sha256")}
                for row in decisions if row["status"] == "target_supported"]
    report |= {"state": "complete", "placements_kept_per_role": 2,
               "competition_scope": "literal exact positions in the named panel; records with ambiguous bases do not give a complete exact witness",
               "windows_exact_placements": {key: placements[word] for key, word in words.items()},
               "window_sequence_sha256": {key: digest(word) for key, word in words.items()},
               "decision_counts": dict(Counter(row["status"] for row in decisions)),
               "records": decisions}
    return report, selected


def write_json(path, obj):
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(obj, stream, indent=2, sort_keys=True)
        stream.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path, help="hash-bound source screen plan JSON")
    parser.add_argument("--out", type=Path, required=True, help="new directory for audit and optional native manifest")
    parser.add_argument("--max-candidates", type=int, default=1024,
                        help="junction-matching records to compare before declaring truncation")
    args = parser.parse_args()
    try:
        require(args.max_candidates > 0, "--max-candidates must be positive")
        require(not args.out.exists(), f"output directory already exists: {args.out}")
        plan = ScreenPlan.load(args.plan)
        report, selected = screen(plan, args.max_candidates)
        plan.check()
        out = args.out.resolve()
        out.mkdir(parents=True, exist_ok=False)
        path = out / "source-screen.json"
        write_json(path, report)
        if selected:
            proof = {"path": str(path), "sha256": digest(path.read_bytes())}
            write_json(out / "read-admission.json", {
                "schema": "ovasm.read-admission.v1", "scope": plan.scope,
                "records": [row | {"target": plan.target, "provenance": proof} for row in selected]})
    except (OSError, ValueError) as error:
        parser.exit(1, f"read source screen: {error}\n")


if __name__ == "__main__":
    main()
