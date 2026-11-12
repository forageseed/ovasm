"""Verify explicit original frontier events before source-model eligibility."""

from collections import defaultdict
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import re

from discover_graph_frontier import (
    MAX_SCAN_K, MISSING, REPORT, events, graph_words, json_object, parse_graph, revcomp,
)
from prepare_read_admission import LockedFile, open_input, read_json
from source_competition import contains


SCHEMA = "ovasm.frontier-event-bindings.v1"
DIGEST = re.compile(r"[0-9a-f]{64}")
HEADER = re.compile(r"(?:[0-9a-f]{2})+")
DNA = re.compile(rb"[ACGT]+")
IDENTITY_FIELDS = ("read_id_hex", "sequence_sha256", "quality_sha256", "molecule")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(value):
    return isinstance(value, str) and bool(DIGEST.fullmatch(value))


def stable_hash(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def match_file(locked, binding):
    return (isinstance(binding, dict) and isinstance(binding.get("path"), str)
            and Path(binding["path"]).is_absolute() and Path(binding["path"]).resolve() == locked.path
            and binding.get("sha256") == locked.sha256)


def nonnegative(value):
    return type(value) is int and value >= 0


@dataclass(frozen=True)
class Registration:
    input_index0: int
    record_index0: int
    read_id_hex: str
    sequence_sha256: str
    quality_sha256: str
    molecule: str
    event_index0: int
    event_sha256: str
    window_id: str
    side: str
    paf_row_sha256: str

    @classmethod
    def parse(cls, obj, windows, k, original_count):
        require(isinstance(obj, dict) and set(obj) == set(cls.__dataclass_fields__),
                "invalid frontier registration fields")
        require(all(nonnegative(obj[name]) for name in ("input_index0", "record_index0", "event_index0"))
                and obj["input_index0"] < original_count, "invalid frontier original or event ordinal")
        require(isinstance(obj["read_id_hex"], str) and HEADER.fullmatch(obj["read_id_hex"]),
                "frontier registration requires the hexadecimal full original FASTQ header")
        require(all(digest(obj[name]) for name in ("sequence_sha256", "quality_sha256",
                    "event_sha256", "paf_row_sha256")), "invalid frontier registration SHA-256")
        require(isinstance(obj["molecule"], str) and obj["molecule"].strip(),
                "frontier registration requires an original molecule identity")
        require(obj["side"] in ("left", "right"), "invalid frontier registration side")
        key = obj["window_id"]
        require(isinstance(key, str) and key in windows, "unknown frontier registration window")
        window = windows[key]
        require(window.mode == "exact_context" and window.require_full_query and window.graph_required
                and set(window.exclusive_sides) == {"left", "right"}
                and window.anchor_bp == k and window.interval[1] - window.interval[0] == k,
                "frontier window must be a full-query exact graph core at the frozen discovery k")
        return cls(**obj)

    @property
    def ordinal(self):
        return self.input_index0, self.record_index0

    @property
    def query(self):
        tokens = bytes.fromhex(self.read_id_hex).split()
        require(bool(tokens), "frontier registration header has no mapper query name")
        return tokens[0].decode("ascii")

    def check_identity(self, record, fields=IDENTITY_FIELDS):
        require(all(record.get(key) == getattr(self, key) for key in fields),
                "frontier registration does not bind the original identity")

    def hash_event(self, event):
        return stable_hash({key: getattr(self, key) for key in IDENTITY_FIELDS}
                           | {"event_index0": self.event_index0, "event": event})


@dataclass(frozen=True)
class Context:
    raw_core: tuple[int, int]
    raw_affected: tuple[int, int]
    raw_footprint: tuple[int, int]
    graph_core: tuple[str, int, int]
    graph_footprint: tuple[str, int, int]
    word: bytes


def proximal_context(sequence, event, side, k, graph_sequences, words):
    affected = tuple(event["affected_raw_interval0"])
    start, stop = event["absent_word_start_interval0"]
    require(affected == (start, stop + k - 1) and contains((0, len(sequence)), affected),
            "frontier affected interval does not bind its absent word starts")
    anchor = event[f"{side}_anchor"]
    require(anchor is not None, "frontier registration has no unique discovery anchor")
    raw_anchor = tuple(anchor["raw_interval0"])
    raw_start = start - 1 if side == "left" else stop
    node = anchor["segment"] + anchor["orientation"]
    graph_start, graph_end = anchor["oriented_interval0"]
    require(raw_anchor == (raw_start, raw_start + k) and contains((0, len(sequence)), raw_anchor)
            and graph_end - graph_start == k
            and words.get(sequence[slice(*raw_anchor)], MISSING) == (node, graph_start),
            "frontier discovery anchor is not unique at its bound graph coordinate")
    a, b = affected
    core = (a - k, a) if side == "left" else (b, b + k)
    require(contains((0, len(sequence)), core), "frontier raw record has no disjoint proximal core")
    footprint = min(raw_anchor[0], core[0]), max(raw_anchor[1], core[1])
    footprint_start = graph_start + footprint[0] - raw_anchor[0]
    footprint_end = footprint_start + footprint[1] - footprint[0]
    word = sequence[slice(*footprint)]
    graph_word = graph_sequences[node]
    require(contains((0, len(graph_word)), (footprint_start, footprint_end)) and DNA.fullmatch(word)
            and word == graph_word[footprint_start:footprint_end],
            "frontier proximal footprint is not a literal contiguous single-node spelling")
    core_start = graph_start + core[0] - raw_anchor[0]
    require(words.get(sequence[slice(*core)], MISSING) == (node, core_start),
            "frontier proximal core does not have a unique graph placement")
    return Context(core, affected, footprint, (node, core_start, core_start + k),
                   (node, footprint_start, footprint_end), min(word, revcomp(word)))


class FrontierBindings:
    """Rebind registrations to a complete discovery, original bytes and one PAF row."""

    @classmethod
    def load(cls, locked_binding: LockedFile, lock, graph, originals, windows):
        raw = read_json(locked_binding)
        require(isinstance(raw, dict) and set(raw) == {"schema", "discovery_report", "discovery_evidence", "records"}
                and raw["schema"] == SCHEMA, "invalid frontier event bindings schema or fields")
        report_file = lock(raw["discovery_report"])
        evidence_file = lock(raw["discovery_evidence"])
        report = read_json(report_file)
        require(isinstance(report, dict) and report.get("schema") == REPORT
                and report.get("raw_scan_status") == "complete", "frontier discovery raw scan is incomplete")
        require(match_file(graph, report.get("graph")), "frontier discovery uses a different graph")
        require(match_file(evidence_file, report.get("evidence")), "frontier discovery evidence binding differs")
        source_reads = report.get("raw_reads")
        require(isinstance(source_reads, list) and len(source_reads) == len(originals)
                and all(match_file(original, registered) for original, registered in zip(originals, source_reads)),
                "frontier discovery original FASTQ order or binding differs")
        k = report.get("scan_k")
        require(type(k) is int and 1 <= k <= MAX_SCAN_K, "invalid frontier discovery k")
        index = report.get("graph_index")
        require(isinstance(index, dict) and index.get("status") == "complete"
                and nonnegative(index.get("traversal_state_budget")), "frontier discovery graph index is incomplete")
        require(isinstance(raw["records"], list), "frontier registrations must be an array")

        registrations = defaultdict(list)
        keys = set()
        query_ordinals = {}
        molecule_ordinals = {}
        for obj in raw["records"]:
            row = Registration.parse(obj, windows, k, len(originals))
            key = row.ordinal, row.event_index0, row.window_id, row.side, row.paf_row_sha256
            require(key not in keys, "duplicate frontier registration")
            keys.add(key)
            require(query_ordinals.setdefault(row.query, row.ordinal) == row.ordinal
                    and molecule_ordinals.setdefault(row.molecule, row.ordinal) == row.ordinal,
                    "frontier query or molecule binds multiple original ordinals")
            registrations[row.ordinal].append(row)

        ledger = {}
        if registrations:
            with open_input(evidence_file) as stream:
                for line in stream:
                    observed = json_object(line, "frontier discovery evidence row")
                    ordinal = observed.get("input_index0"), observed.get("record_index0")
                    require(all(nonnegative(value) for value in ordinal),
                            "invalid frontier discovery evidence ordinal")
                    if ordinal not in registrations:
                        continue
                    require(ordinal not in ledger, "duplicate original ordinal in frontier discovery evidence")
                    require(observed.get("status") == "graph_anchored_candidate"
                            and isinstance(observed.get("events"), list),
                            "frontier discovery row has no observed raw events")
                    for row in registrations[ordinal]:
                        row.check_identity(observed, IDENTITY_FIELDS[:3])
                        require(row.event_index0 < len(observed["events"])
                                and row.hash_event(observed["events"][row.event_index0]) == row.event_sha256,
                                "frontier registration event does not bind to discovery evidence")
                    ledger[ordinal] = observed
        require(ledger.keys() == registrations.keys(), "frontier registered original missing from discovery evidence")

        with open_input(graph) as stream:
            sequences, adjacency = parse_graph(stream)
        words, states, complete = graph_words(sequences, adjacency, k, index["traversal_state_budget"])
        require(complete and states == index.get("traversal_states")
                and len(words) == index.get("represented_oriented_words")
                and sum(value is not None for value in words.values()) == index.get("unique_oriented_node_words"),
                "frontier complete graph-word index does not reproduce discovery")
        return cls(locked_binding, report_file, evidence_file, k, registrations, ledger, sequences, adjacency, words)

    def __init__(self, binding, report_file, evidence_file, k, registrations, ledger, sequences, adjacency, words):
        self._binding = binding
        self._report_file = report_file
        self._evidence_file = evidence_file
        self._k = k
        self._registrations = registrations
        self._ledger = ledger
        self._graph = sequences, adjacency, words
        self._originals = {}
        self._contexts = {}
        self._donors = defaultdict(set)
        self._graph_originals = set()
        self._footprints = {}
        self._finished = False
        self._verified_projections = 0
        self._donor_supported = 0
        self.queries = {row.query for rows in registrations.values() for row in rows}
        self.paf_hashes = {row.paf_row_sha256 for rows in registrations.values() for row in rows}

    def record_original(self, input_index, ordinal, raw):
        require(not self._finished, "frontier original scan already finished")
        key = input_index, ordinal
        if key not in self._registrations:
            return
        require(key not in self._originals, "registered frontier original visited more than once")
        rows = self._registrations[key]
        observed = self._ledger[key]
        original = raw.binding()
        require(raw.query == rows[0].query and len(raw.seq) == observed.get("read_length"),
                "frontier registration query or original length differs")
        for row in rows:
            row.check_identity(original)
        sequence = raw.seq.upper()
        sequences, adjacency, words = self._graph
        recomputed = events(sequence, words, sequences, adjacency, self._k)[0]
        require(recomputed == observed["events"], "frontier raw events do not reproduce against the frozen graph")
        self._originals[key] = raw
        for row in rows:
            context_key = key, row.event_index0, row.side
            if context_key in self._contexts:
                continue
            context = proximal_context(sequence, recomputed[row.event_index0], row.side, self._k, sequences, words)
            self._contexts[context_key] = context
            self._footprints.setdefault(context.word, revcomp(context.word))

    def finish_originals(self):
        require(not self._finished, "frontier original scan already finished")
        require(self._originals.keys() == self._registrations.keys(), "registered frontier original missing from FASTQs")
        self._graph = None
        self._ledger = None
        self._finished = True

    def record_graph_original(self, raw):
        require(self._finished, "finish the original scan before checking frontier donors")
        require(raw.molecule not in self._graph_originals, "duplicate original molecule in frontier graph inputs")
        self._graph_originals.add(raw.molecule)
        sequence = raw.seq.upper()
        for word, reverse in self._footprints.items():
            if word in sequence or reverse in sequence:
                self._donors[word].add(raw.molecule)

    def verify_alignments(self, paf_rows, panel, windows):
        """Return verified (query, window, alignment) joins with independent footprint donors."""
        require(self._finished, "finish the original scan before verifying frontier alignments")
        require(self.paf_hashes <= paf_rows.keys(), "registered frontier PAF row missing from bound input")
        verified = defaultdict(list)
        self._verified_projections = 0
        self._donor_supported = 0
        for key, rows in self._registrations.items():
            raw = self._originals[key]
            for row in rows:
                aln = paf_rows[row.paf_row_sha256]
                window = windows[row.window_id]
                context = self._contexts[key, row.event_index0, row.side]
                require(aln.query == raw.query and aln.reference == window.reference and aln.reference in panel
                        and aln.query_interval == (0, len(raw.seq)) and all(step.op != "N" for step in aln.steps),
                        "frontier PAF row does not bind one continuous full original and registered reference")
                reference = panel[aln.reference]
                union = min(context.raw_core[0], context.raw_affected[0]), max(context.raw_core[1], context.raw_affected[1])
                require(aln.continuous_over(union), "frontier PAF row does not cover the core and affected interval")
                reference_core = aln.project_exact(context.raw_core, raw.seq, reference)
                reference_affected = aln.project_exact(context.raw_affected, raw.seq, reference)
                require(reference_core == tuple(window.interval) and reference_affected is not None
                        and (reference_core[1] <= reference_affected[0] or reference_affected[1] <= reference_core[0]),
                        "frontier PAF row does not project disjoint exact core and affected intervals")
                self._verified_projections += 1
                donors = sorted(self._donors[context.word] - {raw.molecule})
                if not donors:
                    continue
                self._donor_supported += 1
                verified[raw.query, row.window_id, aln].append({
                    "input_index0": row.input_index0, "record_index0": row.record_index0,
                    "event_index0": row.event_index0, "event_sha256": row.event_sha256,
                    "window_id": row.window_id, "side": row.side, "paf_row_sha256": row.paf_row_sha256,
                    "raw_core": context.raw_core, "raw_affected": context.raw_affected,
                    "affected_raw_sha256": sha256(raw.seq[slice(*context.raw_affected)]).hexdigest(),
                    "raw_footprint": context.raw_footprint, "graph_core": context.graph_core,
                    "graph_footprint": context.graph_footprint, "reference_core": reference_core,
                    "reference_affected": reference_affected, "independent_donor_molecules": donors,
                })
        return dict(verified)

    def summary(self):
        return {
            "binding": self._binding.binding(), "discovery_report": self._report_file.binding(),
            "discovery_evidence": self._evidence_file.binding(), "scan_k": self._k,
            "registered_rows": sum(map(len, self._registrations.values())),
            "registered_originals": len(self._registrations), "verified_originals": len(self._originals),
            "registered_paf_rows": len(self.paf_hashes), "verified_projections": self._verified_projections,
            "verified_with_independent_donor": self._donor_supported,
            "verified_without_independent_donor": self._verified_projections - self._donor_supported,
            "graph_input_originals": len(self._graph_originals),
        }
