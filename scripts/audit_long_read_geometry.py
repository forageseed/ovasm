"""Audit long-read geometry and optionally verify whole-read edit competition.

A verified target CIGAR can supply a bound. Complete competition at that bound
comes from searching the named reference bytes, independently of mapper hits.
"""

import argparse
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path
import re
import sys

from prepare_read_admission import LockedFile, fasta, fastq, open_input, read_json
from source_competition import Anchor, QueryScope, contains, parse_paf, verify_witness
from validate_source_exclusion import Abstain, Graph, enumerate_routes, new_output, publish_outputs

ROLES = ("mitochondrion", "plastid", "nuclear")
BITS = {role: 1 << i for i, role in enumerate(ROLES)}
COMP = bytes.maketrans(b"ACGTN", b"TGCAN")
BASE_CODES = bytes.maketrans(b"ACGT", bytes(range(4)))
DNA = re.compile(rb"[ACGT]+")
MOVIE = re.compile(rb"(?<!\S)(m[^/\s]+/[0-9]+)(?:/ccs(?:/[^\s]*)?)?(?!\S)")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def reverse(seq):
    return seq.translate(COMP)[::-1]


def canonical(seq):
    return min(seq, reverse(seq))


def find(seq, word):
    start = seq.find(word)
    while start >= 0:
        yield start
        start = seq.find(word, start + 1)


def word_hits(seq, word):
    rev = reverse(word)
    return [(a, a + len(word), strand)
            for strand, pattern in (("+", word), ("-", rev))
            for a in find(seq, pattern)]


def matches_window(seq, word, guards):
    if word in seq or reverse(word) in seq:
        return True
    left, right = (word_hits(seq, g) for g in guards)
    return any(a[2] == b[2] and (a[1] <= b[0] if a[2] == "+" else b[1] <= a[0])
               for a in left for b in right)


@dataclass(frozen=True)
class Window:
    id: str
    reference: str
    interval: tuple
    anchor_bp: int
    mode: str
    exclusive_sides: tuple
    graph_required: bool
    physical_boundaries: tuple
    require_full_query: bool

    def guard(self, side):
        lo, hi = self.interval
        return (lo - self.anchor_bp, lo) if side == "left" else (hi, hi + self.anchor_bp)

    def graph_interval(self):
        if self.mode == "adjacent_boundary":
            return self.guard("left")[0], self.guard("right")[1]
        return self.interval

    @classmethod
    def load(cls, obj, panel):
        require(isinstance(obj, dict) and set(obj) == {"id", "reference", "interval", "anchor_bp", "mode", "exclusive_sides",
                             "graph_required", "physical_boundaries", "require_full_query"},
                "invalid window keys")
        window = cls(**(obj | {k: tuple(obj[k]) for k in ("interval", "exclusive_sides", "physical_boundaries")}))
        require(isinstance(window.id, str) and window.id and window.reference in panel,
                "invalid window identity")
        require(type(window.anchor_bp) is int and window.anchor_bp > 0, "invalid anchor width")
        require(len(window.interval) == 2 and all(type(x) is int for x in window.interval)
                and contains((0, len(panel[window.reference])), window.interval), "invalid core interval")
        require(window.mode in ("exact_context", "adjacent_boundary", "internal_control"), "invalid window mode")
        require(set(window.exclusive_sides) <= {"left", "right"}
                and len(set(window.exclusive_sides)) == len(window.exclusive_sides), "invalid exclusive sides")
        if window.mode != "internal_control":
            require(window.exclusive_sides, "source geometry needs a declared exclusive side")
        require(type(window.graph_required) is bool and type(window.require_full_query) is bool,
                "invalid window requirement")
        for side in ("left", "right"):
            interval = window.guard(side)
            require(contains((0, len(panel[window.reference])), interval)
                    and DNA.fullmatch(panel[window.reference][slice(*interval)]), "invalid enrollment guard")
        require(DNA.fullmatch(panel[window.reference][slice(*window.interval)]), "ambiguous core")
        if window.physical_boundaries:
            require(window.mode == "internal_control" and len(window.physical_boundaries) == 2
                    and all(type(x) is int for x in window.physical_boundaries)
                    and 0 <= window.physical_boundaries[0] < window.interval[0]
                    < window.interval[1] < window.physical_boundaries[1] <= len(panel[window.reference]),
                    "invalid registered exterior boundary bounds")
        return window


@dataclass(frozen=True)
class Original:
    query: str
    molecule: str
    header: bytes
    seq: bytes
    quality: bytes

    def binding(self):
        return {"query": self.query, "molecule": self.molecule, "read_id_hex": self.header.hex(),
                "sequence_sha256": sha256(self.seq).hexdigest(),
                "quality_sha256": sha256(self.quality).hexdigest()}


def records(source):
    for rid, seq, qual in fastq(source):
        m = MOVIE.findall(rid)
        require(len(m) == 1, "one original PacBio movie/ZMW identity is required per record")
        require(rid.split(), "original FASTQ header has no query name")
        yield Original(rid.split()[0].decode("ascii"), m[0].decode("ascii"), rid, seq, qual)


def enrolled(sources, windows, panel, query_input, frontier):
    words = {w.id: (panel[w.reference][slice(*w.interval)],
                    tuple(panel[w.reference][slice(*w.guard(side))] for side in ("left", "right")))
             for w in windows.values()}
    selected = {}
    identities = {}
    query_ids = set()
    for input_index, source in enumerate(sources):
        for ordinal, raw in enumerate(records(source)):
            require(raw.query not in query_ids and raw.molecule not in identities,
                    "duplicate query ID or original molecule")
            query_ids.add(raw.query)
            identities[raw.molecule] = raw.binding()
            seq = raw.seq.upper()
            if raw.query in query_input:
                require(query_input[raw.query] == seq, "mapper query bytes or enrollment differ from original FASTQs")
            if frontier is not None:
                frontier.record_original(input_index, ordinal, raw)
                if raw.query not in frontier.queries:
                    continue
            hits = [wid for wid, (word, guards) in words.items() if matches_window(seq, word, guards)]
            if hits:
                selected[raw.query] = raw, hits
    if frontier is None:
        require(query_input.keys() == selected.keys(), "mapper query bytes or enrollment differ from original FASTQs")
    else:
        frontier.finish_originals()
        require(selected.keys() == frontier.queries and selected.keys() <= query_input.keys() <= query_ids,
                "frontier enrollment and mapper query names differ from bound originals")
    return selected, identities


def parse_alignments(paf, query_input, selected, panel, requested_hashes):
    rows = defaultdict(set)
    named = {}
    seen = ignored = 0
    with open_input(paf) as stream:
        for line in stream:
            aln = parse_paf(line.decode("ascii"))
            require(aln.query in query_input and aln.reference in panel, "PAF row outside named originals or panel")
            require(aln.query_length == len(query_input[aln.query])
                    and aln.reference_length == len(panel[aln.reference]), "PAF lengths differ from bound bytes")
            seen += 1
            if aln.query in selected:
                rows[aln.query].add(aln)
            else:
                ignored += 1
            if requested_hashes:
                digest = sha256(line).hexdigest()
                if digest in requested_hashes:
                    named[digest] = aln
    return {q: sorted(v, key=alignment_key) for q, v in rows.items()}, named, seen, ignored


def alignment_key(aln):
    return (aln.reference, aln.strand, aln.query_interval, aln.reference_interval,
            tuple((s.oriented_query, s.reference, s.op) for s in aln.steps))


def reference_projection(aln, interval):
    """Inverse of a gapless aligned interval, in forward raw coordinates."""
    lo, hi = interval
    reached = lo
    start = end = None
    for step in aln.steps:
        t0, t1 = step.reference
        a, b = max(lo, t0), min(hi, t1)
        if a >= b:
            continue
        if step.op not in "M=" or a != reached:
            return None
        qa = step.oriented_query[0] + a - t0
        qb = qa + b - a
        if end is not None and end != qa:
            return None
        if start is None:
            start = qa
        reached, end = b, qb
        if reached == hi:
            return aln.orient((start, end))
    return None


def exact_anchor(aln, seq, ref, interval):
    raw = reference_projection(aln, interval)
    if raw is not None and aln.project_exact(raw, seq, ref) == interval:
        return Anchor(raw, interval)
    return None


def packed_codes(codes):
    value = 0
    for base in codes:
        value = (value << 2) | base
    return value


def packed_words(words_by_width):
    """Map forward and reverse-complement encodings back to canonical words."""
    return {width: {packed_codes(orientation.translate(BASE_CODES)): word
                    for word in words for orientation in (word, reverse(word))}
            for width, words in words_by_width.items()}


def exact_word_counts(seq, words_by_width):
    """Count overlapping matches on a single linear record, once per physical position."""
    counts = Counter()
    if not words_by_width:
        return counts
    codes = memoryview(seq.translate(BASE_CODES))
    runs = [m.span() for m in DNA.finditer(seq)]
    for width, words in words_by_width.items():
        mask = (1 << (2 * width)) - 1
        for start, end in runs:
            if end - start < width:
                continue
            prefix_end = start + width - 1
            value = packed_codes(codes[start:prefix_end])
            for base in codes[prefix_end:end]:
                value = ((value << 2) | base) & mask
                word = words.get(value)
                if word is not None:
                    counts[word] += 1
    return counts


class AnchorIndex:
    """Literal role membership, plus occurrence count on each named source record."""

    def __init__(self, windows, selected, alignments, panel, roles):
        ranges = defaultdict(list)
        for w in windows.values():
            for side in w.exclusive_sides:
                if w.mode == "adjacent_boundary":
                    ranges[w.reference, w.anchor_bp].append(w.guard(side))
            if w.mode != "exact_context":
                continue
            relevant = [a.reference_interval for q, (_, ids) in selected.items() if w.id in ids
                        for a in alignments.get(q, ()) if a.reference == w.reference
                        and contains(a.reference_interval, w.interval)]
            if relevant:
                ranges[w.reference, w.anchor_bp].extend(
                    [(min(x[0] for x in relevant), w.interval[0]),
                     (w.interval[1], max(x[1] for x in relevant))])
        words = defaultdict(set)
        for (reference, k), intervals in ranges.items():
            seq = panel[reference]
            for lo, hi in intervals:
                for i in range(lo, hi - k + 1):
                    word = seq[i:i+k]
                    if DNA.fullmatch(word):
                        words[k].add(canonical(word))
        self.role_bits = {word: 0 for keys in words.values() for word in keys}
        self.record_counts = {ref: Counter() for ref, _ in ranges}
        encoded = packed_words(words)
        for ref, seq in panel.items():
            matches = exact_word_counts(seq, encoded)
            for word in matches:
                self.role_bits[word] |= BITS[roles[ref]]
            if ref in self.record_counts:
                self.record_counts[ref].update(matches)

    def exclusive(self, word, reference, role):
        key = canonical(word)
        return self.role_bits.get(key) == BITS[role] and self.record_counts[reference][key] == 1


def nearest_anchor(aln, seq, ref, boundary, side, k, exclusive):
    oriented = seq.upper() if aln.strand == "+" else reverse(seq.upper())
    lo, hi = aln.reference_interval
    diagonals = []
    for step in aln.steps:
        if step.op not in "M=":
            continue
        t0, t1 = step.reference
        q0, q1 = step.oriented_query
        if diagonals and diagonals[-1][1] == t0 and diagonals[-1][3] == q0:
            old_t0, _, old_q0, _ = diagonals[-1]
            diagonals[-1] = old_t0, t1, old_q0, q1
        else:
            diagonals.append((t0, t1, q0, q1))
    steps = reversed(diagonals) if side == "left" else diagonals
    for t0, t1, q0, _ in steps:
        if side == "left":
            positions = range(min(t1, boundary) - k, max(t0, lo) - 1, -1)
        else:
            positions = range(max(t0, boundary), min(t1, hi) - k + 1)
        for t in positions:
            q = q0 + t - t0
            word = ref[t:t+k]
            if oriented[q:q+k] == word and exclusive(word):
                return Anchor(aln.orient((q, q+k)), (t, t+k))
    return None


def cigar_text(alignment):
    return "".join(f"{max(step.oriented_query[1] - step.oriented_query[0], step.reference[1] - step.reference[0])}{step.op}"
                   for step in alignment.steps)


def observe(w, raw, alignments, panel, roles, index, verified_full=None, collect_partial_witnesses=False):
    seq = raw.seq
    ref = panel[w.reference]
    role = roles[w.reference]
    exact_cores = set()
    witnesses = {}
    partial_witnesses = {}
    blocked = Counter()
    for aln in alignments:
        if aln.reference != w.reference or not contains(aln.reference_interval, w.interval):
            continue
        core = exact_anchor(aln, seq, ref, w.interval)
        core_exact = core is not None
        if core_exact:
            exact_cores.add((core.raw, aln.strand))
        if w.mode == "internal_control":
            continue
        partial_query = w.require_full_query and aln.query_interval != (0, len(seq))
        if partial_query:
            blocked["original_not_fully_aligned"] += 1
            if not collect_partial_witnesses:
                continue
        if w.mode == "exact_context" and core is None:
            blocked["exact_core_missing"] += 1
            continue
        def exclusive(word):
            return index.exclusive(word, w.reference, role)

        anchors = []
        for side, boundary in zip(("left", "right"), w.interval):
            is_exclusive_side = side in w.exclusive_sides
            if w.mode == "exact_context":
                hit = nearest_anchor(aln, seq, ref, boundary, side, w.anchor_bp,
                                     exclusive if is_exclusive_side else DNA.fullmatch)
            else:
                guard = w.guard(side)
                hit = exact_anchor(aln, seq, ref, guard)
                if is_exclusive_side and not exclusive(ref[slice(*guard)]):
                    hit = None
            anchors.append(hit)
        l, r = anchors
        if l is None or r is None:
            blocked["exterior_anchor_missing"] += 1
            continue
        if core is None:
            a, b = sorted((l.raw, r.raw))
            core = Anchor((a[1], b[0]), w.interval)
        scope = QueryScope.bind(raw.query, seq, aln.query_interval if w.require_full_query
                                else (min(l.raw[0], r.raw[0]), max(l.raw[1], r.raw[1])))
        proof = verify_witness(scope, aln, seq, ref, role, core, l, r,
                               exact_boundary=core_exact)
        if proof is None:
            blocked["single_cigar_geometry_incompatible"] += 1
            continue
        if partial_query:
            partial_witnesses[proof] = asdict(proof) | {"cigar": cigar_text(aln)}
            continue
        if (verified_full is not None and aln.query_interval == (0, len(seq))
                and all(step.op != "N" for step in aln.steps)):
            # The two bases straddling each junction must be exact and adjacent.
            # This permits edits inside the interval while witnessing both joins.
            if w.mode == "exact_context" or all(
                    exact_anchor(aln, seq, ref, (pos - 1, pos + 1)) for pos in w.interval):
                full_proof = proof
                if scope.interval != aln.query_interval:
                    full_scope = QueryScope.bind(raw.query, seq, aln.query_interval)
                    full_proof = verify_witness(full_scope, aln, seq, ref, role, core, l, r,
                                                exact_boundary=core_exact)
                if full_proof is not None:
                    verified_full.append((aln, full_proof))
            else:
                blocked["exact_junction_missing"] += 1
        roles_covering_scope = sorted({roles[a.reference] for a in alignments if a.continuous_over(scope.interval)})
        witnesses[proof] = asdict(proof) | {"cigar": cigar_text(aln),
            "reported_paf_roles_covering_same_scope": roles_covering_scope}
    distance = None
    if w.physical_boundaries:
        distance = min(w.interval[0] - w.physical_boundaries[0], w.physical_boundaries[1] - w.interval[1])
    result = dict(exact_core_projections=[dict(raw_bed0=i, strand=s) for i, s in sorted(exact_cores)],
                  witnesses=list(witnesses.values()), excluded_alignment_counts=dict(blocked),
                  registered_exterior_distance_bp=distance,
                  read_length_reaches_registered_exterior=(len(seq) > distance if distance is not None else None),
                  competition_state="unresolved_reported_paf_only")
    if partial_witnesses:
        result["partial_witnesses_for_source_exclusion"] = list(partial_witnesses.values())
    return result


def graph_context(graph, graph_reads, windows, panel, originals, frontier):
    gfa = None
    try:
        with open_input(graph) as stream:
            gfa = Graph.parse(stream.read())
    except Abstain as error:
        note = str(error)
    donors = defaultdict(dict)
    graph_ids = set()
    windows = {key: w for key, w in windows.items() if w.graph_required}
    words = {key: panel[w.reference][slice(*w.graph_interval())] for key, w in windows.items()}
    for f in graph_reads:
        for raw in records(f):
            require(raw.molecule not in graph_ids, "graph input repeats an original molecule")
            graph_ids.add(raw.molecule)
            identity = raw.binding()
            if raw.molecule in originals:
                require(identity == originals[raw.molecule], "graph input bytes differ for original molecule")
            if frontier is not None:
                frontier.record_graph_original(raw)
            seq = raw.seq.upper()
            for key, word in words.items():
                if word in seq or reverse(word) in seq:
                    donors[key][raw.molecule] = identity
    routes = {key: (enumerate_routes(gfa, word.decode("ascii"), 200_000, 4096) if gfa is not None
                    else dict(complete=False, reason=note, routes=None)) for key, word in words.items()}
    return routes, donors, graph_ids


def load_named_panel(spec, lock):
    require(isinstance(spec, dict) and set(spec) == {"fasta", "records"}, "invalid named panel fields")
    panel_fasta = lock(spec["fasta"])
    panel = {k: bytes(v) for k, v in fasta(panel_fasta)}
    require(isinstance(spec["records"], dict) and spec["records"].keys() == panel.keys(),
            "named panel FASTA and its record bindings differ")
    roles = {}
    components = defaultdict(lambda: defaultdict(list))
    component_inventory = []
    for k, record in spec["records"].items():
        require(isinstance(record, dict) and set(record) == {"role", "source", "source_record"}
                and isinstance(record["source_record"], str) and record["source_record"]
                and record["role"] in ROLES, "panel record needs a role, source FASTA and original record ID")
        roles[k] = record["role"]
        source = record["source"]
        require(isinstance(source, dict) and set(source) == {"path", "sha256"}
                and all(isinstance(x, str) for x in source.values()),
                "panel component requires path and sha256")
        components[source["path"], source["sha256"]][record["source_record"]].append(k)
    require(set(roles.values()) == set(ROLES), "panel roles are not complete")
    for (path, digest), wanted in components.items():
        verified = set()
        records_in_component = 0
        component = lock(dict(path=path, sha256=digest))
        for source_name, sequence in fasta(component):
            records_in_component += 1
            if source_name not in wanted:
                continue
            for target_name in wanted[source_name]:
                require(sequence == panel[target_name], "combined panel bytes differ from component")
            verified.add(source_name)
        require(verified == wanted.keys(), "source FASTA record missing from named panel")
        component_inventory.append(component.binding() | dict(records_in_component=records_in_component,
                                                               records_in_named_panel=len(verified)))
    return panel, roles, component_inventory


def assess_full_edits(rows, selected, observed, windows, panel, roles, budget, frontier_proofs,
                      exclusion_target=None):
    """Attach scoped whole-read competition to each geometry result."""
    from long_read_edit import cigar_edit_upper_bound, search_competitors

    row_key = "full_edit_competition" if exclusion_target is None else "full_edit_source_exclusion"
    role_key = "target" if exclusion_target is None else "source_role"
    bound_key = "verified_target_upper_bound" if exclusion_target is None else "verified_source_upper_bound"
    by_target = defaultdict(dict)
    decisions = {}
    for row in rows:
        query = row["query"]
        raw = selected[query][0]
        own = observed.get(query, {})
        targets = {roles[windows[key].reference] for key in own}
        if not targets:
            if exclusion_target is not None:
                status = "unresolved_no_full_competing_geometry"
            else:
                status = ("unresolved_no_verified_frontier_target_geometry" if frontier_proofs is not None
                          else "unresolved_no_full_target_geometry")
            decisions[query] = {"status": status}
            continue
        if len(targets) > 1:
            decisions[query] = {"status": "unresolved_conflicting_geometry"}
            continue
        target = next(iter(targets))
        observations = row["windows"]
        if any(roles[windows[key].reference] != target
               and (result["witnesses"] or (exclusion_target is not None
                   and result.get("partial_witnesses_for_source_exclusion")))
               for key, result in observations.items()):
            decisions[query] = {"status": "unresolved_conflicting_geometry", role_key: target}
            continue
        supported = {key: alignments for key, alignments in own.items()
                     if (exclusion_target is not None and not windows[key].graph_required)
                     or (observations[key]["graph_spelling_unique"]
                         and observations[key]["other_graph_input_exact_donors"])}
        if not supported:
            decisions[query] = {"status": "unresolved_no_independent_graph_support", role_key: target}
            continue
        sequence = raw.seq.upper()
        if not DNA.fullmatch(sequence):
            decisions[query] = {"status": "unresolved_unobservable", role_key: target}
            continue
        radius, key, alignment, proof = min(
            ((cigar_edit_upper_bound(aln, sequence, panel[aln.reference]), key, aln, witness)
             for key, alignments in supported.items() for aln, witness in alignments),
            key=lambda item: (item[0], item[1], alignment_key(item[2])),
        )
        decisions[query] = {
            role_key: target, "radius": radius,
            bound_key: {
                "window": key, "reference": alignment.reference, "strand": alignment.strand,
                "interval": alignment.reference_interval, "unit_edit_count": radius,
                "cigar": cigar_text(alignment), "witness": asdict(proof),
            },
        }
        if frontier_proofs is not None:
            decisions[query][bound_key]["frontier_events"] = frontier_proofs[query, key, alignment]
        by_target[target][query] = (sequence, radius)

    rows_by_query = {row["query"]: row for row in rows}
    for target, queries in by_target.items():
        competitors = {name: sequence for name, sequence in panel.items() if roles[name] != target}
        results = search_competitors(queries, competitors, {name: roles[name] for name in competitors}, budget)
        require(results.keys() == queries.keys(), "incomplete full-edit search result")
        for query, result in results.items():
            require(set(result) == set(ROLES) - {target}, "incomplete full-edit comparison roles")
            decision = decisions[query]
            decision["competing_roles"] = {role: asdict(value) for role, value in result.items()}
            if any(value.best_hit is not None for value in result.values()):
                status = "competing_within_edit_bound"
            elif any(value.status == "budget_exhausted" for value in result.values()):
                status = "unresolved_edit_budget"
            elif any(value.status != "complete" for value in result.values()):
                status = "unresolved_unobservable"
            elif (exclusion_target is None
                  and rows_by_query[query]["windows"][decision[bound_key]["window"]]["graph_input_original"]):
                status = "already_graph_input"
            else:
                status = "qualified_within_edit_model"
            decision["status"] = status

    for row in rows:
        row[row_key] = decisions[row["query"]]
    counts = Counter(decision["status"] for decision in decisions.values())
    result = {
        "scope": "full-query unit Levenshtein to contiguous substrings of named literal competing references",
        "budget": asdict(budget),
        "model_searches": sum(map(len, by_target.values())),
        "decision_counts": counts,
        "eligible_records": counts["qualified_within_edit_model"],
    }
    if exclusion_target is None:
        result["target_bound"] = "unit edits in one verified whole-original target alignment for a graph-supported window"
    else:
        result["excluded_target"] = exclusion_target
        result["source_bound"] = "unit edits in one verified whole-original alignment with role-exclusive bilateral anchors"
    return result


def admission_records(result):
    rows = [row for row in result["records"]
            if row["full_edit_competition"]["status"] == "qualified_within_edit_model"]
    if "frontier_event_binding" in result:
        require(all(row["full_edit_competition"]["verified_target_upper_bound"].get("frontier_events")
                    for row in rows), "frontier-qualified original is missing its event proof")
    # A stable alphanumeric query order keeps appended records reproducible across streams.
    def query_order(row):
        name = row["query"]
        return tuple(int(p) if i % 2 else p for i, p in enumerate(re.split(r"([0-9]+)", name))), name

    return [{key: row[key] for key in ("read_id_hex", "sequence_sha256", "quality_sha256")}
            | {"target": row["full_edit_competition"]["target"]}
            for row in sorted(rows, key=query_order)]


def exclusion_records(result):
    target = result["full_edit_source_exclusions"]["excluded_target"]
    return [{key: row[key] for key in ("read_id_hex", "sequence_sha256", "quality_sha256")} | {"target": target}
            for row in result["records"]
            if row["full_edit_source_exclusion"]["status"] == "qualified_within_edit_model"]


def report(plan_path, edit_budget=None, exclusion_target=None):
    locked_inputs = []
    root = plan_path.resolve().parent

    def lock(spec):
        locked = LockedFile.lock(root, spec)
        locked_inputs.append(locked)
        return locked

    raw_plan = lock(dict(path=str(plan_path.resolve()),
                         sha256=sha256(plan_path.read_bytes()).hexdigest()))
    plan = read_json(raw_plan)
    require(isinstance(plan, dict) and plan.get("schema") == "ovasm.long-read-geometry.v1",
            "incorrect geometry plan schema")
    required = {"schema", "scope", "panel", "windows", "original_fastqs", "mapper_queries",
                "paf", "graph", "window_provenance"}
    require(set(plan) - {"frontier_events"} == required and isinstance(plan["scope"], str) and plan["scope"].strip(),
            "invalid geometry plan fields or empty scope")
    require(isinstance(plan["windows"], list) and plan["windows"]
            and isinstance(plan["original_fastqs"], list) and plan["original_fastqs"],
            "geometry plan requires registered windows and original FASTQs")
    require(isinstance(plan["paf"], dict) and set(plan["paf"]) == {"data", "provenance"},
            "PAF requires a data binding and mapper provenance")
    graph = plan["graph"]
    require(isinstance(graph, dict) and set(graph) == {"gfa", "reads", "provenance"}
            and isinstance(graph["reads"], list) and graph["reads"],
            "graph requires a GFA, original graph FASTQs and provenance")
    panel, roles, component_inventory = load_named_panel(plan["panel"], lock)
    windows = {}
    for obj in plan["windows"]:
        w = Window.load(obj, panel)
        require(w.id not in windows, "duplicate window ID")
        windows[w.id] = w
    originals = [lock(p) for p in plan["original_fastqs"]]
    frontier = graph_file = None
    if "frontier_events" in plan:
        from frontier_event_binding import FrontierBindings

        graph_file = lock(graph["gfa"])
        frontier = FrontierBindings.load(lock(plan["frontier_events"]), lock, graph_file, originals, windows)
    query_file = lock(plan["mapper_queries"])
    with open_input(query_file) as stream:
        query_empty = not stream.read(1)
    query_input = {} if query_empty else {k: bytes(v) for k, v in fasta(query_file)}
    selected, identities = enrolled(originals, windows, panel, query_input, frontier)
    alignments, named_paf, paf_rows, ignored_paf = parse_alignments(
        lock(plan["paf"]["data"]), query_input, selected, panel, frontier.paf_hashes if frontier is not None else ())
    lock(plan["paf"]["provenance"])
    lock(graph["provenance"])
    routes, donors, graph_ids = graph_context(graph_file or lock(graph["gfa"]), [lock(x) for x in graph["reads"]],
                                              windows, panel, identities, frontier)
    lock(plan["window_provenance"])
    frontier_proofs = frontier.verify_alignments(named_paf, panel, windows) if frontier is not None else None
    index = AnchorIndex(windows, selected, alignments, panel, roles)
    rows = []
    counts = defaultdict(Counter)
    observed = defaultdict(dict) if edit_budget else None
    competing_observed = defaultdict(dict) if exclusion_target is not None else None
    for query, (raw, ids) in sorted(selected.items()):
        results = {}
        for key in ids:
            w = windows[key]
            role = roles[w.reference]
            bilateral = w.mode in ("exact_context", "adjacent_boundary") and set(w.exclusive_sides) == {"left", "right"}
            admission_geometry = bilateral and observed is not None and w.require_full_query and role in ROLES[:2]
            exclusion_geometry = bilateral and competing_observed is not None and role != exclusion_target
            full = [] if admission_geometry or exclusion_geometry else None
            result = observe(w, raw, alignments.get(query, []), panel, roles, index, full,
                             collect_partial_witnesses=competing_observed is not None)
            if full and frontier_proofs is not None:
                full = [(aln, proof) for aln, proof in full if (query, key, aln) in frontier_proofs]
            if full and admission_geometry:
                observed[query][key] = full
            if full and exclusion_geometry:
                competing_observed[query][key] = full
            result["graph_input_original"] = raw.molecule in graph_ids
            result["other_graph_input_exact_donors"] = sorted(set(donors[key]) - {raw.molecule})
            route = routes.get(key, {})
            result["graph_spelling_unique"] = bool(route.get("complete") and len(route["routes"]) == 1)
            no_reported_conflict = any(x["reported_paf_roles_covering_same_scope"] == [role]
                                       for x in result["witnesses"])
            result["graph_corroborated_observation"] = bool(no_reported_conflict and result["graph_spelling_unique"]
                                                       and result["other_graph_input_exact_donors"])
            counts[key]["enrolled"] += 1
            counts[key]["exact_core_originals"] += bool(result["exact_core_projections"])
            counts[key]["dual_anchor_originals"] += bool(result["witnesses"])
            counts[key]["exact_dual_anchor_originals"] += any(w["exact_boundary"] for w in result["witnesses"])
            counts[key]["gapped_dual_anchor_originals"] += any(not w["exact_boundary"] for w in result["witnesses"])
            counts[key]["graph_corroborated_observations"] += result["graph_corroborated_observation"]
            counts[key]["new_graph_corroborated_observations"] += (result["graph_corroborated_observation"]
                                                              and raw.molecule not in graph_ids)
            results[key] = result
        rows.append(raw.binding() | {"windows": results})
    edit_report = (assess_full_edits(rows, selected, observed, windows, panel, roles, edit_budget, frontier_proofs)
                   if edit_budget else None)
    exclusion_report = (assess_full_edits(rows, selected, competing_observed, windows, panel, roles,
                                          edit_budget, frontier_proofs, exclusion_target)
                        if competing_observed is not None else None)
    for artifact in locked_inputs:
        artifact.check()
    implementation = {}
    modules = [sys.modules[LockedFile.__module__], sys.modules[Graph.__module__],
               sys.modules[verify_witness.__module__], sys.modules[__name__]]
    if edit_budget:
        modules.append(sys.modules["long_read_edit"])
    if frontier is not None:
        modules.extend(sys.modules[name] for name in ("frontier_event_binding", "discover_graph_frontier"))
    for module in modules:
        path = Path(module.__file__).resolve()
        implementation[path.name] = dict(path=str(path), sha256=sha256(path.read_bytes()).hexdigest())
    panel_summary = {role: dict(records=sum(r == role for r in roles.values()),
                                bp=sum(len(panel[name]) for name, r in roles.items() if r == role))
                     for role in ROLES}
    result = dict(schema="ovasm.long-read-geometry-report.v1", plan=raw_plan.binding(), scope=plan["scope"],
                competition_scope="Reported PAF geometry and literal external anchors in the named panel.",
                competition_complete=False, native_admission_created=False, paf_rows=paf_rows,
                original_molecules_scanned=len(identities), original_molecules_enrolled=len(selected),
                graph_input_molecules=len(graph_ids), graph_routes=routes,
                graph_input_donors=donors, counts=counts, records=rows,
                panel_records_by_role=panel_summary, panel_component_inventory=component_inventory,
                implementation=implementation, immutable_inputs=[x.binding() for x in locked_inputs])
    if frontier is not None:
        result["frontier_event_binding"] = frontier.summary() | {
            "ignored_mapper_queries": len(query_input) - len(selected),
            "ignored_paf_rows": ignored_paf,
        }
    if edit_report is not None:
        result["full_edit_competition"] = edit_report
    if exclusion_report is not None:
        result["full_edit_source_exclusions"] = exclusion_report
    return result


def json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path)
    parser.add_argument("--out", type=Path, required=True, help="New JSON audit path")
    parser.add_argument("--full-edit-competition", action="store_true",
                        help="Search named competitors within each verified full-read target edit bound")
    parser.add_argument("--out-admission", type=Path,
                        help="New native manifest path for qualifying graph-supported original records")
    parser.add_argument("--exclude-target", choices=ROLES[:2],
                        help="Organelle assignment prohibited by verified evidence of a competing physical source")
    parser.add_argument("--out-exclusions", type=Path,
                        help="New native manifest path for whole-read source exclusions")
    parser.add_argument("--max-edit-seed-hits", type=int, help="Exact seed hits per query and role (default 500000)")
    parser.add_argument("--max-edit-text-bases", type=int, help="Verification bases per query and role (default 2000000)")
    args = parser.parse_args()
    try:
        out = new_output(args.out)
        admission = new_output(args.out_admission) if args.out_admission else None
        exclusion = new_output(args.out_exclusions) if args.out_exclusions else None
        paths = [p for p in (out, admission, exclusion) if p is not None]
        require(len(set(paths)) == len(paths), "audit, admission and exclusion need different output paths")
        require((exclusion is None) == (args.exclude_target is None),
                "--out-exclusions and --exclude-target must be provided together")
        require(args.full_edit_competition or all(value is None for value in
                (admission, exclusion, args.max_edit_seed_hits, args.max_edit_text_bases)),
                "admission, exclusions and edit budgets require --full-edit-competition")
        budget = None
        if args.full_edit_competition:
            from long_read_edit import EditBudget
            require(all(value is None or value > 0 for value in
                        (args.max_edit_seed_hits, args.max_edit_text_bases)),
                    "edit budgets must be positive")
            budget = EditBudget(max_seed_hits=args.max_edit_seed_hits or 500_000,
                                max_text_bases=args.max_edit_text_bases or 2_000_000)
        result = report(args.plan, budget, args.exclude_target)
        admitted = admission_records(result) if admission else []
        excluded = exclusion_records(result) if exclusion else []
        result["native_admission_created"] = bool(admitted)
        if exclusion:
            result["native_exclusions_created"] = bool(excluded)
        audit_data = json_bytes(result)
        outputs = [(out, audit_data)]
        proof = {"path": str(out), "sha256": sha256(audit_data).hexdigest()}
        for path, records, schema in ((admission, admitted, "ovasm.read-admission.v1"),
                                      (exclusion, excluded, "ovasm.read-source-exclusions.v1")):
            if records:
                manifest = {"schema": schema, "scope": result["scope"],
                            "records": [row | {"provenance": proof} for row in records]}
                outputs.append((path, json_bytes(manifest)))
        publish_outputs(outputs)
    except (OSError, ValueError) as error:
        parser.exit(1, f"long-read geometry audit: {error}\n")
    summary = f'{result["original_molecules_scanned"]} originals, {result["original_molecules_enrolled"]} enrolled'
    if budget:
        model = result["full_edit_competition"]
        summary += f', {model["model_searches"]} full-edit searches, {model["eligible_records"]} eligible'
        if exclusion:
            summary += f', {len(excluded)} source exclusions'
    else:
        summary += ', PAF competition unresolved'
    print(summary)


if __name__ == "__main__":
    main()
