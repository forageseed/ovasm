"""Audit same-original long-read geometry from a frozen, reported PAF.

This adapter emits a geometry report. It never treats a heuristic mapper's
missing hit as an excluded source and never writes a native admission manifest.
"""

import argparse
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from hashlib import sha256
from pathlib import Path
import re
import sys

from prepare_read_admission import LockedFile, fasta, fastq, open_input, read_json, write_json
from source_competition import Anchor, QueryScope, contains, parse_paf, verify_witness
from validate_source_exclusion import Abstain, Graph, enumerate_routes

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


def enrolled(sources, windows, panel):
    words = {w.id: (panel[w.reference][slice(*w.interval)],
                    tuple(panel[w.reference][slice(*w.guard(side))] for side in ("left", "right")))
             for w in windows.values()}
    selected = {}
    identities = {}
    query_ids = set()
    for source in sources:
        for raw in records(source):
            require(raw.query not in query_ids and raw.molecule not in identities,
                    "duplicate query ID or original molecule")
            query_ids.add(raw.query)
            identities[raw.molecule] = raw.binding()
            seq = raw.seq.upper()
            hits = [wid for wid, (word, guards) in words.items() if matches_window(seq, word, guards)]
            if hits:
                selected[raw.query] = raw, hits
    return selected, identities


def parse_alignments(paf, selected, panel):
    rows = defaultdict(set)
    seen = 0
    with open_input(paf) as stream:
        for line in stream:
            aln = parse_paf(line.decode("ascii"))
            require(aln.query in selected and aln.reference in panel, "PAF row outside named originals or panel")
            require(aln.query_length == len(selected[aln.query][0].seq)
                    and aln.reference_length == len(panel[aln.reference]), "PAF lengths differ from bound bytes")
            rows[aln.query].add(aln)
            seen += 1
    return {q: sorted(v, key=alignment_key) for q, v in rows.items()}, seen


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


def observe(w, raw, alignments, panel, roles, index):
    seq = raw.seq
    ref = panel[w.reference]
    role = roles[w.reference]
    exact_cores = set()
    witnesses = {}
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
        if w.require_full_query and aln.query_interval != (0, len(seq)):
            blocked["original_not_fully_aligned"] += 1
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
        roles_covering_scope = sorted({roles[a.reference] for a in alignments if a.continuous_over(scope.interval)})
        witnesses[proof] = asdict(proof) | {"cigar": "".join(
            f"{max(x.oriented_query[1]-x.oriented_query[0], x.reference[1]-x.reference[0])}{x.op}"
            for x in aln.steps),
            "reported_paf_roles_covering_same_scope": roles_covering_scope}
    distance = None
    if w.physical_boundaries:
        distance = min(w.interval[0] - w.physical_boundaries[0], w.physical_boundaries[1] - w.interval[1])
    return dict(exact_core_projections=[dict(raw_bed0=i, strand=s) for i, s in sorted(exact_cores)],
                witnesses=list(witnesses.values()), excluded_alignment_counts=dict(blocked),
                registered_exterior_distance_bp=distance,
                read_length_reaches_registered_exterior=(len(seq) > distance if distance is not None else None),
                competition_state="unresolved_reported_paf_only")


def graph_context(graph, graph_reads, windows, panel, originals):
    gfa = None
    try:
        with open_input(graph) as stream:
            gfa = Graph.parse(stream.read())
    except Abstain as error:
        note = str(error)
    donors = defaultdict(dict)
    graph_ids = set()
    windows = {key: w for key, w in windows.items() if w.graph_required}
    words = {key: panel[w.reference][slice(*w.interval)] for key, w in windows.items()}
    for f in graph_reads:
        for raw in records(f):
            require(raw.molecule not in graph_ids, "graph input repeats an original molecule")
            graph_ids.add(raw.molecule)
            identity = raw.binding()
            if raw.molecule in originals:
                require(identity == originals[raw.molecule], "graph input bytes differ for original molecule")
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


def report(plan_path):
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
    require(set(plan) == required and isinstance(plan["scope"], str) and plan["scope"].strip(),
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
    selected, identities = enrolled(originals, windows, panel)
    query_file = lock(plan["mapper_queries"])
    with open_input(query_file) as stream:
        query_empty = not stream.read(1)
    query_input = {} if query_empty else {k: bytes(v) for k, v in fasta(query_file)}
    require(query_input.keys() == selected.keys() and all(query_input[k] == v[0].seq.upper()
            for k, v in selected.items()), "mapper query bytes or enrollment differ from original FASTQs")
    alignments, paf_rows = parse_alignments(lock(plan["paf"]["data"]), selected, panel)
    lock(plan["paf"]["provenance"])
    lock(graph["provenance"])
    routes, donors, graph_ids = graph_context(lock(graph["gfa"]), [lock(x) for x in graph["reads"]],
                                              windows, panel, identities)
    lock(plan["window_provenance"])
    index = AnchorIndex(windows, selected, alignments, panel, roles)
    rows = []
    counts = defaultdict(Counter)
    for query, (raw, ids) in sorted(selected.items()):
        results = {}
        for key in ids:
            w = windows[key]
            result = observe(w, raw, alignments.get(query, []), panel, roles, index)
            result["graph_input_original"] = raw.molecule in graph_ids
            result["other_graph_input_exact_donors"] = sorted(set(donors[key]) - {raw.molecule})
            route = routes.get(key, {})
            result["graph_spelling_unique"] = bool(route.get("complete") and len(route["routes"]) == 1)
            role = roles[w.reference]
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
    for artifact in locked_inputs:
        artifact.check()
    implementation = {}
    for module in (sys.modules[LockedFile.__module__], sys.modules[Graph.__module__],
                   sys.modules[verify_witness.__module__], sys.modules[__name__]):
        path = Path(module.__file__).resolve()
        implementation[path.name] = dict(path=str(path), sha256=sha256(path.read_bytes()).hexdigest())
    panel_summary = {role: dict(records=sum(r == role for r in roles.values()),
                                bp=sum(len(panel[name]) for name, r in roles.items() if r == role))
                     for role in ROLES}
    return dict(schema="ovasm.long-read-geometry-report.v1", plan=raw_plan.binding(), scope=plan["scope"],
                competition_scope="Reported PAF geometry and literal external anchors in the named panel.",
                competition_complete=False, native_admission_created=False, paf_rows=paf_rows,
                original_molecules_scanned=len(identities), original_molecules_enrolled=len(selected),
                graph_input_molecules=len(graph_ids), graph_routes=routes,
                graph_input_donors=donors, counts=counts, records=rows,
                panel_records_by_role=panel_summary, panel_component_inventory=component_inventory,
                implementation=implementation, immutable_inputs=[x.binding() for x in locked_inputs])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    try:
        require(not args.out.exists(), "output report already exists")
        result = report(args.plan)
        write_json(args.out, result)
    except (OSError, ValueError) as error:
        parser.exit(1, f"long-read geometry audit: {error}\n")
    print(f'{result["original_molecules_scanned"]} originals, '
          f'{result["original_molecules_enrolled"]} enrolled; competition unresolved')


if __name__ == "__main__":
    main()
