"""Complete, bounded literal edit competition for plant organelle long reads.

A result applies to the named, linear reference records and unit-cost edits.
Panel provenance, biological source interpretation, and output are caller work.
"""

from collections import defaultdict
from dataclasses import dataclass, field
import re

from source_competition import Interval, PafAlignment

ROLES = frozenset({"mitochondrion", "plastid", "nuclear"})
DNA = re.compile(rb"[ACGT]+")
BASE_CODES = bytes.maketrans(b"ACGT", bytes(range(4)))
COMPLEMENT = bytes.maketrans(b"ACGTRYSWKMBDHVN", b"TGCAYRSWMKVHDBN")
MAX_PREFIX = 31


def reverse_complement(sequence: bytes) -> bytes:
    return sequence.translate(COMPLEMENT)[::-1]


def cigar_edit_upper_bound(alignment: PafAlignment, raw: bytes, reference: bytes) -> int:
    """Count literal substitutions and indel bases on one full original CIGAR."""
    if (alignment.query_interval != (0, len(raw)) or alignment.query_length != len(raw)
            or alignment.reference_length != len(reference)):
        raise ValueError("CIGAR must align the full original query against the named reference")
    if alignment.strand not in ("+", "-") or not alignment.steps:
        raise ValueError("invalid CIGAR strand or empty alignment")
    query = raw.upper()
    if alignment.strand == "-":
        query = reverse_complement(query)
    qpos = 0
    rpos, rend = alignment.reference_interval
    if not 0 <= rpos < rend <= len(reference):
        raise ValueError("invalid CIGAR reference interval")
    edits = 0
    for step in alignment.steps:
        a, b = step.oriented_query
        c, d = step.reference
        if (a, c) != (qpos, rpos) or not (a <= b <= len(query) and c <= d <= rend):
            raise ValueError("CIGAR does not describe a continuous edit path")
        if step.op in ("M", "=", "X") and b - a == d - c > 0:
            edits += sum(x != y for x, y in zip(query[a:b], reference[c:d].upper()))
        elif step.op == "I" and b > a and c == d:
            edits += b - a
        elif step.op == "D" and d > c and a == b:
            edits += d - c
        else:
            raise ValueError("CIGAR contains a skip or unsupported edit operation")
        qpos, rpos = b, d
    if (qpos, rpos) != (len(query), rend):
        raise ValueError("CIGAR does not consume the full declared alignment")
    return edits


def disjoint_blocks(query: bytes, radius: int) -> tuple[tuple[int, bytes], ...]:
    """Partition the oriented query into radius + 1 nonempty exact blocks."""
    if type(radius) is not int or not 0 <= radius < len(query) or not DNA.fullmatch(query):
        raise ValueError("nonempty ACGT query and radius smaller than its length are required")
    n = radius + 1
    edges = [i * len(query) // n for i in range(n + 1)]
    return tuple((lo, query[lo:hi]) for lo, hi in zip(edges, edges[1:]))


def candidate_interval(read_length: int, radius: int, offset: int, hit_start: int,
                       ref_length: int) -> Interval | None:
    """Bound a full-query alignment around one exact seed start."""
    if (any(type(x) is not int for x in (read_length, radius, offset, hit_start, ref_length))
            or not (0 <= radius < read_length and 0 <= offset < read_length and 0 <= hit_start < ref_length)):
        raise ValueError("invalid seed coordinates or edit radius")
    lo = max(0, hit_start - offset - radius)
    hi = min(ref_length, hit_start - offset + read_length + radius)
    return (lo, hi) if hi - lo >= read_length - radius else None


def _compile_pattern(pattern: bytes) -> tuple[int, dict[int, int]]:
    if not pattern:
        raise ValueError("edit pattern must be nonempty")
    masks = {}
    for i, base in enumerate(pattern):
        masks[base] = masks.get(base, 0) | (1 << i)
    return len(pattern), masks


def _scores(compiled, text, free_prefix):
    width, equals = compiled
    mask = (1 << width) - 1
    high_bit = 1 << (width - 1)
    positive, negative, distance = mask, 0, width
    first_bit = 0 if free_prefix else 1
    for char in text:
        equal = equals.get(char, 0)
        vertical = equal | negative
        horizontal = (((equal & positive) + positive) ^ positive) | equal
        plus = (negative | ~(horizontal | positive)) & mask
        minus = positive & horizontal
        distance += bool(plus & high_bit) - bool(minus & high_bit)
        plus = ((plus << 1) | first_bit) & mask
        minus = (minus << 1) & mask
        positive = (minus | ~(vertical | plus)) & mask
        negative = plus & vertical
        yield distance


def edit_scores(pattern: bytes, text: bytes, *, free_prefix: bool = True):
    """Yield full-pattern edit distances for every nonempty text endpoint."""
    return _scores(_compile_pattern(pattern), text, free_prefix)


def _merge(intervals):
    merged = []
    for lo, hi in sorted(intervals):
        if merged and lo <= merged[-1][1]:
            merged[-1] = merged[-1][0], max(merged[-1][1], hi)
        else:
            merged.append((lo, hi))
    return merged


@dataclass(frozen=True)
class EditBudget:
    """Per-query, per-role limits. The optional panel limit counts reference bases."""

    max_seed_hits: int
    max_text_bases: int
    max_panel_bases: int | None = None

    def __post_init__(self):
        for limit in (self.max_seed_hits, self.max_text_bases):
            if type(limit) is not int or limit < 0:
                raise ValueError("edit search budgets must be nonnegative integers")
        limit = self.max_panel_bases
        if limit is not None and (type(limit) is not int or limit < 0):
            raise ValueError("panel base budget must be a nonnegative integer")


@dataclass(frozen=True)
class EditHit:
    reference: str
    strand: str
    interval: Interval
    edit_distance: int


@dataclass(frozen=True)
class RoleResult:
    """A complete role either has its best hit within radius or no qualifying hit."""

    status: str
    reason: str | None
    panel_bases: int
    seed_hits: int
    windows_examined: int
    examined_bases: int
    best_hit: EditHit | None


@dataclass(frozen=True)
class _Seed:
    query_id: str
    strand: str
    offset: int
    read_length: int
    radius: int
    word: bytes


@dataclass
class _Work:
    panel_bases: int
    status: str = "complete"
    reason: str | None = None
    seed_hits: int = 0
    windows: dict = field(default_factory=lambda: defaultdict(set))

    def exhaust(self, reason: str):
        self.status = "budget_exhausted"
        self.reason = reason


def _packed(codes):
    value = 0
    for base in codes:
        value = (value << 2) | base
    return value


def _index_seeds(queries, oriented):
    indexes = defaultdict(lambda: defaultdict(list))
    for query_id, (query, radius) in queries.items():
        if query_id not in oriented:
            continue
        for strand, word in oriented[query_id].items():
            for offset, block in disjoint_blocks(word, radius):
                width = min(len(block), MAX_PREFIX)
                code = _packed(block[:width].translate(BASE_CODES))
                indexes[width][code].append(_Seed(query_id, strand, offset, len(query), radius, block))
    return indexes


def _seed_hits(sequence, indexes):
    """Enumerate overlapping matches to supplied oriented blocks."""
    if not indexes:
        return
    codes = memoryview(sequence.translate(BASE_CODES))
    for width, prefixes in sorted(indexes.items()):
        mask = (1 << (2 * width)) - 1
        for run in DNA.finditer(sequence):
            lo, hi = run.span()
            if hi - lo < width:
                continue
            start = lo + width - 1
            code = _packed(codes[lo:start])
            for i, base in enumerate(codes[start:hi], lo):
                code = ((code << 2) | base) & mask
                bucket = prefixes.get(code)
                if bucket is not None:
                    for seed in bucket:
                        if sequence.startswith(seed.word, i):
                            yield seed, i


def _seed_windows(panel, roles, indexes, work, budget):
    for name, raw in sorted(panel.items()):
        role = roles[name]
        if not any(c[role].status == "complete" for c in work.values()):
            continue
        sequence = raw.upper()
        for seed, start in _seed_hits(sequence, indexes):
            current = work[seed.query_id][role]
            if current.status != "complete":
                continue
            if current.seed_hits == budget.max_seed_hits:
                current.exhaust("seed_hits")
                current.windows.clear()
                continue
            current.seed_hits += 1
            interval = candidate_interval(seed.read_length, seed.radius, seed.offset, start, len(sequence))
            if interval:
                current.windows[seed.strand, name].add(interval)


def _search_windows(work, oriented, radius, panel, budget):
    best = None
    examined = windows_examined = 0
    compiled = {strand: _compile_pattern(query) for strand, query in oriented.items()}
    for (strand, ref), intervals in sorted(work.windows.items()):
        if work.status != "complete":
            break
        for lo, hi in _merge(intervals):
            if examined + hi - lo > budget.max_text_bases:
                work.exhaust("edit_text_bases")
                break
            sequence = panel[ref][lo:hi].upper()
            examined += len(sequence)
            windows_examined += 1
            for end, distance in enumerate(_scores(compiled[strand], sequence, True), 1):
                if distance <= radius and (best is None or distance < best[0]):
                    best = (distance, ref, strand, lo, lo + end)
    hit = None
    if best is not None:
        distance, ref, strand, lo, end = best
        query = oriented[strand]
        if distance == 0:
            hit = EditHit(ref, strand, (end - len(query), end), distance)
        else:
            full_length = min(end - lo, len(query) + radius)
            length = min(full_length, budget.max_text_bases - examined)
            text = panel[ref][end-length:end].upper()[::-1]
            for i, score in enumerate(edit_scores(query[::-1], text, free_prefix=False), 1):
                examined += 1
                if score == distance:
                    hit = EditHit(ref, strand, (end - i, end), distance)
                    break
            if hit is None:
                if length == full_length:
                    raise RuntimeError("semi-global edit result has no valid start")
                work.exhaust("edit_traceback_bases")
    return RoleResult(work.status, work.reason, work.panel_bases, work.seed_hits,
                      windows_examined, examined, hit)


def search_competitors(queries, panel, roles, budget: EditBudget) -> dict[str, dict[str, RoleResult]]:
    """Search both query orientations in all supplied named competitor records.

    Queries map identifier to (raw sequence, maximum edits). Panel and roles must
    identify the same nonempty set of records. A complete no-hit result excludes
    only alignments at or below that radius in those supplied linear references.
    """
    if (not panel or panel.keys() != roles.keys() or any(role not in ROLES for role in roles.values())
            or any(not isinstance(seq, (bytes, bytearray)) or not seq for seq in panel.values())):
        raise ValueError("supply nonempty named sequences and a valid role for every competitor record")
    if any(not isinstance(name, str) or not name for name in panel):
        raise ValueError("invalid competitor record name")
    role_sizes = defaultdict(int)
    for name, sequence in panel.items():
        role_sizes[roles[name]] += len(sequence)
    work = {}
    oriented = {}
    normalized = {}
    for query_id, (raw, radius) in queries.items():
        if not isinstance(query_id, str) or not query_id or not isinstance(raw, (bytes, bytearray)):
            raise ValueError("invalid query identifier or raw sequence")
        if type(radius) is not int or radius < 0:
            raise ValueError("edit radius must be a nonnegative integer")
        query = bytes(raw).upper()
        states = {role: _Work(size) for role, size in sorted(role_sizes.items())}
        work[query_id] = states
        if not DNA.fullmatch(query) or radius >= len(query):
            for state in states.values():
                state.status = "inapplicable"
                state.reason = "nonempty_acgt_partition_unavailable"
            continue
        for state in states.values():
            if budget.max_panel_bases is not None and state.panel_bases > budget.max_panel_bases:
                state.exhaust("panel_bases")
        if any(state.status == "complete" for state in states.values()):
            oriented[query_id] = {"+": query, "-": reverse_complement(query)}
            normalized[query_id] = query, radius
    indexes = _index_seeds(normalized, oriented)
    _seed_windows(panel, roles, indexes, work, budget)
    return {
        query_id: {
            role: (_search_windows(state, oriented[query_id], normalized[query_id][1], panel, budget)
                   if state.status == "complete" else RoleResult(
                       state.status, state.reason, state.panel_bases, state.seed_hits, 0, 0, None))
            for role, state in states.items()
        }
        for query_id, states in work.items()
    }
