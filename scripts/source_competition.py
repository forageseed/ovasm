"""Pure, conservative geometry and decision helpers for plant organelle source audits.

The caller binds original records, reference identities, and complete competing
searches. PAF alone does not certify an absent placement or a physical source.
All intervals are zero-based and half-open. Raw query coordinates are never
reported in reverse-complement space.
"""

from dataclasses import dataclass
from hashlib import sha256
import re

Interval = tuple[int, int]
ROLES = frozenset({"mitochondrion", "plastid", "nuclear"})
ORGANELLES = ROLES - {"nuclear"}
CIGAR_OP = re.compile(r"([1-9][0-9]*)([M=XIDN])")
DNA = frozenset(b"ACGT")
COMPLEMENT = bytes.maketrans(b"ACGTacgt", b"TGCAtgca")


def contains(outer: Interval, inner: Interval) -> bool:
    return outer[0] <= inner[0] < inner[1] <= outer[1]


@dataclass(frozen=True)
class CigarStep:
    oriented_query: Interval
    reference: Interval
    op: str


@dataclass(frozen=True)
class PafAlignment:
    query: str
    query_length: int
    query_interval: Interval
    strand: str
    reference: str
    reference_length: int
    reference_interval: Interval
    steps: tuple[CigarStep, ...]

    def orient(self, query_interval: Interval) -> Interval:
        if self.strand == "+":
            return query_interval
        start, end = query_interval
        return self.query_length - end, self.query_length - start

    def continuous_over(self, query_interval: Interval) -> bool:
        """Check that one CIGAR covers the span without a skipped reference region."""
        if not contains(self.query_interval, query_interval):
            return False
        start, end = self.orient(query_interval)
        return not any(step.op == "N" and start < step.oriented_query[0] < end
                       for step in self.steps)

    def projects_boundary(self, anchor: "Anchor") -> bool:
        """Verify both coordinate pairs and a nonempty aligned interior on one path."""
        q0, q1 = self.orient(anchor.raw)
        t0, t1 = anchor.reference
        start = end = interior = False
        for step in self.steps:
            a, b = step.oriented_query
            c, d = step.reference
            if step.op in "M=X":
                start |= a <= q0 <= b and t0 - c == q0 - a
                end |= a <= q1 <= b and t1 - c == q1 - a
                interior |= max(q0, a, a + t0 - c) < min(q1, b, a + t1 - c)
            elif step.op in "ID":
                start |= a <= q0 <= b and c <= t0 <= d
                end |= a <= q1 <= b and c <= t1 <= d
        return start and end and interior

    def project_exact(self, query_interval: Interval, raw: bytes, reference: bytes) -> Interval | None:
        """Project a gapless, literal A/C/G/T anchor on this one alignment."""
        if len(raw) != self.query_length or len(reference) != self.reference_length:
            raise ValueError("PAF lengths differ from original or reference")
        if not contains(self.query_interval, query_interval):
            return None
        start, end = self.orient(query_interval)
        covered = start
        ref_start = ref_end = None
        for step in self.steps:
            q0, q1 = step.oriented_query
            a, b = max(start, q0), min(end, q1)
            if a >= b:
                continue
            if step.op not in "M=" or a != covered:
                return None
            ta = step.reference[0] + a - q0
            tb = ta + b - a
            if ref_end is not None and ta != ref_end:
                return None
            if ref_start is None:
                ref_start = ta
            covered, ref_end = b, tb
            if covered == end:
                break
        if covered != end:
            return None
        raw_word = raw[slice(*query_interval)].upper()
        if self.strand == "-":
            raw_word = raw_word.translate(COMPLEMENT)[::-1]
        if not set(raw_word) <= DNA or raw_word != reference[ref_start:ref_end].upper():
            return None
        return ref_start, ref_end


def parse_paf(line: str) -> PafAlignment:
    """Read one PAF row with one cg:Z tag and consume its CIGAR exactly."""
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 12:
        raise ValueError("malformed PAF row")
    query, qlen, qa, qb, strand, ref, tlen, ta, tb, matches, columns, mapq = fields[:12]
    qlen, qa, qb, tlen, ta, tb, matches, columns, mapq = map(
        int, (qlen, qa, qb, tlen, ta, tb, matches, columns, mapq))
    if (not query or not ref or strand not in ("+", "-")
            or not contains((0, qlen), (qa, qb)) or not contains((0, tlen), (ta, tb))
            or not 0 <= matches <= min(columns, qb - qa, tb - ta)
            or columns <= 0 or not 0 <= mapq <= 255):
        raise ValueError("invalid PAF coordinates, strand or counts")
    tags = [tag for tag in fields[12:] if tag.startswith("cg:")]
    if len(tags) != 1 or not tags[0].startswith("cg:Z:"):
        raise ValueError("exactly one cg:Z PAF CIGAR is required")
    cigar = tags[0][5:]
    x, y = (qa if strand == "+" else qlen - qb), ta
    cursor = 0
    steps = []
    for match in CIGAR_OP.finditer(cigar):
        if match.start() != cursor:
            raise ValueError("malformed or unsupported PAF CIGAR")
        n, op = int(match[1]), match[2]
        qnext = x + n if op in "M=XI" else x
        tnext = y + n if op in "M=XDN" else y
        steps.append(CigarStep((x, qnext), (y, tnext), op))
        x, y, cursor = qnext, tnext, match.end()
    if not steps or cursor != len(cigar):
        raise ValueError("malformed or unsupported PAF CIGAR")
    if x != (qb if strand == "+" else qlen - qa) or y != tb:
        raise ValueError("PAF CIGAR does not consume the advertised coordinates")
    return PafAlignment(query, qlen, (qa, qb), strand, ref, tlen, (ta, tb), tuple(steps))


@dataclass(frozen=True)
class QueryScope:
    query: str
    read_sha256: str
    read_length: int
    interval: Interval

    @classmethod
    def bind(cls, query: str, raw: bytes, interval: Interval):
        if not query or not contains((0, len(raw)), interval):
            raise ValueError("invalid original query scope")
        return cls(query, sha256(raw).hexdigest(), len(raw), interval)


@dataclass(frozen=True, order=True)
class Placement:
    role: str
    reference: str
    boundary: Interval
    strand: str

    def __post_init__(self):
        if (self.role not in ROLES or not self.reference or self.strand not in ("+", "-")
                or not 0 <= self.boundary[0] < self.boundary[1]):
            raise ValueError("invalid physical placement")


@dataclass(frozen=True)
class Anchor:
    raw: Interval
    reference: Interval


@dataclass(frozen=True)
class Witness:
    scope: QueryScope
    placement: Placement
    boundary: Anchor
    left: Anchor
    right: Anchor
    exact_boundary: bool


def verify_witness(
    scope: QueryScope,
    alignment: PafAlignment,
    raw: bytes,
    reference: bytes,
    role: str,
    boundary: Anchor,
    left: Anchor,
    right: Anchor,
    *,
    exact_boundary: bool,
) -> Witness | None:
    """Check one same-original alignment across two exterior anchors and a boundary.

    Left and right are ordered on the forward reference, regardless of PAF strand.
    Exterior anchors are literal and gapless. An indel in between may be accepted
    when the caller has registered a gapped boundary policy. An N is never bridged.
    """
    if role not in ROLES:
        raise ValueError("unknown reference role")
    if (scope != QueryScope.bind(alignment.query, raw, scope.interval)
            or len(raw) != alignment.query_length or len(reference) != alignment.reference_length):
        raise ValueError("alignment, sequence and scope do not bind to the same original")
    if not alignment.continuous_over(scope.interval):
        return None
    points = left, boundary, right
    if any(not contains(scope.interval, p.raw)
           or not contains(alignment.reference_interval, p.reference) for p in points):
        return None
    for first, second in zip(points, points[1:]):
        if (first.reference[1] > second.reference[0]
                or alignment.orient(first.raw)[1] > alignment.orient(second.raw)[0]):
            return None
    if any(alignment.project_exact(p.raw, raw, reference) != p.reference for p in (left, right)):
        return None
    if exact_boundary:
        if alignment.project_exact(boundary.raw, raw, reference) != boundary.reference:
            return None
    elif not alignment.projects_boundary(boundary):
        return None
    placement = Placement(role, alignment.reference, boundary.reference, alignment.strand)
    return Witness(scope, placement, boundary, left, right, exact_boundary)


@dataclass(frozen=True)
class RoleSearch:
    """Caller-certified same-scope search, with plausible physical placements.

    complete requires a verified source-panel scope and a documented search policy.
    In particular, omission from a capped or heuristic PAF is not complete absence.
    """

    scope: QueryScope
    role: str
    complete: bool
    placements: tuple[Placement, ...]


@dataclass(frozen=True)
class Decision:
    state: str
    reason: str
    placement: Placement | None = None


def decide(
    scope: QueryScope,
    target_role: str,
    searches: tuple[RoleSearch, ...],
    witnesses: tuple[Witness, ...],
    *,
    aliases: dict[Placement, Placement] | None = None,
) -> Decision:
    """Require one witnessed physical placement over complete mt/PT/nuclear searches.

    Alias equivalence must be independently certified upstream. It identifies the
    same physical locus even when coordinate or strand differs across references.
    Unmapped aliases remain separate possible placements.
    """
    if target_role not in ORGANELLES:
        raise ValueError("target role must be a plant organelle")
    aliases = {} if aliases is None else aliases
    for a, b in aliases.items():
        if a == b or b in aliases or a.role != b.role:
            raise ValueError("aliases must point to one direct canonical placement of the same role")
    by_role = {}
    for search in searches:
        if search.role not in ROLES or search.role in by_role:
            raise ValueError("unexpected or repeated search role")
        if search.scope != scope:
            return Decision("unresolved", "competition_scope_mismatch")
        if any(p.role != search.role for p in search.placements):
            raise ValueError("search placement has the wrong role")
        by_role[search.role] = search
    if by_role.keys() != ROLES or any(s.complete is not True for s in by_role.values()):
        return Decision("unresolved", "incomplete_competition")
    if any(w.scope != scope for w in witnesses):
        return Decision("unresolved", "witness_scope_mismatch")
    placements = {aliases.get(p, p) for s in searches for p in s.placements}
    supported = {aliases.get(w.placement, w.placement) for w in witnesses}
    if not supported <= placements:
        return Decision("unresolved", "witness_outside_search")
    if len({p.role for p in placements}) > 1:
        return Decision("unresolved", "competing_roles")
    if len(placements) > 1:
        return Decision("unresolved", "multiple_placements")
    if not supported:
        return Decision("unresolved", "no_verified_geometry")
    placement, = supported
    state = "target_supported_local" if placement.role == target_role else "competing_supported_local"
    return Decision(state, "unique_witnessed_placement_in_completed_panel_scope", placement)
