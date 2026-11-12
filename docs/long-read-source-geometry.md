# Audit original long-read source geometry

`python3 scripts/audit_long_read_geometry.py PLAN.json --out NEW_REPORT.json` reproduces observed same-original junction geometry. It binds raw FASTQs, mapper queries and PAF, named reference records, registered windows, and a GFA with the original molecules used to build it. This invocation produces an audit report.

## Input

Use `ovasm.long-read-geometry.v1` as the plan's `schema` and provide a nonempty `scope`. Every file binding has `path` and lowercase `sha256`. Relative paths resolve from the plan directory.

| Plan field | Value |
| --- | --- |
| `original_fastqs` | One or more four-line FASTQ bindings. The first header field is the mapper query name. Each full header must also contain one native PacBio movie/ZMW identity, such as `m012345_123456/17/ccs`. |
| `mapper_queries` | FASTA binding containing exactly the automatically enrolled original query names and sequence bytes, normalized to uppercase. |
| `paf` | `data` binding for PAF with `cg:Z` and a `provenance` binding for mapper settings or receipt. PAF coordinates are zero-based and half-open. |
| `panel` | `fasta` binding for the mapper's combined reference and `records` mapping from each PAF reference name to `role`, `source` FASTA binding and original `source_record` name. Include mitochondrial, plastid and nuclear roles. Multi-record source FASTAs are supported. |
| `windows` | Registered windows as described below. |
| `window_provenance` | Binding for the source of the registered locations and physical boundaries. |
| `graph` | `gfa` binding, original FASTQ bindings under `reads`, and a `provenance` binding for their assembly relation. |

A window declares `id`, `reference`, the half-open `interval`, a positive `anchor_bp` width, `mode`, `exclusive_sides`, `require_full_query`, `graph_required` and `physical_boundaries`. Use an empty array for unknown physical boundaries. Set either both or one of `left` and `right` in `exclusive_sides` for source geometry according to the registered context.

- `exact_context` requires exact raw projection of the complete interval. It finds the nearest verified exterior anchors on either side through one CIGAR. Set `require_full_query` for whole-read context.
- `adjacent_boundary` uses the two adjacent exact guards. A gapped interior must have both registered endpoint pairs on the same CIGAR and aligned bases in between. A reference skip, split alignment or incorrect strand cannot provide the witness.
- `internal_control` records exact core projections. If two enclosing `physical_boundaries` are known, it reports the distance to the closest registered exterior, including whether the read is physically long enough to reach it.

## Interpretation

The audit enrolls reads from actual FASTQ sequence by exact whole interval or an ordered same-strand adjacent guard pair. It verifies all query identities and mapper query bytes before inspecting PAF. An exclusive anchor must occur once on its named reference record and only within its declared role in the combined linear panel. Both orientations and overlapping occurrences are searched in every record. Reference aliases or circular-origin placements must be represented explicitly by the input design.

When requested, graph corroboration requires one exhaustively spelled graph route, an exact sequence match in a different graph-input movie/ZMW, and no reported full competing role for the candidate's witness interval. The graph sequence is the registered core for `exact_context` and the core with both adjacent guards for `adjacent_boundary`.

Geometry-only records retain `unresolved_reported_paf_only` for competition. Alignment heuristics, secondary caps and partial reference coverage can hide competitors. The optional whole-read comparison establishes a separate result under the named reference and edit model. The audit never promotes a missing PAF placement to proven absence.

## Bind graph-frontier candidates

For graph-derived windows, add a `frontier_events` file binding to the plan.
Its schema is `ovasm.frontier-event-bindings.v1` with `discovery_report`,
`discovery_evidence`, and `records` fields. Resolve all file bindings from the
geometry plan directory. Use the complete output of `discover_graph_frontier.py`
on the same GFA and ordered original FASTQs. A mapper FASTA and PAF superset
may supply geometry, but only registered originals enter this audit. Their
other reported alignments and matching diagnostic windows can still oppose a
target claim.

Each registration contains:

| Fields | Value |
| --- | --- |
| `input_index0`, `record_index0` | Zero-based FASTQ input index and record ordinal. |
| `read_id_hex`, `sequence_sha256`, `quality_sha256`, `molecule` | Complete original header, raw sequence and quality digests, and movie/ZMW. |
| `event_index0`, `event_sha256` | Zero-based event index and digest of compact sorted JSON with the four original identity fields above, `event_index0`, and `event` from the discovery ledger. Omit the final newline. |
| `window_id`, `side` | Registered full-query `exact_context` graph window of length discovery k, with k-base exclusive anchors on both sides, and the left or right side of the raw event from which to derive its disjoint graph core. |
| `paf_row_sha256` | Digest of the complete literal PAF line including its line ending. Register each distinct line separately. |

The audit recomputes events for those originals and verifies the graph core and
adjacent discovery anchor on one oriented node. An independent graph-input
original must spell that entire footprint. One registered full-query CIGAR must
project both the core and the affected raw interval exactly. Only that CIGAR
can supply a source-competition bound. The donor supports the graph footprint.
The affected sequence is recorded separately and receives no graph-path claim
from that donor.

## Whole-read edit competition

Use `--full-edit-competition` to verify competing placements directly in the named
reference bytes. The bound for each original read is its literal number of
substitutions, insertions and deletions in a verified whole-read target CIGAR.
The same registered window must have role-exclusive exterior anchors on both
sides, a unique exact route through the graph interval, and a different
original movie/ZMW supporting that same sequence in the graph input.
For an `exact_context` window, the candidate also supports the full core
exactly. An `adjacent_boundary` window allows edits between the two junctions
when the two bases straddling each junction are exact and adjacent in the
same whole-read CIGAR. The gapped witness records the candidate's aligned
interval alongside the exact graph spelling.

The search partitions each oriented read into one more exact block than the edit
bound. Any continuous competing placement within the bound contains at least
one exact block. The program scans the competing references for all such blocks,
then computes full-query unit Levenshtein distance to contiguous substrings
in the resulting ranges. Equal or better placements withhold a candidate.
Verified geometry from another role also withholds the record.

Each record's `full_edit_competition` gives its target bound, the exact selected
CIGAR and witness, competing placements, and the decision. `--max-edit-seed-hits` and `--max-edit-text-bases` set work limits
per read and role. Exhausted searches stay unresolved. Use `--out-admission`
with this mode to export qualifying originals that were absent from the graph
input. The manifest binds exact original header, sequence and quality bytes to
the audit SHA-256. The existing `competition_complete` field continues to
describe reported PAF competition.

This comparison uses a linear named reference panel and unit edit costs.
Missing reference sequence, unrepresented arrangements, correlated errors and
source mixing are outside the search scope. Keep the model decision with its
reference bindings and graph provenance.

## Exclude a verified competing source

Use `--full-edit-competition --exclude-target mitochondrion --out-exclusions EXCLUSIONS.json`
to prepare a source exclusion for that recruitment destination. `plastid`
is also supported. The physical source comes from a registered window in
a different role. A single alignment must span the entire original read,
cover both role-exclusive exterior anchors and verify the registered
context or both exact junctions. The decision uses the literal edit count
of that alignment and searches every other supplied role at the same bound.
A verified witness in another role withholds exclusion, including a partial
witness at a registered window. A continuous equal or better placement or
an incomplete search also withholds exclusion.

A window marked `graph_required` also needs a unique route and an independent
original donor in that graph. An independently registered competing reference
window can omit that requirement. Exclusion applies only to the named
recruitment destination. The native manifest binds each exact original
to the audit, and an empty decision set produces only the audit. The result
describes the named source role and panel, with possible additional
placements within the same role.
