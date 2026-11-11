# Audit original long-read source geometry

`python3 scripts/audit_long_read_geometry.py PLAN.json --out NEW_REPORT.json` reproduces observed same-original junction geometry. It binds raw FASTQs, mapper queries and PAF, named reference records, registered windows, and a GFA with the original molecules used to build it. The only output is an audit report. The command has no native read-admission export.

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

When requested, graph corroboration requires one exhaustively spelled graph route, a different actual graph-input movie/ZMW whose read contains the exact core, and no reported full competing role for the candidate's witness interval. A graph observation is not evidence that a physical compartment is uniquely determined.

All records retain `unresolved_reported_paf_only` for competition. Alignment heuristics, secondary caps and partial reference coverage can hide competitors. A future qualification step needs independently justified evidence for the same original interval and explicit search completeness under its error and reference scope. The audit never promotes a missing PAF placement to proven absence.
