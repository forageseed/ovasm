# Retained source and cross-organelle candidates

Seed recruitment preserves its original `recruit/numt_candidates.fastq` unchanged.
It also exports `recruit/source_candidates/`:

- `mitochondrion.fastq`: NUMT candidates compatible with mitochondrial bait.
- `plastid.fastq`: NUPT candidates compatible with plastid bait.
- `nuclear-insertion-candidates.jsonl`: read identifiers, length, best bait target,
  original classification scores and each target's compatible read intervals.
- `mt-pt-shared-reads.fastq`: accepted reads tied between MT and PT, before optional
  final target-depth subsampling.
- `summary.json`: counts, paths, classification parameters and interpretation scope.

A candidate may occur in multiple target files; do not sum these as unique reads.
Read intervals are 0-based half-open, gap-bridged bait-compatible blocks, not
reference alignments. Partial matching does not prove nuclear origin. Flanks
remain unassigned and `nuclear_origin_certified` is false. No nuclear chromosome
or insertion junction is inferred without independent nuclear-context evidence.
Short reads below `numt_min_block` cannot enter this pool; zero candidates do not
exclude NUMTs/NUPTs. These exports do not change recruitment or assembly choices.

After a both-organelle run, top-level `source_candidates/` contains:

- `mt-pt-shared-words.jsonl`: exact shared canonical words, **every** occurrence on
  both GFA graphs, canonical orientation, and per-target recruited-record support.
- `mt-pt-unique-seed-blocks.jsonl`: overlapping collinear exact words unique in
  each graph, merged into exact blocks with both segment coordinates and relative
  orientation. Ambiguous repeat occurrences remain in the word table.
- `mt-pt-sharing.json`: input paths, coordinate/read-support scopes and counts.

This detects exact sharing. It does not exhaust divergent homology, identify an
insertion junction or phase a repeat. A record containing a shared word supports
that word, not a particular genomic occurrence. MT and PT recruitment files may
overlap; their counts must not be summed. Missing recruited inputs yield null
support, while an existing empty file yields zero. Both graph files are required;
otherwise the cross-organelle report says `not_evaluated`.

`transfer_direction` is null and `physical_origin_certified` remains false. Exact
sharing can reflect transfer, repeats or contamination; it cannot distinguish
MT-to-PT from PT-to-MT. The assembly's graph, sequence and physical-topology claims
are not changed by this report.

To add sharing evidence to an existing both-organelle assembly without reassembly:

```sh
ovasm sharing --assembly existing-assembly --out evidence-report --k 31
```

The input graphs and recruited FASTQ are read only; reports go under the chosen
output's `source_candidates/` directory. New `run --organelle both` results include
the sharing report automatically.
