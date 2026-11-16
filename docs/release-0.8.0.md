# OVASM 0.8.0

## Self-seeded recruitment across plants, animals and fungi

`run hifi` now extracts organelle protein-coding anchors directly from the sample's own whole-genome
reads across land plants, animals (`-x animal`, vertebrate and invertebrate mitochondrial genetic
codes) and fungi (`-x fungi`), with a reversed-sequence false-positive control (`--no-protein-seeds`
falls back to the built-in seed library only).

- Lineage-specific fallback seeds and shallow animal mitochondrial graph pruning keep divergent phyla
  from collapsing or retaining nuclear background.
- Leftover mitochondrial bootstrap is bounded by what its reads stand for at target depth, recovering
  missed multipartite circles and single-copy crossroads without running away into nuclear satellites.

## Audited seed provenance and cross-organelle candidates

- Normal runs use audited seed source views (`source_candidates/`), reporting exact shared words and
  collinear unique seed blocks between mitochondrial and plastid graphs (`mt-pt-shared-words.jsonl`,
  `mt-pt-unique-seed-blocks.jsonl`) without forcing unresolved repeats into a single transfer call.
- Graph cover search indexes bridge endpoints and deduplicates candidate covers, preserving depth
  uncertainty and partial-path evidence under bounded solver budgets.

## Everyday CLI options and `result.json`

- Short options follow PMAT2 where they overlap: `-i` (`--reads`), `-G mt|pt|all`, `-x 0|1|2`
  (`plant`, `animal`, `fungi`), `-g` (`--genome-size`), `-t` / `-T` (`--threads`), and `-o` (`--out`).
- Every `ovasm run` writes `result.json` (`ovasm.result.v1`; `--json` prints it to stdout) with run
  `status` (`complete`, `partial`, `no_molecule`, `failed`), molecule lengths and circularity, output
  paths, warnings, wall time and peak memory.
