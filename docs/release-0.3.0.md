# OVASM 0.3.0

HiFi recruitment in `ovasm run` is gated by read depth by default and no longer runs into the
nuclear genome on species far from the built-in seed library.

## Recruitment (`ovasm run`, HiFi reads)

- Gated by default: `run` recruits HiFi reads with `ovasm recruit --extend-gate --bootstrap`;
  `--ungated-recruit` restores 0.2.1. Ungated, iterative baiting from the built-in seed library
  ran into the nucleus on Carex echinata (686,022 "plastid" reads, 8.5 Gb) and assembly took
  301 GB before the run failed after 8 hours.
- The depth gate reads the nuclear depth from per-read depths: the mode of the reads' median
  sketch depths in quarter-log2 bins, accepted when its reads stand for a nuclear-sized genome
  (20 Mb), else the k-mer estimate. The mode of distinct k-mers, used until now, fell to 2x on
  three of four benchmark libraries, where a few percent of shallow reads hold more distinct
  k-mers at count 2 than any count of the nuclear peak (organelle threshold 6x; about 62% of
  Carex reads passed). `recruit.json` reports both estimates.
- Reads are cached for recruitment in `<out>/recruit/input-cache` when that disk has room, and
  removed when done (`--recruit-cache`, `--no-recruit-cache`): passes after the first and the
  depth sketch replay the cache instead of decompressing the input again.
- Each organelle's recruited reads are capped at about 250x read depth (`--recruit-max-depth`)
  and 1 Gb (`--recruit-max-bases`) before assembly.
- Up to 60 passes (`--recruit-iterations`). A read joins a target only when the bait covers 80%
  of it, so a pass extends the bait by about a tenth of a read length past each edge, and a
  large mitochondrion grows by about as many reads every pass; with the gate such a target is
  not taken for saturated.
- A high-copy repeat recruited with a target sets neither the target's depth nor its extension
  (the Lophophora williamsii plastid: most first-pass reads were a repeat at ~700,000x).
- Extension of very deep targets keeps recurrent errors out, and subsamples its reads to 200x by
  the target's own read depth.

## Linearization

- An automatic cover envelope (`linearize` without `--max-molecules`, as `run` calls it) that
  holds no cover widens one molecule at a time until one does. A branched component needing
  more molecules than the envelope went entirely to the partial paths; components with a cover
  inside the envelope are searched as before.

## Outputs

- Same layout as 0.2.1: `recruit/<organelle>.fastq`, `<organelle>/molecules.fasta`,
  `<organelle>/assembly.gfa`, `summary.json`. `recruit/<organelle>.fastq` holds the capped reads
  (`recruit.json`: `recruited_reads`, `kept_reads`, `sampling_fraction`, `sampled_by`).
- `recruit.json` adds `depth_gate`, per-pass read depths and implied sizes, the reads left out
  of extension (`extension_gated_out`, `extension_shallow_out`, `extension_deep_out`) and each
  target's `stop_reason`.
- The built-in seed library is unchanged (v20261007).

## Validation

Four benchmark HiFi libraries, `ovasm run --organelle both`, 16 threads, 2026-10-22. 0.2.1
gated: the benchmark's `recruit --extend-gate --bootstrap` and `run --read-set target` per
organelle.

| Library | Version | Wall | Peak memory | Mitochondrion | Plastid |
|---|---|---|---|---|---|
| Carex echinata | 0.3.0 | 67 min | 1.5 GB | 1 molecule, 1,296 kb (whole-library PMAT2: 1,332 kb) | 179 kb circle |
| | 0.2.1 gated | 37 min | 1.6 GB | 6 molecules, 1,162 kb | 179 kb circle |
| | 0.2.1 default | failed, 8.0 h | 301 GB | – | – |
| Lophophora williamsii | 0.3.0 | 72 min | 3.7 GB | 7 molecules, 1,181 kb (graph 1,175 kb) | 116 kb circle |
| | 0.2.1 gated | 55 min | 3.2 GB | 7 molecules, 535 kb | 116 kb circle |
| | 0.2.1 default | timeout, 12 h | 20 GB | 7 molecules, 586 kb | – |
| Medicago polymorpha | 0.3.0 | 14 min | 1.3 GB | 1 molecule, 227 kb (graph 267 kb) | 124 kb circle |
| | 0.2.1 gated | 67 min | 3.2 GB | 236 molecules, 8.1 Mb (29% of its reads over 10x the mitochondrion's depth) | 124 kb circle |
| Oryza sativa Nipponbare | 0.3.0 | 32 min | 1.5 GB | graph 418 kb; gold 16/16 junctions, 8/8 branches | 135 kb |
| | 0.2.1 gated | 2.1 h | 1.4 GB | graph 418 kb; gold 16/16, 8/8 | 135 kb |

- Medicago's mitochondrion from the same recruited reads by the benchmark's other assemblers:
  Oatk 234 kb (1 contig), PMAT2 238 kb, HiMT 362 kb.
- Linux: 245 tests pass, 0 fail, 4 ignored, at the release source.

## Limits

- Large mitochondria take many passes. The Lophophora williamsii mitochondrion (2.19 Mb by
  whole-library HiMT) was still growing after 60 passes, its reads standing for 1.11 Mb;
  `--recruit-iterations` raises the limit.
- Windows packages are not yet built or tested for 0.3.0.
