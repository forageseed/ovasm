# OVASM 0.2.1

Rust-only release: the Syncasm C source, its build script and FFI are removed. Graph
construction, correction, rescue and evidence are Rust on Linux and Windows.

## Assembly engine

- One engine (`assemble::run_engine`) with two coordinate encodings. Base coordinates are
  the default. `--rust-hpc` selects homopolymer-run coordinates (experimental): each read is
  compressed as it is read in (no compressed-read file), the same backbone, edge rule,
  cleaning and gap rescue build the graph, run lengths are restored from original reads, and
  the shared original-read junction stage follows. `--multichrom-graph` remains as a hidden
  legacy name and does not choose a chromosome count.
- `--rust-hpc` takes the same `--min-count` / `--min-count-fraction`; `--k` counts runs
  (default 1001). Give a read-derived `--min-count` for whole-genome input: the automatic
  value (node depth / 10, at least 3) suits recruited organelle reads and on whole-genome
  Schoenoplectus reads resolved to 3 and kept nuclear sequence.
- An edge between two kept nodes is also kept when its reads reach 0.35 of the shallower
  node's count (at least three reads). Short isolated pieces survive the first cleaning pass
  so read-agreed gap rescue can join them. Every unitig removed by cleaning is written to
  `<out>.cleaning.fa`.

## Repeat phases (opt-in; defaults unchanged)

- `ovasm unzip`: per repeat node, complete original reads per entry/exit phase, classified
  resolved / partial / conflicted / mixed / unresolved at the run's `min_count`. Never
  changes the input graph. `--out-gfa` adds supported `P` records; `--materialize-gfa`
  splits only fully resolved one-to-one repeats.
- `primary --phases` and `linearize --phases` forbid only zero-read pairings at repeats
  whose every end has a supported phase.
- `assemble --cleaning-audit-only` records cleaning proposals without deleting unitigs.

## Packages

- Precompiled x86_64 packages for Linux (glibc 2.34+) and Windows with checksum-verifying
  installers. Short bilingual README; full benchmark records are in `docs/`.

## Validation and limits

- Both platforms: 233 tests pass, 0 fail, 3 ignored, at the release source.
- Ten recruited inputs x two encodings: graphs of identical parameters are byte-identical
  to the previous release build; the opt-in phase stages left every default output
  unchanged.
- Schoenoplectus full input remains 0/8 record-for-record against the Oatk-derived
  reference (all eight records at least 99.97% covered, identity 1.0); that reference is
  not independent truth.
- Read correction is not implemented in the run encoding. Graph cycles and FASTA record
  counts do not certify physical chromosomes.
