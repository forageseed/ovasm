# OVASM 0.2.1

- Rust graph implementations on Linux and Windows, including the experimental
  homopolymer-coordinate multichrom graph.
- Eight-chromosome synthetic recovery and original-read consensus accounting.
- Precompiled x86_64 packages and checksum-verifying installers.
- Short bilingual README; full benchmark records remain in `docs/`.

- One assembly engine (`assemble::run_engine`): homopolymer-run coordinates are an
  internal encoding of the same backbone, edge rule, audited cleaning, gap rescue and
  original-read junction stage, not a separate multichromosome mode.
- An edge between two kept nodes is also kept when its reads reach 0.35 of the shallower
  node's count (at least three reads). The default-route graphs of ten paired regression
  inputs were byte-identical with and without this rule.
- Short isolated pieces survive the first cleaning pass so read-agreed gap rescue can join
  them; every unitig removed by cleaning is written to `<out>.cleaning.fa`.

- Run compression now happens as reads are read in, inside the same backbone loop: the
  11.6 GB `reads.hpc.fasta` intermediate (Schoenoplectus full input) is no longer written.
  `--rust-hpc` takes the same `--min-count` / `--min-count-fraction` as base coordinates and
  `--k` defaults to 1001 runs. On ten paired inputs the graphs of identical parameters are
  byte-identical; an omitted `--min-count` is resolved from the run-coordinate node depth
  (different from the base-coordinate value) and has not been validated on whole-genome input.

HPC mode does not yet implement read correction or general repeat unzipping. Graph cycles
and FASTA record counts do not certify physical chromosomes.
