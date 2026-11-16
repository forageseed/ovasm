# OVASM 0.9.0

## Paired-end short-read recruitment and fragment filtering

- Short-read recruitment (`find recruit`, `run short`) reaches past bait edges using paired-end read
  mates (`MateReach`) that join in pieces closing open graph ends, stopping pair expansion if rapid
  nuclear growth appears rather than aborting the organelle target.
- Unplaced short-read graph fragments shorter than three read lengths without organelle seed support
  are withheld from `organelle.gfa` and `molecules.fasta` instead of being emitted as certified
  organelle chromosomes.

## Independent long-read junction and traversal validation

- Graph links and repeat traversals can be audited against independent long reads (`--dump-reads`
  lists the supporting or tied reads per link and traversal).
- Flanks lying inside repeat overlaps are skipped, and reads tied across alternative repeat probes
  leave a traversal ambiguous rather than falsely contradicted.

## Homopolymer polishing and multipartite linearization

- `run ont` polishes assembled homopolymer runs against the uncorrected reads' own length-dependent
  error model, and `run hybrid` judges homopolymer run lengths directly from the paired short reads
  whenever short-read coverage is present.
- `solve linearize` counts single-copy bases outside repeat overlaps, treats a mostly single-copy
  unit as a subcircle rather than a tandem copy, and closes linear walks when reads cross a link a
  genome apart.
