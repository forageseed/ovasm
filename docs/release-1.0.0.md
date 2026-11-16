# OVASM 1.0.0

## Species pan-graph reconstruction (`pan reconstruct` v4) and `pan cohort`

`ovasm pan reconstruct` and the new `ovasm pan cohort` bring individual organelle genome
reconstruction and population-scale pan-genome analysis to thousand-sample cohorts:

- **Two-layer individual reconstruction**: writes both `<organelle>.reads.fasta` (joins cut across
  repeats that the target's own reads could span yet do not support, recorded in `junctions.tsv`)
  and `<organelle>.fasta` (representative configuration on the species panel, ranking all candidate
  configurations and spelling the top 20).
- **One-shot panel recruitment**: the species panel recruits in a single pass by default, with target
  reads extending only past open panel edges (`--iterative-recruitment` restores multi-pass baiting);
  fragmented components receive partial matchings instead of dropping out.
- **Cohort summary (`ovasm pan cohort`)**: aggregates hundreds to thousands of `pan reconstruct` runs
  on a shared panel into segment presence/absence, read-supported junction tables, variant alleles
  and closest panel-member weights.

## Faster, safer short-read recruitment

- Paired-end `MateReach` probing runs in parallel alongside read classification, and 2-bit replay
  passes lazily decode only the reads re-examined in later iterations (cutting 6-pass ZM4 short-read
  wall time from 34:10 to 28:43 with byte-identical outputs).
- Targets placed at or under `max_read_depth` by the k-mer gate (`kmer_read_depth`) are never
  downsampled when plastid-derived shared reads inflate the per-read depth sketch (`cap_depth_from =
  "kmer_under_cap"`).
