# OVASM 1.1.0

## Plant organelle methylation profiling and epigenetic NUMT/NUPT separation

New pure-Rust base-modification engine (`ovasm check methylation` and `--methylation` /
`--methylation-bam` / `--max-read-methylation` in `run hifi`, `run ont` and `find recruit`), reading
standard PacBio HiFi and ONT Dorado modBAMs (`MM`/`ML` tags) directly with zero external C/htslib
dependencies:

- Calls per-site and per-segment **`CG`**, **`CHG`**, **`CHH`** 5mC and **`6mA`** methylation across
  assembled `molecules.fasta` and `organelle.gfa` (optional `organelle.methyl.gfa` tags `MC:f`,
  `MP:i`, `MD:f`, `MA:f`), verifying projected `CG`/`CHG`/`CHH` motif contexts directly against the
  assembled reference contig.
- Circular-wrap and symmetric repeat/IR-preserving colinear 31-mer placement retains **100% of
  genomic CpG sites** across plant plastid inverted repeats and large mitochondrial repeats where
  reference-linear `MAPQ` filters collapse.
- Adaptive background-calibrated read classification separates unmethylated true organelle reads from
  hypermethylated nuclear NUMT/NUPT insertions (`--max-read-methylation` filters hypermethylated reads
  during recruitment: on *Astragalus alpinus* HiFi modBAM, removing 31,610 hypermethylated NUMT reads
  restores the complete 313,827 bp circular mitogenome; on ONT R10.4.1 Dorado `sup,5mC_5hmC,6mA`
  modBAMs across six *Medicago* T2T species, `run ont --methylation --max-read-methylation` rescues
  *M. arabica* mtDNA from 247,231 bp to the complete 334,332 bp circular mitogenome with all 33 genes).

## Native ONT R10 Simplex phased overlap error correction

`run ont` and `assemble correct` now include a 4-stage `hifiasm (--ont)`-style phased overlap error
corrector (`--chemistry auto|r10|legacy`, `--r10`):

- Homopolymer-compressed (HPC) minimizer all-vs-all overlaps, unbridged chimera clipping, directional
  strand-bias CpG methylation and homopolymer site pruning, `O(n^2)` four-gamete dynamic-programming
  repeat phasing, and Phred-weighted base + run-length consensus, followed by `k = 501` graph probing.
- Under `--chemistry auto` (default), calibrated R10 Simplex FASTQ reads (`Q_eff` 15–25) select `r10`
  (validated across six *Medicago* ONT R10.4.1 SUP + `5mC`/`6mA` assemblies), while Dorado R9, PacBio
  CLR, PacBio HiFi/CCS reads, `run hybrid` and FASTA-derived constant qualities automatically keep the
  legacy path with byte-identical outputs across all eight regression datasets.

## Pan-genome, plastid layout and low-depth accuracy fixes

- Circular plastid molecules are oriented into canonical `LSC-IRb-SSC-IRa` order before downstream
  methylation coordinate projection, and plastids whose ends repeat around a boundary homopolymer close
  cleanly in `pan reconstruct`.
- `depth_peak` normalizes histogram bin widths by integer count support and takes the plateau median,
  eliminating the count-21 bin artifact on shallow (8–50×) targets and re-measuring close competing
  targets on non-shared bait k-mers.
- `pan reconstruct` reports library read length, fragment span and k-mer gate depth provenance, marks
  low-own-kmer segments uncallable in `genotypes.tsv`, and withholds configurations when an assembly
  is untrustworthy or belongs to the competing organelle.
- `pan cohort` applies read-length gating (`--min-read-length`), scores presence from segment-own
  k-mers (`>= 20`), reports `nearest.tsv` panel weights, and counts zero-molecule organelle runs
  cleanly; `KmerCounts::add` retries lost CAS updates under multi-threaded counting.
