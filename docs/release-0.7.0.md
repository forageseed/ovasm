# OVASM 0.7.0

## Faster recruitment

Recruiting organelle reads from a whole HiFi library is the slowest part of `run hifi`. It is faster
now, and the output is byte-identical to 0.6.0 (38 output files compared on Potentilla sterilis and
Eleocharis quinqueflora):

- read lookups go through a Bloom filter and an exact table with prefetching; reads whose k-mers cannot
  change the result are skipped; the depth sketch is built during the first pass; encoding and
  classification run in parallel;
- replay passes scan the cached reads in their 2-bit form;
- `OVASM_GPU=1` checks the changed k-mers of replay passes on an NVIDIA GPU (prototype). The kernel
  is shipped as PTX and loaded through the driver; without a GPU, or if loading fails, OVASM falls back
  to the CPU and the output is the same.

| Same server, 8 threads | Wall time | Peak memory |
|---|---:|---:|
| Potentilla sterilis, 0.6.0 | 15:48 | – |
| Potentilla sterilis, 0.7.0 CPU | 12:10 | 2.2 GB |
| Potentilla sterilis, 0.7.0 `OVASM_GPU=1` | 9:50 | 2.2 GB |
| Eleocharis quinqueflora, 0.6.0 | 51:40 | – |
| Eleocharis quinqueflora, 0.7.0 `OVASM_GPU=1` | 36:07 | 2.6 GB |

`OVASM_PROFILE=1` prints the time spent in each recruitment stage.

## Depth gate for repeat-rich and polyploid genomes

`run hifi` only lets reads extend the mitochondrial bait when they look deeper than the nucleus. The
nuclear depth was taken from the peak of per-read median k-mer depth. In genomes full of repeats, or
polyploids, that peak sits at twice the real nuclear depth, so the gate shut out the mitochondrion
too. In *Lespedeza davurica* (allotetraploid) the read-median peak was 56× against a k-mer peak of
28×; the mitochondrion, at about 72×, stopped growing at 250 kb of 379 kb.

- Extension now follows a target that is only a little deeper than the nucleus.
- When the read-median peak is 1.7 to 2.5 times the k-mer peak (and the k-mer peak is at least 10×),
  extension measures against the k-mer peak and a warning says so.
- New options: `--genome-size` (for example `1.26g`) or `--nuclear-depth` to give the nuclear depth
  yourself; `--gate-fold`; for `run hifi`, `--extend-gate-fold` and `--no-extend-gate`.

*Lespedeza davurica* now assembles its complete 378,667 bp mitochondrion with default settings
(34:43, 2.4 GB). Carex and Lophophora give the same recruited reads as before.

## ONT reads

New `run ont` assembles from ONT reads alone and `run hybrid` from ONT plus Illumina reads of the
same sample ([how it works](long-reads.md)). On 38 Gb of alfalfa ZM4 ONT reads, `run ont` gave the
plastid and both mitochondrial circles of the reference (52 min, 3.4 GiB). PacBio CLR reads also go
through `run ont`, but it has not yet been checked on real CLR data.
