# Noisy long reads: `run ont` and `run hybrid`

`ovasm run ont` assembles organelles from ONT or PacBio CLR reads alone; `ovasm run hybrid`
from such long reads together with short reads (Illumina) of the same sample. Both end in the
same stages and files as `run hifi`.

```sh
ovasm run ont    --reads sample.ont.fastq.gz --organelle both --threads 8 --out ont
ovasm run hybrid --long sample.ont.fastq.gz \
                 --short sample_1.fastq.gz --short sample_2.fastq.gz \
                 --organelle both --threads 8 --out hybrid
```

Seeds are the built-in land-plant library unless `--seeds ORGANELLE=PATH` is given, as for
`run hifi`. `--read-set target` takes long reads that are the organelle's already (one
organelle) and skips their recruitment; `run hybrid` still recruits the short reads.

## Why the flow is not `run hifi` on noisy reads

- **No extension from raw reads.** An ONT read's error k-mers are short random words that the
  nuclear genome holds too, so a bait extended from raw reads runs into the nucleus. `find
  recruit --preset ont` makes one pass and never extends; the bait grows only from corrected
  reads.
- **k 21, not 15.** At the ONT preset's k 15 the built-in library (tens of Mb of k-mers, about a
  twentieth of all 15-mers) covers every read by chance: ZM4's ONT library (CRR2374622) gave
  61 of its 131 Gb as "plastid" in an earlier batch run. At k 21 a chance hit is about a hundred
  thousand times rarer, and a read with 5% errors keeps a third of its 21-mers.
- **Partly covered reads grow the bait.** A read joins a target when the bait covers 80% of it.
  The reads that reach past the bait's edge are therefore only a fifth of the depth there, below
  the solid count of correction and of `--seed-reads`; with a multi-species library a whole
  ONT read of a plant mitochondrion (genes are a few percent of it) is rarely 80% seed. The
  reads holding a 2 kb block of the bait (`numt_candidates.fastq` of the pass) keep the
  organelle's depth along their whole length, so they are corrected with the recruited reads to
  make the next round's bait. NUMT reads are among them, but at nuclear depth, so their flanks
  stay out of the solid k-mers; and they are never assembled.

## Stages

1. `find recruit --preset ont --k 21 --iterations 1` with the seeds (`recruit/round1/`).
2. Rounds 2 to `--recruit-rounds` (default 2): the previous round's recruited and partly covered
   reads, self-corrected (`correct/round<N>/<organelle>.bait.fasta`), join the seeds as
   `--seed-reads`. A round that adds less than 5% of the recruited bases is not used and ends the
   rounds. Each round reads the whole input again.
3. The last used round's recruited reads, self-corrected on their own (`assemble correct`, k 17
   then 25): `correct/<organelle>.fasta`.
4. `run hybrid` only: those corrected reads seed one short-read pass (`find recruit --preset sr
   --seed-reads --mates-together`, `recruit/short/`), and the raw long reads are corrected again
   against the recruited short reads, then self rounds for what the short reads do not cover
   (`correct/hybrid/<organelle>.fasta`). Short reads cannot come first: seeds alone recruit the
   short reads of the conserved genes only.
5. The HiFi stages on the corrected reads (`run hifi --read-set target`): sparse graph at the
   automatic k with the k ladder, `check evidence`, `solve linearize`, `assemble unify`, `solve
   chromosomes`, `check assess`.
6. `run hybrid` only: `check evidence --short-read-adjacency` with the short reads on the same
   graph: `evidence.short.json`.

`--recruit-max-bases` (400 Mb) and `--short-max-bases` (1 Gb) bound what correction holds in
memory.

## summary.json

That of `run hifi`, with `read_type` `ont` or `ont_hybrid`, `recruited_reads/bases` of the last
used round, and:

- `long_read.rounds`: per round the bait (seeds, or seeds + corrected reads), recruited and kept
  reads and bases, whether the round was used, and its seconds;
- `long_read.corrections`: each correction (`purpose` bait or assembly, hybrid or not, rounds,
  reads out);
- `long_read.short_reads` (hybrid): recruited and kept short reads;
- `short_read_evidence` (hybrid): links and short-read-supported links of the graph;
- `total_elapsed_seconds` and `peak_rss_mb` of the whole run.

## ZM4 trial (2026-10-30/31)

Alfalfa cultivar ZM-4, existing reads only (nothing downloaded): ONT CRR2374622 (CRA034439,
R10.4.1 HAC, 131 Gb, reads of about 95 kb), HiFi CRR2374621 of the same project, Illumina
SRR13333747 (cultivar ZM-4; not shown to be the same plant). Prefixes of the files: ONT the
first 130,000 reads (12.3 Gb) or 400,000 (38 Gb), HiFi the first 800,000 (15.8 Gb), Illumina
the first 20 M pairs (6 Gb). Built-in seeds, `--organelle both --threads 6`, CPUs 24-27,80-83,
64 GiB limit, ovasm e5c2a3c (e92095b in the second table). Reference: the full-input ZM4 HiFi
assembly (CRA032787: mitochondrion 160,953 + 146,073 bp, a 307,026 bp single cycle; plastid
125,509 bp); molecules aligned with minimap2 asm5, identity = matching bases / alignment
columns.

| Run | Organelle | Molecules | Length (bp) | Circular | Reference covered | Identity | Wall | Peak RSS |
|---|---|---:|---|---|---:|---:|---:|---:|
| `run ont` 12 Gb | plastid | 1 | 125,438 | yes | 100% | 99.943% | 13.5 min | 3.04 GiB |
| | mitochondrion | 6 | 307,562 in all, longest 207,572 | 0 | 99.95% | 99.987% | | |
| `run hybrid` 12 Gb + 6 Gb | plastid | 1 | 125,498 | yes | 100% | 99.991% | 29.4 min | 3.74 GiB |
| | mitochondrion | 2 | 298,803 + 1,215 | 0 | 97.3% | 100.00% | | |
| `run ont` 38 Gb | plastid | 1 | 125,438 | yes | 100% | 99.943% | 52.1 min | 3.44 GiB |
| | mitochondrion | 2 | 160,947 + 146,059 | 2 | 100% | 99.993% | | |
| `run hifi` 15.8 Gb (control) | plastid | 2 | 90,083 in all | no | 71.8% | 99.999% | 20.8 min | 1.49 GiB |
| | mitochondrion | 3 | 199,416 in all | no | 65.0% | 100.00% | | |

ovasm e92095b on the same inputs:

| Run | Organelle | Molecules | Length (bp) | Reference covered | Identity | Wall | Peak RSS |
|---|---|---:|---|---:|---:|---:|---:|
| `run ont` 12 Gb | plastid | 1 | 125,456 circle | 100% | 99.957% | 41.4 min | 4.09 GiB |
| | mitochondrion | 73 | 460,768 in all | 91.6% | 99.61% | | |
| `run hybrid` 12 Gb + 6 Gb | plastid | 1 | 125,498 circle | 100% | 99.991% | 42.6 min | 5.26 GiB |
| | mitochondrion | 3 | 374,992 circle + 6,330 + 4,640 | 100% | 100.00% | | |
| `run ont` 38 Gb | plastid | 1 | 125,432 circle | 100% | 99.938% | 316.0 min | 11.18 GiB |
| | mitochondrion | 1 | 307,000 | 100% | 99.991% | | |

- Plastid, ONT alone: one circle 71 bp shorter than the reference, recruited 1,683 reads in
  round 1 (131 Mb, about 1,000x); what is left is the residual ONT error (about 1 in 1,800,
  homopolymer-like). Hybrid correction takes it to 1 in 11,000. Round 2 added under 5% and was
  not used. `check assess` calls it simple, without an inverted repeat, as the IRLC legumes are.
- Mitochondrion, ONT alone: no read was 80% covered by the library's seeds (round 1: 0
  reads); the partly covered reads at the mitochondrion's own depth made the round-2 bait, which
  recruited 79 reads (6.6 Mb, about 21x). At 12 Gb they assemble into 6 pieces covering 99.95%
  of the reference; at 38 Gb 266 reads (20.8 Mb, about 68x) give the reference's two circles at
  k 1001 (12 of 12 links supported; decisive false, as the two are also one 307 kb cycle).
- Mitochondrion, hybrid: the 79 corrected long reads recruit the short reads that correct them;
  the graph at k 101 has a 298,803 bp piece, and the short reads support 78 of its 100 links.
  8 kb of the reference (2.7%) is missing; e92095b, from 177 noisier reads, had covered it all
  in a longer circle.
- HiFi at 15.8 Gb is not a complete control: its recruitment stopped at 199 kb of the
  mitochondrion and the plastid was downsampled to 931 reads (read-depth cap 250x) in two
  pieces. The full-input HiFi assembly is the reference above.

What changed from e92095b to e5c2a3c, and why:

- **Partly covered reads are gated at their own depth when nothing was recruited.** The
  mitochondrion's round 1 recruits nothing, so the e92095b gate (a quarter of the recruited
  reads' depth) kept every partly covered read, nuclear ones too. Without recruited reads they
  are now compared with one another and kept at a quarter of the 90th percentile of their
  depths. At 12 Gb that keeps 64 of 72 mitochondrial reads and 14 of 60 others; at 38 Gb 224 of
  255 and 35 of 190.
- **Correction depth is capped at 300x** (`depth_cap` of each correction). The plastid's reads
  are about 3,000x at 38 Gb, where systematic ONT errors pass the solid count; above 300x a
  hash-chosen share of the reads is corrected. With the cleaner bait this took the 38 Gb run
  from 5.3 h and 11.2 GiB to 52 min and 3.4 GiB.
- **An organelle whose round was sampled at `--recruit-max-bases` takes no further round**
  (`capped_before` in its next round's entry); the rounds end when no organelle below the cap
  grew 5%. In ZM4 the plastid's 38 Gb round 1 stopped at 399,979,844 bases, just under the cap,
  so it still ran round 2.
- **Short tandem repeats are masked in baits** (`tandem_masked_bases`): runs of period 1 to 6
  and at least 21 bases are N in the next round's bait and in the long-read seeds of the
  short-read pass, not in the reads that are assembled. Through (GA)n k-mers of the bait the
  mitochondrion had recruited two nuclear microsatellite reads; hybrid correction lengthened
  their runs and the automatic k took that k-mer as the node depth, leaving no segment.
- Earlier fixes from these runs: k 15 recruitment (every read), the round-1 refusal of an
  organelle no read was 80% covered for, and partly covered mitochondrial reads with
  plastid-derived stretches growing the plastid's bait (1.3 M k-mers for a 125 kb genome).
