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

## ZM4 trial (2026-10-30)

Alfalfa cultivar ZM-4, existing reads only (nothing downloaded): ONT CRR2374622 (CRA034439,
R10.4.1 HAC, 131 Gb, reads of about 95 kb), HiFi CRR2374621 of the same project, Illumina
SRR13333747 (cultivar ZM-4; not shown to be the same plant). Prefixes of the files: ONT the
first 130,000 reads (12.3 Gb), HiFi the first 800,000 (15.8 Gb), Illumina the first 20 M pairs
(6 Gb). Built-in seeds, `--organelle both --threads 6`, CPUs 24-27,80-83, 64 GiB limit, ovasm
e92095b. Reference: the full-input ZM4 HiFi assembly (CRA032787: mitochondrion 160,953 +
146,073 bp, a 307,026 bp single cycle; plastid 125,509 bp); molecules aligned with minimap2
asm5, identity = matching bases / alignment columns.

| Run | Organelle | Molecules | Length (bp) | Circular | Reference covered | Identity | Wall | Peak RSS |
|---|---|---:|---|---|---:|---:|---:|---:|
| `run ont` 12 Gb | plastid | 1 | 125,456 | yes | 100% | 99.957% | 41.4 min | 4.09 GiB |
| | mitochondrion | 73 | 460,768 in all, longest 178,656 | 4 | 91.6% | 99.61% | | |
| `run hybrid` 12 Gb + 6 Gb | plastid | 1 | 125,498 | yes | 100% | 99.991% | 42.6 min | 5.26 GiB |
| | mitochondrion | 3 | 374,992 + 6,330 + 4,640 | 2 | 100% | 100.00% | | |
| `run hifi` 15.8 Gb (control) | plastid | 2 | 90,083 in all | no | 71.8% | 99.999% | 20.8 min | 1.49 GiB |
| | mitochondrion | 3 | 199,416 in all | no | 65.0% | 100.00% | | |

- Plastid, ONT alone: one circle of the reference's length (53 bp short), recruited 1,683 reads
  in round 1 (131 Mb, about 1,000x); what is left is the residual ONT error (about 1 in 2,300,
  homopolymer-like). Hybrid correction takes it to 1 in 11,000. Round 2 added under 5% and was
  not used. `check assess` calls it simple, without an inverted repeat, as the IRLC legumes are.
- Mitochondrion, ONT alone: no read was 80% covered by the library's seeds (round 1: 0
  reads); 142 partly covered reads made the round-2 bait, which recruited 177 reads (14 Mb,
  about 46x). Self-correction at that depth anchors 83% of the bases and leaves the graph in 73
  pieces that cover 92% of the reference. Depth, not recruitment, is what ONT alone lacks here:
  see the 38 Gb run.
- Mitochondrion, hybrid: the 177 corrected long reads recruited 2.6 M short reads (394 Mb); after
  hybrid correction the same reads assemble at k 1001 into a 374,992 bp circle that covers the
  whole reference without a mismatch in the aligned columns (its extra length is a further copy
  of repeats the reads do not resolve; decisive false), and the short reads support 26 of the 36
  graph links.
- HiFi at 15.8 Gb is not a complete control: its recruitment stopped at 199 kb of the
  mitochondrion and the plastid was downsampled to 931 reads (read-depth cap 250x) in two
  pieces. The full-input HiFi assembly is the reference above.
- Bugs these runs found and fixed: k 15 recruitment (every read), the round-1 refusal of an
  organelle no read was 80% covered for, and partly covered mitochondrial reads with
  plastid-derived stretches growing the plastid's bait (1.3 M k-mers for a 125 kb genome; round
  2 then took 920 Mb as plastid): partly covered reads now join a bait only at a quarter of its
  recruited reads' depth (909 of 976 left out).
