# Eight-run HiFi panel and short-read target validation

The graph-conditioned reconstruction command was validated on an existing alfalfa
eight-run HiFi panel and the usable paired prefix of ZM4 WGS SRR13333747. This is
development validation on known data, not a new blind test or independent genome truth.
The eight source runs are SRR31693270, SRR31693266, SRR31693271, SRR31693262,
SRR31693276, SRR31693269, SRR31693268 and SRR31693272. Run labels do not certify
eight different biological individuals. Neither ZM4 nor CRR2374621 was included in
this reference panel. All source mitochondrial graphs retain unresolved configurations.

The original WGS files have official MD5 mismatches. Their use was explicitly
authorized, but a subsequent FASTQ parser found a malformed mate-2 record at pair
77,777,447. The retained input contains exactly the preceding 77,777,446 complete
matching pairs: 155,554,892 records and 23,333,233,800 bases. No read sequence was
cropped, and the malformed record and everything beyond it were not claimed as valid
input. This prefix is not the complete run. Cultivar labels do not establish a
same-individual pairing with an independent HiFi truth.

## Actual staged results

All measurements and sequencing computations ran on forageseed-235 with eight threads.
The first all-panel raw-input recruitment/projection stage took 1,270.446 seconds and
peaked at 0.763 GiB in sampled process-group RSS (GNU maximum RSS 799,092 KiB).
It kept 350,083 mitochondrial and 128,390 plastid records. Both logical recruitment
rounds scanned all 155,554,892 records and 23,333,233,800 bases. Its de novo stage
failed because the new wrapper passed 10 replicates to the established 1..5 boundary.
The first run and its reports were retained, including its then-incorrect exit-zero
workflow interpretation. The corrected command now preserves the error report and
exits with an error if a de novo stage fails.

The corrected final target stage reused those exact retained files and all eight
reference paths; it did not rescan the full WGS prefix. It took 47.154 seconds, with
0.630 GiB sampled process-group RSS and GNU maximum RSS 660,376 KiB. These two
stages constitute a measured staged validation, not a whole-input rerun of the final
binary. Original single-reference and failed attempts were retained separately.

| Target | Certain exact reference-node bases | Individual GFA stored bases | Connected components | Largest component bases |
| --- | ---: | ---: | ---: | ---: |
| Mitochondrion | 291,760 | 301,574 | 43 | 290,784 |
| Plastid | 124,999 | 128,970 | 1 | 128,970 |

Stored node bases include alternative branches; they are not a physical genome
length. The single-reference plastid de novo graph had 14 disconnected components.
The eight-panel individual graph has one connected component, but long-repeat phase
remains unresolved. Its nonbranching FASTA has 663 records; the mitochondrial FASTA
has 1,078. These are graph contigs, not chromosomes. The mitochondrial de novo main
component remains unsolved and has 15 partial paths; the individual graph retains its
sequence, unbridged evidence and small separate components. Neither target is reported
as a complete physical assembly.

The real data exposed two general graph/path defects in addition to the replicate
boundary: duplicate one-base allele representations could not join through a node-
internal 63-mer, and ten path ends required at least five covering paths while the
automatic search stopped at four. Exact sequence plus both oriented neighbour sets
now coalesces identical parallel representations (400 MT and 164 PT nodes in this
case); all input sequences and paths remain verified. The automatic search now includes
the necessary path-end bound, with explicit positive limits and expansion budgets
preserved. On this dataset, searching size five still found no full MT cover; the
report retains that failure and does not promote the small 1,455-bp separate record
to a complete mitochondrial result.

The final read-only audit verified every individual path edge, exact FASTA/path
spelling, every input segment after coalescence, and complete reference-interval
partitions. Unsupported and ambiguous-only reference intervals were not imported.
The raw target graphs, all biological evidence fields and individual FASTA/GFA were
unchanged by the cover-bound fix; only the recorded search range changed.

## Checks and remaining limits

The server scoped Rust suite passed 184 tests, with zero failures, three existing
ignored tests and three correction tests excluded from this scoped run. Windows
passed the nine new graph-reconstruction tests and the new five-path cover-bound test.
The OrganelleVerse Python entry point `ov.asm.reconstruct` completed the small real
native two-organelle, second-source-allele fixture in raw-input and explicit-target
modes: eight API tests passed on Windows and eight on the server. Existing evidence,
pipeline and facade checks contributed another eight passing tests on each host.
An early Windows check lacking `httpx` failed two desktop-boundary tests; the same
checks passed with the existing backend environment. First failures are retained in
the task evidence. No new default whole-project suite was claimed as passing; the
earlier unrelated default-suite file-descriptor exhaustion remains a recorded limit.

Exact matching is a local compatibility test, not a source or transfer-direction
certificate. Sequencing errors and alleles absent from all source paths use the
separate de novo sequence channel. A new recombinant configuration absent from source
P paths is not inferred by copying a reference path. Mate insert-size phasing is not
implemented by this command. Strict and relaxed scores remain unevaluated without
independent compatible truth, physical molecule count remains null, and physical
topology remains unidentified. No desktop/browser reconstruction workflow was certified
by these native and Python checks.

Terminal run IDs are `ovasm-eight-pan-WGS-v1`, `ovasm-eight-pan-WGS-target-v2` and
`ovasm-eight-pan-WGS-target-v3`. Final native source was
`d84b60f78811a97e8519134acbbf097ac0bd521f`. Reports, manifests, resource samples and
the independently inspected terminal outputs are retained under the server project's
`pan-wgs-validation` task directory.

## Exact short-read adjacency and repeat classification repair

The retained-read follow-up `ovasm-mt-pan-wgs-adjacency-repeat-v1` used source
788d4d8 and the same target files and eight reference runs. Two dense-overlap
edges had 35 and 23 unique exact adjacency records but zero strict long-span
crossings. Repeat classification previously ignored those adjacency records,
misclassifying a branching 360-bp, depth-42.3 segment as a single-copy anchor.
Unique adjacency now contributes to identifying a branch; it never contributes
to read traversals or to selecting a repeat entry/exit pairing.

The main MT component now has 28 complete anchor-cover candidates instead of
zero. Its selected five linear representatives have lengths 94,370, 72,825,
72,733, 44,881 and 6,022 bp; the independent 1,455-bp component remains. The
combined representative length is 292,286 bp. All input graph sequence and
reference projections are unchanged, and the plastid representative FASTA is
byte-identical. MT `decisive=false` remains: the selected cover has ten
topological bridges and the alternative scores are nearly tied. This repair
resolves the anchor-cover failure; it does not identify physical chromosomes.

The independent terminal audit verified input preservation and every final
path/FASTA spelling. Scoped server checks passed 24 evidence tests, 36 path
tests (one existing ignored), nine reconstruction tests and eight public API
tests. Windows passed 24 evidence tests and eight public API tests. The first
server build lacked the existing zlib include path; it was preserved and the
build succeeded with the existing py313 headers/libraries. No default full
Python suite was rerun or reclassified as passing.

A held-out candidate comparison finds 15,212 bp absent from these representatives.
This is a development sequence comparison, not independent truth. A separate
full-prefix raw-read audit scanned all 155,554,892 records in 565.182 seconds and
found exact compatibility across all those missing bases. The two large intervals
(8,490 and 6,645 bp) had 34,024 and 18,147 compatible query names, respectively.
Those names can denote fragments where mates share a name; counts are not source
certificates. The eight-source panel matches only 561 bp of the first interval
and none of the second. The original seed recruitment stopped after two of its
20 allowed rounds, with net growth below the existing five-percent rule; this
was not a fixed-point test. No reads were rescued using held-out coordinates.
The separate reference-free WGS novelty-channel pilot is still under evaluation;
complete reconstruction has not been claimed.

## Recruitment to a fixed point, and the frontier walk

Recruitment was first rerun at f6d1745 with `--saturation 0 --iterations 60` and the
same bait (as `find recruit` with `--max-bases` 150 x 350,586). It never settled: the
mitochondrial reads went from 856,248 (pass 2) to 940,026 (pass 4) and 4,630,842
(pass 5), then to 28.8 M at pass 60, and the plastid's from 7.35 M to 48.4 M. The
150x subsample was then mostly nuclear; the reconstructed mitochondrion fell apart
(largest component 9,753 bp, 230 components, de novo molecules at most 1.8 kb) and
neither interval was closed (4:37:36, 7.4 GB peak).

A one-pass diagnostic recruitment on the held-out candidate's 20.6 kb around the two
intervals (diagnostic only, never used by `pan reconstruct`) found them in the library at
the depth of their flanks: 399 and 390x over all primary alignments, 368 and 357x for
reads within 1% of the candidate along 90% of their length, against 369-370x for the
2 kb on either side, with no 250 bp window below 293x and no base without such reads.

The frontier walk (step 0 of the reference reconstruction workflow) was then run on the
same inputs and command (`--target-depth 150`). Counting 155,554,892 records (23.3 Gb)
took 646 s. The mitochondrion's depth came out at 271 (31-mer counts), its walk began at
2,204 places, 86 found new k-mers, 2,109 steps rejoined the reference, and it added 15,501
31-mers in 93 pieces (18,291 bp); the plastid (depth 6,932) added 66. Recruitment then
took 891,167 mitochondrial reads (856,256 before) and stopped at the same pass.

| | f6d1745 | frontier walk |
|---|---|---|
| Wall clock / peak RSS | 7:41 / 1.18 GB | 18:25 / 4.73 GB |
| MT stored bases / components / largest | 301,574 / 43 / 290,784 | 314,620 / 41 / 305,875 |
| MT de novo molecules (bp) | 94,370, 72,825, 72,733, 44,881, 6,022, 1,455 | 238,651, 66,913, 454 |
| Interval 57245+8490 in the individual graph / molecules | 590 / 0 bp | 6,890 / 7,010 bp |
| Interval 67190+6645 in the individual graph / molecules | 0 / 0 bp | 6,645 / 6,645 bp |
| MT molecule bases matching neither the candidate nor the plastid | 95 | 95 |
| Plastid graph and molecules | | identical |

Interval coverage is exact matching to the held-out candidate (minimap2 `=` runs), a
development comparison, not independent truth. Of the walked bait, 5,655 bp matched
neither the candidate nor the plastid, in pieces under 500 bp; the individual graph's
sequence matching neither stayed at 33,364 bp (33,365 before). 1,480 bp of the first
interval are still missing from the molecules. `decisive=false` remains.

### Walk past short shared stretches; several k; a final record

The first frontier walk stopped 1,480 bp short inside the first interval: about 30 bp
there are held elsewhere in the reference, and the walk took them for a rejoin from both
sides. Those reads were then at the depth of their flank (381x against 386x for the 1.5 kb
before them, one recruit pass each), not a nuclear repeat. A walk now counts as rejoined
after 300 reference k-mers in a row. Letting it walk that far with the first thresholds
(a quarter of the depth, comparable branches at a quarter) took it into the nucleus at
NUMT edges (67,280 new k-mers in 4,221 pieces; 2.34 M mitochondrial reads recruited); with
0.4 of the depth, branches at half and dead-end pieces under 150 bp dropped, it added
17,177 k-mers in 91 pieces and recruitment kept 895,480 reads (saturated at pass 2).

The short-read de novo assembly at k 127 left the plastid in four pieces at up to 1,000x
of reads and in twelve at 150x (its 127-mers ~21x deep). At k 99, 77 and 55 it was one
molecule of 124.9-125.1 kb (a 136 bp piece beside it at 99), so short reads now try those
k as well (see the reconstruction workflow).

| `--target-depth` | 1000 | 150 |
|---|---|---|
| Wall clock / peak RSS | 18:08 / 4.73 GB | 20:03 / 4.73 GB |
| MT de novo molecules (bp) | 307,479 (k 127) | 307,113, 146 (k 99) |
| Both intervals in the MT molecules | 8,490 + 6,645 of 15,135 bp | 8,490 + 6,645 |
| Plastid molecules (bp) | 125,122, 136 (k 99) | 102,855, 22,001 |

The held-out candidate's two molecules total 307,026 bp. The 150x mitochondrion's
comparison first reported 42,754 bp matching neither candidate nor plastid; aligned on
their own, those two stretches match the candidate exactly (38,284 of 38,284 and 4,459 of
4,469 bp), so that figure was the whole-molecule alignment's, not the assembly's.

Every molecule is linear and `decisive=false`. The mitochondrion has an 8,268 bp repeat at
twice its depth (and repeats of 135-837 bp); no read spans it, so short reads cannot say
whether the genome is one molecule or the candidate's two. The plastid (no inverted
repeat) stops at an AT-rich inverted repeat (hairpin) whose reads thin to about a tenth of
the plastid's depth. `reconstruction.json` now reports each organelle's `final` status
with these reasons, beside `<organelle>.gfa` and `<organelle>.fasta`.

### Closing the plastid with read pairs

The plastid's linear molecule (125,122 bp) ends at an AT-rich inverted repeat. With mates
kept together, 28 read pairs at 1,000x join its end to its start (insert median 333 bp, the
ends 185 bp apart by them); the reads anchored at either end and their mates spell 194 bp
between them, and the plastid is one circle of 125,300 bp (`complete`, no loose ends, read
pairs only). Against the held-out HiFi plastid (125,509 bp) it aligns over its whole length
in two blocks split at the circle's origin, 114,017 of 114,127 and 11,283 of 11,382 bp
matching. The first full run sampled reads by their whole names, so mates were kept apart:
at the plastid's 11.6%, one joining pair remained and the molecule stayed linear.

This check reused recruited reads (`--target-reads`: the v6 mitochondrial reads and plastid
reads of the earlier recruitment downsampled by read name to 1,000x; 5:59, 1.25 GB). The
mitochondrion (307,479 bp) stays linear: 25 pairs join its ends but put them 569 bp into
each other, without an exact overlap, which is the 8.3 kb repeat's two copies at its ends
rather than missing sequence. No reference bases were filled in either organelle.

### Every configuration of the mitochondrion

With the same reads (`--target-reads`, 3:13, 1.20 GB), the mitochondrion's anchor graph has
five repeats at twice its depth (8,268, 837, 360, 183 and 135 bp) and 32 covers, all with the
same score: reads decide none of the pairings. All 32 are listed. Among them:

| Configurations | Molecules (bp) |
|---|---|
| 1-8 | one linear molecule of 307,353 |
| 9-10 | 160,953 circular + 146,400 linear |
| 17-24 | one circle of 306,900 (a master circle) + 453 linear |
| 25-26 | 160,953 circular + 145,947 circular + 453 linear |
| 11-16, 27-32 | other pairs (94,352 / 127,083 / 128,222 circles with the rest) |

Configuration 25 is the held-out HiFi assembly's two molecules: its 160,953 bp circle matches
HiFi `linear.1` over its whole length base for base. The 453 bp linear piece is the two short
dead-end tips beside the 360 bp repeat with that repeat; it is the loose end that keeps the
largest molecule of the other configurations from closing. Configuration 17 joins both HiFi
molecules into one circle through the repeat.

These lengths were first spelled from the organelle graph's blunt segments, and its second
circle came out 145,947 bp, aligning to HiFi `linear.2` (146,073 bp) with 126 bp missing. The
360 bp repeat has two links with 126 bp overlaps, to the two tips; making the graph blunt takes
those bases off the repeat's segment for every walk, also those that pass it without the tips.
The molecules are now spelled on the assembly graph, each used link's overlap taken off once;
the GFA walks still name the blunt segments, and `gfa_walk_bp` gives what they spell. Now
configuration 25 is the HiFi assembly exactly: 160,953 bp identical to `linear.1`, 146,073 bp
identical to `linear.2` (two blocks split at the origin, 97,490 and 48,583 bp, no difference),
beside the 453 bp tip piece; configuration 17's master circle is their sum, 307,026 bp.

### Ranking the configurations by the reference runs

Of the eight reference runs, the pairings of ZM4's repeats (runs holding each / runs holding
its entry's flank):

| Repeat | Pairings |
|---|---|
| 135 bp (u16) | u3-/u5+ 8 of 8; u3-/u4- 0 of 8; u1-/u4- 3 of 6; u1-/u5+ 0 of 3 |
| 183 bp (u15) | u0+/u4+ 4 of 8; u0+/u5- 3 of 8; u2-/u4+ 2 of 8; u2-/u5- 0 of 7 |
| 360 bp (u13) | all four pairings 8 of 8 (its flanks include the two short tips) |
| 837 bp (u10) | all four pairings 0 of 7-8 |
| 8,268 bp (u9) | u1+/u7+ 5 of 8; u0-/u1- 5 of 7; u0-/u3- 0 of 7; u3+/u7+ 0 of 8 |

The held-out HiFi assembly's configuration (two circles, 160,953 and 146,073 bp) uses the
135 bp repeat's pairings the runs hold (8 of 8 and 3 of 6) and ranks 9th-14th (four
configurations give those two molecules); the 16 ranked 17th-32nd all use the 135 bp pairing
no run holds and none is the HiFi one. At the 8.3 kb repeat the HiFi pairings (5 of 7 and
0 of 8) score as the alternative's (5 of 8 and 0 of 7): no run holds ZM4's arrangement of it,
and the runs do not choose between the two. The 837 bp repeat's pairings are in no run at all.
So the reference prior halves the candidates here and keeps the HiFi configuration among the
remaining 16, but does not pick it. With a 200 bp distance tolerance both pairings of the 135
and 360 bp repeats were "held by 8 of 8"; the tolerance is now 5% (at least 20 bp).

### Genotypes at the reference graph's sites

With the same reads (`--target-reads`, 4:38, 1.24 GB for the whole run), against the held-out
HiFi assembly (an allele is its truth when the HiFi sequence holds that allele's spelling and
no other's):

| | Mitochondrion | Plastid |
|---|---|---|
| Sites (simple bubbles) | 660 | 247 |
| Called | 584 | 231 |
| Major allele = HiFi allele | 569 of 570 (99.8%) | 229 of 229 |
| Sites with a minor allele (>= 10%) | 6 | 0 |
| Segments called present / of them in HiFi | 786 / 777 | 338 / 337 |
| Segments called absent / of them not in HiFi | 170 / 169 | 45 / 45 |

The one discordant mitochondrial site is split 46:54 between its alleles (218 and 253): a
minor allele by the reads, and the HiFi assembly holds the other. Sites where HiFi held
neither allele (14 and 2) or both were not scored.

### The mosaic over the reference runs

| | Mitochondrion | Plastid |
|---|---|---|
| Sites on the backbone walk | 627 of 660 | 243 of 247 |
| Closest runs (weight) | SRR31693268 0.945 | SRR31693269 0.691, SRR31693276 0.309 |
| Viterbi stretches | one (SRR31693268) | one (SRR31693269) |
| Called sites predicted without their call, agreeing | 529 of 584 (90.6%) | 230 of 231 (99.6%) |
| The same against the held-out HiFi | 513 of 564 (the reads' calls: 563 of 564) | 223 of 224 (224 of 224) |
| Imputed sites (mean probability) | 76 (0.90) | 16 (0.84) |
| Imputed sites the HiFi shows, right | 5 of 5 | 7 of 8 |

ZM4's mitochondrion is closest to one run throughout, with no switch, but differs from it at
about a tenth of the sites the reads call: the reads' own calls are right where the closest
run's are not, and imputation from the runs is right about nine times in ten there. Its plastid
is all but predicted by two near-identical runs. 4:06, 1.29 GB for the whole `--target-reads`
run.

### Circular representative, VCF and write-back

Same reads (`--target-reads`, 11:03, 1.38 GB). The mitochondrion's representative is now
configuration 3 of 32: one circle of 307,026 bp, the fewest circular molecules; it aligns to the
held-out HiFi assembly in four blocks covering both its molecules, every base identical (HiFi
`linear.1` 0-45,952 and 45,769-160,953, `linear.2` 0-97,673 and 97,490-146,073), i.e. the two HiFi
molecules joined through the 8.3 kb repeat; the 453 bp tip piece is listed unplaced. All 32
configuration walks of the GFA spell their FASTA lengths (the first version took a repeat's
pieces by name and missed the one its overlapping neighbour owns; `pan` refused the walk). The
VCF has 154 mitochondrial and 181 plastid records. Written back, the mitochondrial graph grows
from 350,586 to 374,365 bp (23,779 bp of ZM4 sequence no reference run holds, the 15 kb in
two pieces among it) and the plastid graph from 128,321 to 133,668 bp.

### End to end from the WGS reads (e4c7eec)

`pan reconstruct` on ZM4's 155,554,892 WGS records with `--target-depth 1000`: 1:07:10, 4.73 GB
peak (counting the 31-mers took 2,734 s on a slow disk; 646 s before). Recruitment kept every
mitochondrial read (895,480) and 856,699 of 7,354,090 plastid reads, mates together.

| | Mitochondrion | Plastid |
|---|---|---|
| Status | `circular_configuration_chosen` | `complete` |
| Representative | one circle of 307,026 bp (configuration 3 of 32) + 453 bp tip piece | one circle of 125,243 bp (de novo at k 77) + 92 bp piece |
| Against the held-out HiFi | four blocks covering both HiFi molecules, every base identical | two blocks, 120,227 of 120,373 and 5,016 of 5,079 bp matching (HiFi 125,509 bp) |
| Sites called / major = HiFi | 584 / 569 of 570 | 232 / 230 of 230 |
| Segments present in HiFi / absent not in HiFi | 777 of 786 / 169 of 170 | 336 of 337 / 45 of 45 |
| Mosaic: agreement / imputed right where HiFi shows | 90.6% / 5 of 5 | 99.6% / 6 of 7 |
| VCF records | 154 | 181 |
| New bases written back into the species graph | 23,779 | 5,307 |

The plastid closed in the de novo graph itself this time (the reads with mates kept together;
no read-pair closure was needed) at k 77, 266 bp shorter than the HiFi plastid: the smaller k
may have compressed its AT-rich inverted repeat. No reference bases were filled anywhere.

### The plastid's 266 bp

The end-to-end plastid circle (125,243 bp) lacks 266 bp of the held-out HiFi plastid at four
AT-rich places: 74 bp at HiFi 102,965 (microsatellite with short inverted repeats), 72 bp at
103,708 (two of four 36 bp tandem copies), 63 bp at 120,842 (three of seven 21 bp copies) and 57 bp
at 120,373-120,430. Each place, with the nearest unique 20 bp on either side, spans 128-247 bp:
no 150 bp read holds one whole (one read at the 128 bp place). Unit 17-mer depth over the flanks'
gives 3.6 copies at the 36 bp array (HiFi 4, ours 2) but 2.1 at the 21 bp array (HiFi 7, ours 4),
where AT-rich coverage is low. Read pairs spanning each place (288-1,063 of them) have, on our
circle, median inserts 6-17 bp shorter than all pairs' (328 bp) at three places and 16 bp longer
at the fourth, while two control places give 10 and 18 bp longer: against the controls our circle
is some 20-30 bp short at three places and not short at the fourth, not the 57-74 bp HiFi has. The
short reads do not measure these arrays' lengths; the pairs suggest ZM4's WGS DNA has somewhat
more than our circle and less than the HiFi assembly there (plastid microsatellites vary in
length between plants and within one).

### Marking stretches of uncertain length

On the end-to-end reads (`--target-reads` of that run, 2:23, 1.21 GB) the plastid circle has
7 marked stretches, among them all four places where it differs from the HiFi plastid (circle
20-187, 412-496, 107,981-108,116 and 108,541-108,800); the other three (176-298 bp) are as long
as in HiFi but no read spans them. The mitochondrion has 3 (two inside the 8.3 kb repeat and a
112 bp one). A first version also marked stretches whose spanning pairs had 15 bp longer inserts
and those spanned at under a fifth of the expected rate: 15 plastid stretches, 11 of them as long
as in HiFi.
