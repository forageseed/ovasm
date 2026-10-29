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
same bait (as `recruit` with `--max-bases` 150 x 350,586). It never settled: the
mitochondrial reads went from 856,248 (pass 2) to 940,026 (pass 4) and 4,630,842
(pass 5), then to 28.8 M at pass 60, and the plastid's from 7.35 M to 48.4 M. The
150x subsample was then mostly nuclear; the reconstructed mitochondrion fell apart
(largest component 9,753 bp, 230 components, de novo molecules at most 1.8 kb) and
neither interval was closed (4:37:36, 7.4 GB peak).

A one-pass diagnostic recruitment on the held-out candidate's 20.6 kb around the two
intervals (diagnostic only, never used by `reconstruct`) found them in the library at
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
