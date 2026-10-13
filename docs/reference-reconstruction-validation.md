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
