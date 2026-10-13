# ZM4 HiFi: OATK and OVASM configurations

The completed same-input runs used 7,352,539 complete original HiFi records,
134,060,822,653 bases and eight threads. They ran sequentially on a shared server;
these are measured runs, not an isolated randomized performance experiment.

| Measurement | OATK 1.0 | Released OVASM 74e0cf7 |
| --- | ---: | ---: |
| Command/controller wall time | 56,540.957 s | Public command 12,612.001 s |
| Sampled sum of process-group RSS | 310.416 GiB | 1.983 GiB |
| GNU maximum RSS | 186,622,060 KiB | 2,216,348 KiB |
| Main mitochondrial FASTA | One 307,026-bp graph cycle | Two graph cycles: 160,953 and 146,073 bp |
| Main plastid FASTA | 125,509 bp | 125,509 bp |
| Additional plastid-labelled record | 22,834-bp linear sequence | Absent from primary FASTA |

OATK includes HMM annotation of its whole-input graph. Its first annotation log
arrived 13,296.578 seconds after the start. Thus the observed pre-annotation
interval was about 3 h 42 min, close to OVASM's entire 3 h 30 min command. Buffered
log observations cannot establish exact exclusive internal stage times. Summed
process RSS and GNU maximum RSS are different measurements; sums may count shared
pages in more than one process and are not physical memory consumption.

OVASM also retained a **307,026-bp single-cycle alternative** in the same
`linearization.json`. Spelling that alternative in the original graph gives the
OATK mitochondrial sequence exactly, modulo circular rotation/reverse complement.
No second assembly or ranking change was required to export it.

The graph has unique segments A=152,685 bp and B=137,805 bp and a shared repeat
R=8,268 bp. The alternatives are A-R plus B-R, or A-R-B-R. The following are local
HiFi-supported repeat contexts, not observations of entire physical chromosomes:

| Context | Reads |
| --- | ---: |
| A-R-A | 398 |
| B-R-B | 380 |
| A-R-B | 322 |
| B-R-A | 411 |

OVASM's internal path scores were -1.324517492 for the two-cycle configuration
and -1.457995990 for the single-cycle configuration, a margin of 0.133478498.
The report correctly retained `decisive=false`. All four junction types have
hundreds of supporting reads. The mitochondrial canonical 31-mer multisets,
including graph-cycle boundary words and copy multiplicities, were equal across
the two programs' primary outputs. This is sequence-content compatibility; it
does not select one physical configuration. The main plastid sequences were
exact circular rotations. OATK's extra 22,834-bp record contributes 4,104 distinct
31-mers absent from OVASM's primary plastid; its source is not established here.

The separate packed input-replay development candidate e8a85b3 also completed
on the same full HiFi input. All compared outputs were equivalent, but its
public command took 23,549.547 seconds (6 h 32 min), versus 12,612.001 seconds
for the released baseline: **1.867 times longer**. This candidate did not meet
the acceleration objective and was not integrated into the released assembler.
Prefix speed measurements must not be extrapolated to this whole-input run.

Neither assembler is an independent biological truth for the other. Graph cycles
do not establish physical molecule count or topology; strict and relaxed truth
scores remain unevaluated.
