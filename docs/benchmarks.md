# OVASM benchmark records

## Compared with other tools

<p align="center"><img src="docs/benchmark.svg" alt="ovasm compared with Oatk, HiMT and PMAT" width="900"></p>

| | ovasm | Oatk | HiMT | PMAT |
|---|---|---|---|---|
| Success (63 chloroplasts) | 63/63 | 63/63 | 0/63, see note | 62/63 |
| Accuracy (Col-0 mitochondrion) | complete circle, 100% | 99.0% in 4 pieces | 99.8% in 15 pieces | not run |
| Time (Col-0, full reads) | about 3 min | 31 min | 98 min | – |
| Memory (63 chloroplasts, median) | 58 MB | 1,093 MB | 297 MB | 1,129 MB |
| Platforms | one program; built and tested on Linux and Windows | conda environment | conda environment | container |

All numbers were measured by us on the named data. Historical summaries remain below;
new full-input resource and per-job WGS tables are included in this repository. Read them with these limits in mind:

- **ovasm is not shown to be better across species.** In four rounds of frozen blind tests on new species
  (pre-registered thresholds, identity at least 0.9999, 12 tasks each: mitochondrion and chloroplast of 6 species),
  ovasm passed 3/12, 6/12, 8/12 and 7/12; the first round (3 species) passed 2/6. The other tools were not scored
  under the same thresholds on those species, so there is no head-to-head blind comparison. Fixes made while looking
  at those data raise the replays (for example 10/12 on one round), but replays are not blind tests. The figure uses
  data we had already seen.
- **Col-0 mitochondrion:** the seed was Col-0's own reference genome, which makes recruitment easy; it affects which
  reads are recruited, not how the graph is built. When Oatk and HiMT were given the 150x reads that ovasm recruited,
  both found 8/8 junctions in 42 s and 84 s (plus about 200 s of recruitment); Oatk then gave 3 pieces and HiMT 14.
- **Chloroplast cohort:** the "truth" is the cohort's own published assembly, not an independent gold standard, and it
  was a post-hoc comparison. HiMT is a mitochondrial assembler: its 63 runs failed because it did not produce the
  chloroplast output we asked for, which says nothing about its quality. ovasm failed on 2 samples at its first k and
  was fixed by retrying with a smaller k, a rule that is now built into `ovasm run`.
- **Time and memory:** times were measured on a busy shared server and include contention. The historical
  Col-0 timing does not describe ZM4's iterative recruitment; full ZM4 costs are reported below. The historical
  Col-0 run peaks at 229 MB for recruitment and stays under 170 MB for the later
  steps; the heaviest steps are seed-free `discover` (6.9 GB) and hybrid correction of ONT reads (2.3 GB).
- **Platforms (historical check):** on Windows the program builds in 93 s (17 MB), its 141 tests pass, and a full run on the example reads
  takes 5 s. macOS has not been tested.

### Existing Illumina WGS cohort (first attempts)

Of the original 25-job queue, the **20 WGS jobs are terminal: 14/20 executed
successfully, 6/20 failed** (three at 96 GiB sampled RSS, three at 24 h).
For completed commands, median wall time was **49.34 min** and median peak
RSS sum **2.458 GiB**. These are execution statistics, not complete-assembly
success: **0/14 MT and 1/14 PT were internally decisive**, and independent
truth scores were not evaluated. Failed jobs remain in the denominator and
in the [per-job time/memory table](docs/existing-wgs-benchmark.md).
This is the frozen `b529d67` first-attempt WGS subset, not a new full-queue
or latest-commit benchmark.

### Full-input alfalfa ZM4 HiFi: time, memory and graph evidence

Same **7,352,539 complete HiFi reads (134.061 Gb), eight threads**, shared server;
both runs completed. These are observed run costs, not isolated randomized timings.

| Measurement | OVASM `74e0cf7` | OATK 1.0 |
|---|---:|---:|
| Full public command / OATK command-controller wall | **210.20 min (3 h 30 min)** | **942.35 min (15 h 42 min)** |
| OVASM controller including post-run auditing | 214.23 min | included in OATK wall above |
| Sampled peak sum of process-group RSS | **1.983 GiB** | **310.416 GiB** |
| GNU maximum RSS | 2,216,348 KiB | 186,622,060 KiB |
| Main MT sequence representatives | 160,953 + 146,073 bp; 307,026-bp alternative retained | 307,026 bp |
| Main PT sequence | 125,509 bp | 125,509 bp; additional 22,834-bp linear record of unresolved origin |

OATK includes whole-input HMM annotation. Start to first annotation-load log was
**221.61 min**, close to OVASM's complete command. **942.35 min is not pure graph-building time.**
RSS sums can count shared pages more than once; GNU maximum RSS has a different scope.
The prior 96-GiB OATK failure remains a separate attempt.
[Stages, provenance and limits](docs/hifi-configuration-comparison.md).

These costs belong to the named historical commits. The latest junction repair retained
supported alternatives without changing the principal sequences, but no new complete
134-Gb timing was run. Historical opt-in cache candidate `e8a85b3` took **392.49 min**,
**1.867 times longer** than the baseline, with equivalent outputs; prefix speedups
must not be extrapolated to full input. `--recruit-cache` is opt-in.

![Homology-aligned OVASM and OATK graphs with HiFi, ONT and nuclear-context evidence](docs/figures/OVASM-OATK-homology-matched.svg)

[SVG](docs/figures/OVASM-OATK-homology-matched.svg) · [PDF](docs/figures/OVASM-OATK-homology-matched.pdf) ·
[Per-link read support](docs/figures/zm4-link-read-support.tsv) · [Homology coordinates](docs/figures/zm4-homology-correspondence.tsv) ·
[Figure source and methods](docs/figures/README.md#zm4-hifi-ont-supported-gfa-comparison)

Both repaired graphs retain the shared backbone and local alternatives J1/J2, supported
by HiFi and ONT. OATK additionally retains nuclear-context branch u12151: 22/32
branch-compatible complete ONT records meet the strict bilateral nuclear-flank criterion.
Local edge support is not source certification. Without independent same-individual truth,
relaxed/strict scores remain unevaluated; graph cycles do not certify physical topology.

