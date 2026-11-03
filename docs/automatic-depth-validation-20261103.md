# Automatic HiFi nuclear-depth validation, 2026-11-03

The Lathyrus development case now runs without a supplied nuclear depth. This is not a universal low-coverage guarantee or independent validation of the true nuclear coverage.

## Actual end-to-end result

- Input: the same single ERR13245275.fastq.gz used by the prior default and manual-depth runs; no partially downloaded runs included.
- Formal source location: /home/jiazc/data24t/work-live/PhytoOrga_AI_c/ovasm.
- Automatic depth estimate: 5; nuclear source: reads; organelle gate: 15.
- Mitochondrial output: one circular 446,897-bp representative, decisive=false.
- CIGAR-based windowed Oatk comparator reference-span coverage: 446,899/446,899 = 100%; aligned identity: 99.999553%. The initial approximate PAF comparison (without minimap2 -c) reported 99.991049% and is superseded.
- The mitochondrial FASTA is byte-identical to the earlier --nuclear-depth 4 output.
- No uncovered reference-span intervals remain in the CIGAR-based reanalysis. This is span coverage, not proof of base identity or correct structural conformation. Evidence: ../.agent-work/comparison-cigar-auto.json.
- Wall time 9:14.84; peak RSS 1,949,660 KiB; exit status 0. Shared-machine timing is descriptive.
- Tested executable SHA-256: f1ccd141b12a96161daa0509db44354d5a47e9a11cc44132f227c5c451bee92f.

The tested formal-source snapshot included concurrent uncommitted bait/recruit/pipeline work. Its exact source and binary hashes, original command, status and comparison are retained under:

/home/jiazc/data24t/work-live/PhytoOrga_AI_c/.agent-work/ovasm-auto-depth-20261103

This is a development validation snapshot, not the final immutable comparative-paper release. Later fixes and any final benchmark version need their own recorded checksums and passing tests.

## What changed and what can be replayed cheaply

The main estimator now limits candidate windows to local maxima rather than allowing a broad distribution shoulder to pass only because it clears the genome-size floor. On the recorded integer read-depth distribution, the old estimator gives 85 (window-implied size 20,030,541 bp), while the local-maximum estimator gives 5 (212,305,806 bp). A separate guarded low-depth integer-peak fallback exists but was not triggered in this run. Therefore this result must not be described as an observed fallback to 4 or preservation of judged_peak_depth=85.

The real integer frequency/bases fixture is tests/fixtures/lathyrus-nuclear-depth.json. scripts/replay_nuclear_depth.py replays the old and current estimates without rereading the raw FASTQ or using a reference assembly. Its Rust counterpart checks the production estimator. This is the appropriate first regression for future estimator changes; only changed recruitment behavior requires another end-to-end read run.

The initial audit of 173 old reports cannot establish that all other samples are unchanged. Local-maximum filtering necessarily changes the candidate outcome for the recorded Polytrichum formosum case as well as Lathyrus; the remaining reports lack enough integer-depth/bases data for complete replay. Broader low-depth validation remains necessary.

## Paper experiment

The planned nuclear-depth benchmark compares OVASM with PMAT1 and Oatk, using independent genome-size denominators, matched subsets and repeated seeds. PMAT1 already established the low-coverage HiFi use case. The current Gb-based Anagallis pilot is engineering validation, not evidence of a nuclear-depth limit. See the local manuscript LOW-DEPTH-HIFI-PLAN.md and the server benchmark inventory for the preregistered endpoints and known installation/data limitations.
