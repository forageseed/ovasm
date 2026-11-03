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

The initial audit of 173 old reports cannot establish that all other samples are unchanged. Uncorrected local-maximum filtering necessarily changed the recorded Polytrichum formosum case as well as Lathyrus. The empty-bin boundary correction and actual full-input test below address that specific regression. The other old reports lack enough integer-depth/bases data for complete replay; no claim of universal invariance follows from their binned histograms.

## Paper experiment

The planned nuclear-depth benchmark compares OVASM with PMAT1 and Oatk, using independent genome-size denominators, matched subsets and repeated seeds. PMAT1 already established the low-coverage HiFi use case. The current Gb-based Anagallis pilot is engineering validation, not evidence of a nuclear-depth limit. See the local manuscript LOW-DEPTH-HIFI-PLAN.md and the server benchmark inventory for the preregistered endpoints and known installation/data limitations.

## Cross-species boundary and strong-peak checks

The unified formal-source executable for these checks has SHA-256
`449c7f6acaa6c88029caf4aff760814e9e1681fcbb891f55e397a5b1ea57aaee`.
Its source snapshot and complete commands are recorded in
`../.agent-work/ovasm-five-fixes-20261103/`. These are whole-input runs with default
nuclear-depth parameters, not downsampled extrapolations.

### Polytrichum formosum: censored low-depth boundary

- Input: ERR12721079.fastq.gz, 340,868,549 sequenced bases.
- Both the historical estimator and corrected estimator replay to 3x and a
  23,670,683-bp window-implied size on the exact new integer histogram.
- The actual default run reports nuclear depth 3, organelle threshold 9, source
  `low-depth read boundary`, and a warning that lower depths are censored. This
  preserves the threshold without claiming independently known nuclear depth.
- The empty log bin between integer depths 3 and 4 no longer suppresses a decreasing
  boundary population merely because its smoothed count includes both neighbours.
- Mitochondrial recruitment: 59 reads / 798,284 bases; plastid: 576 reads / 6,762,887
  bases. Both recruited FASTQ files are byte-identical to the historical hifi-v1
  outputs (SHA-256: mt `9c6ea1aedf78bdadd0dec77bbfaa08c50d18200d3fa7c56fc5579df37e420abe`,
  pt `2d9ff7277dc7e336525d0063b02fe632881eb66c00b88eb180f1c6296718b81f`).
  Thus this boundary fix did not add recruited reads in this regression case.
- End-to-end exit 0, wall 1:00.04, peak RSS 879,740 KiB. Output: 60,750-bp linear
  mitochondrial representative and 126,243-bp circular plastid; both decisive=false.
- Exact fixture: `tests/fixtures/polytrichum-formosum-nuclear-depth.json`;
  replay with `python3 scripts/replay_nuclear_depth.py` followed by that path.

Assembly/reference coverage is not used as truth for nuclear depth. Read-set
identity is a stronger regression check here than an unchanged assembly length,
but it still does not certify the physical origin of every historical read.

### Polytrichum commune: strong read peak despite a 2x k-mer mode

- Input: ERR9863245.fastq.gz, 10,185,557,573 sequenced bases.
- Historical and corrected estimators both replay to 24x and 264,538,883 bp of
  window-implied genome size. The actual run reports source `reads`, nuclear 24,
  threshold 72, k-mer mode 2, and no nuclear-depth warning. The 2x k-mer mode does
  not replace a supported read-depth estimate.
- Mitochondrial recruitment: 3,007 reads / 28,542,545 bases; plastid: 12,175 reads /
  119,060,334 bases; both saturated. The recruited files are byte-identical to
  historical hifi-v1 (SHA-256 mt
  `d745aea7a9f2b835811440252fe64286111dba25c8e8c475773c91dc200fe4fc`,
  pt `a0434cd2bba6c18585d21002cff81b17b559d80cb57312830744b27b1f0a607b`).
- End-to-end exit 0, wall 6:01.48, peak RSS 1,690,512 KiB. Mitochondrion:
  114,831-bp circle, decisive=true; plastid: 126,254-bp circle, decisive=false.
- Exact fixture: `tests/fixtures/polytrichum-commune-nuclear-depth.json`, usable
  with the same replay script. Together with Lathyrus and P. formosum these cover
  the repaired broad shoulder, the censored low-depth boundary, and a supported
  read peak that must not follow the low k-mer mode. They are three regression
  cases, not proof of unchanged behavior across all 173 historical reports.


### Reference-span checks for the two regression species

CIGAR-based comparison of both completed runs is recorded in ../.agent-work/ovasm-five-fixes-20261103/comparison.json, with explicit reference IDs cbPolForo1 and cbPolComm4. P. commune covers 100% of both Oatk reference spans. P. formosum covers 52.9162% of the mitochondrial span and 100% of the plastid span. Its identical historical recruitment therefore validates preservation of the low-depth boundary behavior, not complete mitochondrial recovery at 3x. Reference-span coverage is not independent structural certification.
