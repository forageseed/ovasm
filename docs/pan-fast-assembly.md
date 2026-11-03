# Optional fast assembly from a deterministic 20% WGS subset

This is an experimental **primary-assembly** workflow, not the default and not a replacement for whole-input minor-allele or low-coverage analysis. It makes no claim of being faster than GetOrganelle: the measured OVASM workflow reconstructs both organelles, whereas the available GetOrganelle timing covers a plastid workflow with different input selection. The separately measured assembly-only ratio must not be relabelled whole-pipeline performance.

## Reproduce the two timed stages

On Linux with g++ and the installed zlib runtime, compile the standalone sampler. Its small declaration block uses the stable zlib public C ABI because the validated machine lacked development headers. No new dependency or data download is required there.

```sh
g++ -O3 -std=c++17 scripts/pan_pair_sample.cpp -Wl,-l:libz.so.1 -o /absolute/new-work/pan_pair_sample
/usr/bin/time -v -o /absolute/new-work/sampling.time.txt \
  /absolute/new-work/pan_pair_sample R1.fastq.gz R2.fastq.gz \
  /absolute/new-work/pair20_R1.fastq /absolute/new-work/pair20_R2.fastq \
  /absolute/new-work/sampling.json
/usr/bin/time -v -o /absolute/new-work/pipeline.time.txt \
  /absolute/formal-build/ovasm pan reconstruct \
  --sample sample_pair20 --threads 8 --target-depth 1000 \
  --reads /absolute/new-work/pair20_R1.fastq \
  --reads /absolute/new-work/pair20_R2.fastq \
  --reference mitochondrion=/absolute/panel/mitochondrion.gfa \
  --reference plastid=/absolute/panel/plastid.gfa \
  --out /absolute/new-work/assembly
```

The work directory must already exist. Sampler output files must not exist; creation is exclusive and existing files are never overwritten. On any parse, pair, input, or output error the sampler exits nonzero. Created output files can then be partial and are retained with an explicit stderr warning; move them aside before retrying. Require **both exit 0 and a successfully parsed report**, not merely the presence of FASTQ files.

Both stages use the same paired input source. The sampler checks every pair name after taking the header token and removing terminal `/1` or `/2`, verifies FASTQ line structure and sequence/quality lengths, and selects both mates together. It writes their original four record lines unchanged. Selection is fixed: FNV-1a 64-bit initialized with `14695981039346656037 XOR 20261103`, followed by the mix64 finalizer, and retain when `hash % 5 == 0`. Salt is **20261103**, fraction **one in five**, deliberately not a new configurable product default. The two mates must be synchronized four-line FASTQ, not wrapped FASTQ. Gzip input is handled with zlib.

Report fields: `pairs_scanned`, `pairs_kept`, `bases_scanned`, `bases_kept` are integer totals across the validated pair stream; `sampling_seconds` is sampler elapsed time; `salt` and `rule` identify the selection; `pair_names_checked_exactly` and `read_lines_unchanged` describe validation. GNU time additionally records elapsed time and RSS. **Total workflow time is sampling wall time plus pipeline wall time.** Peak RSS for these sequential stages is their maximum, not their sum. Compilation and diagnostic alignment are reported separately from assembly runtime.

## Compare sequence quality

```sh
python3 -I scripts/pan_pair_sample_evaluate.py \
  --run pair20=/absolute/new-work/assembly \
  --run full=/absolute/full-input/assembly \
  --reference mitochondrion=/absolute/diagnostic/mitochondrion.fasta \
  --reference plastid=/absolute/diagnostic/plastid.fasta \
  --minimap2 /absolute/minimap2 --threads 4 \
  --out /absolute/new-work/comparison
```

This helper refuses an existing output directory and records minimap2 PAF/logs and `quality-comparison.json`. It uses asm5, secondary alignments, and exact CIGAR operations: reference aligned/exact base counts are interval unions, so overlapping alignments are not counted twice. Identity is alignment-weighted and can include secondary alignments. Optional `--gaps-json` accepts the existing diagnostic schema `existing_HiFi_candidate_sequence_comparison -> reference_name -> unmatched_intervals`, to quantify exact recovery of historical missing intervals. A missing reconstruction is reported pending, not passed. A HiFi assembly of another individual is a diagnostic reference, not independent truth for every base or topology. Examine final status, uncertain regions, representative/alternative configurations and genotype tables as well as coverage.

## Measured ZM4 experiment, 2026-11-04 (Asia/Shanghai)

The same frozen formal-source binary, eight threads, both organelles, the same eight-run pan, `--target-depth 1000`. Whole input is the authorized usable prefix of SRR13333747: 77,777,446 pairs, 23.333 Gb. The full and sampled runs ran on a shared server; these are single observations, not a randomized replicated timing study.

| Measure | Full WGS | 20% experiment |
|---|---:|---:|
| Input pairs to reconstruction | 77,777,446 | 15,556,683 |
| Extra sampling wall cost | — | 290.28 s |
| Reconstruction wall | 1,304.85 s | 263.83 s |
| Total including sampling | **21m45s** | **9m14s** |
| Pipeline max RSS | 4,726,864 KiB | 4,726,380 KiB |
| Sampling max RSS | — | 9,064 KiB |
| Internal 31-mer counting | 743.20 s | 100.38 s |
| MT representative | 307,026 bp + 453 bp tip | 307,026 bp + 116 bp tip |
| Both diagnostic MT references, exact union | 100% | 100% |
| Historical 8,490 + 6,645 bp gaps | all recovered | all recovered |
| PT representative circle | 125,243 bp | 125,300 bp |
| PT exact reference bases | 125,243 | 125,300 |
| PT uncertain length regions | 7 | 7 |
| MT heteroplasmic sites | 6 | 8 |
| PT called sites | 232 | 231 |

Observed speed ratio including sampling is **2.35x** relative to the same OVASM full-input run. It is not a cross-tool speed claim. MT remains nondecisive; PT remains `closed_length_uncertain`. Neither fills sequence from the reference. MT minor alleles crossed the 10% decision boundary in both directions; PT site238 changed from a major call to no-call. MT segment calls also changed (absent170→185, low34→19). Therefore the result is **not noninferior for all downstream analyses**. Preserve full-depth input and conclusions for those analyses. The PT's increased length is not proof all repeat lengths are correct.

The original measured sampler lacked exclusive output creation; the packaged sampler fixes that behavior and passed exact output equivalence on the validated fixture, deterministic repetition, existing-output refusal without modification, and mismatch rejection with partial-file warning. The full 23Gb experiment was not repeated merely to measure this I/O safety change. The packaged evaluator exactly reproduced the experiment's metrics.

Evidence directory: `/home/jiazc/data10t/ovasm-dev-work-20261103/pan-wgs-pair20-experiment/`. It retains input sizes/mtimes, binary SHA256 and formal source manifest, native source, commands, stage exit codes, time/RSS, sampling totals, complete PAFs, interval coverage, every changed genotype call, and tests. The full control remains at `/home/jiazc/data10t/ovasm-dev-work-20261103/pan-evidence-validation/wgs`.

Frozen experimental binary SHA256: `449c7f6acaa6c88029caf4aff760814e9e1681fcbb891f55e397a5b1ea57aaee`. Formal source hashes: `/home/jiazc/data10t/ovasm-dev-work-20261103/pan-evidence-validation/source-after.json`.
Packaged sampler source SHA256: `d53aadeb553be106beccf5e4a81f56411345080231f674a08e04da73725fdd24`.
