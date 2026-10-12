# Recruitment lookup performance

OVASM uses a compact presence index in front of the exact bait k-mer map. An unset bit proves that a key is absent. A set bit still requires the original exact target-mask lookup; collisions cannot recruit a read or change its target. Changes to the bait invalidate the index, which is rebuilt before the next scan. Read classification, NUMT flagging, extension, saturation, downsampling and output order are unchanged.

## Server development measurements

Measurements used complete original ZM4 HiFi FASTQ records and identical seed libraries, with eight worker threads and the same release build profile on `forageseed-235`. Baseline source: `b529d67`. ABBA execution order provided two one-pass measurements per version and subset. These are development checks on a shared server, not blind validation or confidence intervals.

| Input | Baseline median wall time | Indexed median wall time | Baseline/indexed ratio | Baseline/indexed maximum RSS |
|---|---:|---:|---:|---:|
| First 100,000 original reads, 1,727,183,659 bp, uncompressed | 43.12 s | 28.775 s | 1.499 | 520.16 / 536.76 MiB |
| Same 100,000 complete records, gzip reencoded | 30.485 s | 27.580 s | 1.105 | 521.07 / 536.38 MiB |
| First 10,000 original reads, 172,437,825 bp, uncompressed | 8.035 s | 7.935 s | 1.013 | 519.32 / 534.96 MiB |
| 5,000 recruited MT + 5,000 recruited PT reads, 166,397,873 bp, uncompressed | 9.64 s | 11.31 s | 0.852 | 563.93 / 608.89 MiB |

The larger original-read prefix showed about 33% shorter wall time. The small target-rich subset was slower: index construction and the extra probe can cost more than they save when most reads match. Preserve this negative result when assessing suitability for enriched inputs. The extra index used about 16.6 MiB in the larger-prefix comparison.

The gzip-reader check showed about 9.5% shorter median wall time. Baseline repetitions were 34.29 and 26.68 s; indexed repetitions were 28.35 and 26.81 s. The substantial baseline variation and only two repetitions per version limit the inference. The same selected complete FASTQ bytes were reencoded with gzip level 6; this exercises gzip reading but does not reproduce the original file's compressed blocks. Input preparation took 407.805 s, separate from the repeated scan measurements. Uncompressed and gzip experiments ran at different times on a shared server; their absolute timings must not be used to infer that gzip is faster than uncompressed reading.

All twenty subruns produced byte-identical final FASTQ files and equal scientific recruitment reports, excluding only output-directory paths, measured time and RSS. The smaller-subset checks also requested three iterations, exercising bait extension and saturation; one subset saturated after two iterations in both versions. Regression tests covered all eight target masks, shared sequences, foreign flanks, reverse complements, Ns and several k-mer sizes: 172 passed, zero failed, three existing environment-dependent diagnostics ignored.

Run IDs: `ovasm-lookup-real-read-comparison-v1`, `ovasm-lookup-real-read-large-prefix-v2` and `ovasm-lookup-real-read-gzip-prefix-v3`. Inputs, commands, time/RSS measurements, output comparisons, test logs and the tested source diff are retained under the server project task directory `ovasm-batch-existing20/performance-recruit`.

The default depth gate was disabled in both timing arms of all subset experiments. The gzip experiment includes decompression of the fixed prefix; none of these measurements includes full-input depth inference or whole-pipeline assembly. The completed original ZM4 run took 17,373.18 s: recruitment accounted for 17,338.74 s and the subsequent MT/PT stages for 28.15 s. It scanned the full input seven times for recruitment; the first depth-gate construction added two full-input passes. This lookup change accelerates a query path while retaining those scans. No complete-input acceleration factor is established.

The original 25-job benchmark remains frozen at its original source and binary. Its first-round results, wall times and memory measurements are not replaced by these development comparisons.
