# Existing WGS first-attempt benchmark

Twenty complete existing Illumina WGS inputs from the original 25-job queue.
Frozen OVASM b529d67, eight threads/job, two concurrent jobs; shared server.
First-attempt limits: 96 GiB sampled process-group RSS or 24 h wall/job.
No replacement by later repairs or selected successful reruns. The other five
HiFi/ONT/Pore-C entries are outside this WGS subset; this is not a terminal
summary or success rate for all 25 jobs.

| Job | Outcome | Wall min | Sampled peak RSS sum GiB | GNU max RSS KiB | MT / PT records |
|---|---|---:|---:|---:|---:|
| wgs_01_CRR2683883 | completed | 51.45 | 2.449 | 2615388 | 49 / 3 |
| wgs_02_CRR2683762 | completed | 35.22 | 2.364 | 2523872 | 46 / 4 |
| wgs_03_CRR2683760 | completed | 43.99 | 2.400 | 2561244 | 51 / 3 |
| wgs_04_CRR2683763 | completed | 47.17 | 2.467 | 2630324 | 49 / 4 |
| wgs_05_PI502453 | failed / resource_limit_time | 1440.00 | 6.968 | not recorded | — / — |
| wgs_06_CRR1151719__CRR1151719 | failed / resource_limit_memory | 925.64 | 96.063 | not recorded | — / — |
| wgs_07_F152__CRR1143781 | failed / resource_limit_memory | 574.72 | 96.005 | not recorded | — / — |
| wgs_08_CRR2683993 | completed | 25.81 | 2.132 | 2303496 | 45 / 3 |
| wgs_09_CRR2683802 | completed | 32.57 | 1.809 | 1914500 | 48 / 3 |
| wgs_10_CRR2683983 | completed | 85.49 | 1.917 | 2078832 | 49 / 5 |
| wgs_11_CRR2683756 | completed | 36.32 | 2.554 | 2816112 | 49 / 3 |
| wgs_12_WL_354 | failed / resource_limit_time | 1440.01 | 86.312 | not recorded | — / — |
| wgs_13_CRR1151718__CRR1151718 | failed / resource_limit_memory | 613.46 | 96.002 | not recorded | — / — |
| wgs_14_CRR2683867 | completed | 107.62 | 6.737 | 7063760 | 54 / 3 |
| wgs_15_CRR2683921 | completed | 56.11 | 6.212 | 6742232 | 40 / 4 |
| wgs_16_PI478450_1 | failed / resource_limit_time | 1440.01 | 11.117 | not recorded | — / — |
| wgs_17_CRR2683963 | completed | 47.23 | 2.448 | 2619760 | 51 / 4 |
| wgs_18_CRR2684001 | completed | 78.42 | 7.547 | 7981160 | 46 / 37 |
| wgs_19_WL_237 | completed | 189.59 | 5.062 | 5786696 | 54 / 6 |
| wgs_20_CRR2683845 | completed | 55.40 | 2.597 | 2853732 | 44 / 3 |

Execution completed: 14/20 (70%); failed: 6/20. Among completed commands, median
wall 49.34 min and median sampled peak RSS sum 2.458 GiB. Three time-limit
failures and three memory-limit failures remain in the denominator. Most WGS
outputs are fragmented: MT decisive 0/14, PT decisive 1/14. An accepted internal
k or nonempty FASTA is not complete-assembly acceptance. Record counts do not
establish physical chromosome counts. No same-sample independent truth was
supplied, so strict/relaxed truth scores remain null.

RSS-sum sampling is once per second and can miss short peaks or double-count
shared pages. GNU maximum RSS is a different per-process measurement. Wall
time includes contention on the shared server. Failed jobs did not always
write a GNU-time terminal record; those cells stay unrecorded.

[Machine-readable measurements and original receipt SHA256](zm4-benchmark-resources.json)
retain all first-attempt costs and the independent historical full-HiFi runs.
