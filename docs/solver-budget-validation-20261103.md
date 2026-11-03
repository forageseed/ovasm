# Solver budget and fragmented-graph validation (2026-11-03)

Source: the actual server repository `ovasm/src/linearize.rs`; tests and builds bind to that directory. No source worktree was used. Other concurrent bait/recruit changes were present during compilation, so the generated executable is a working-tree build, not a clean release of one commit. Existing background jobs and their binaries were preserved. Windows/browser verification is not applicable to this server CLI profiling exercise and was not performed.

Changes:
- Stop DFS work promptly after the expansion budget expires; do not keep expanding the automatic cover-size envelope after a truncated search. Preserve truncation/non-decisive metadata.
- Precompute bridge endpoint support, outgoing endpoint membership and partial-path closing candidates; preserve self-loop counting and candidate ordering.
- Reuse the search adjacency index across components and cover sizes.

The 30 linearization tests passed after every optimization. The added zero-budget regression checks that the automatic envelope stops after its first attempted cover size and reports uncertainty.

Real input: `ovasm-worktrees/trials/sr-cohort/runaway-fe948ba/wgs_01_CRR2683883/mitochondrion/{assembly.gfa,evidence.json}` (41,747 segments). GDB launch-and-interrupt samples located repeated full bridge scans in endpoint support, path-end detection, partial-path closing and repeated search-index construction. Attaching to the existing long-running process was denied by ptrace policy; its exact hot stack remains unconfirmed.

An optimized binary built from the formal repository completed `solve linearize --max-steps 1000` on that input in **9.99 seconds**, peak RSS **213,784 KiB**, exit 0. Its report has 43 components, 48 molecules, `search_truncated=true`, `decisive=false`, and `cover_size_limited=true`: execution completion does not establish biological completeness.

A separate run at the default 5,000,000-step budget hit the 60-second external timeout (exit 124), without a final report. The default-budget pathological case therefore remains unresolved; this is not a claim that the previous ~30-hour task is fixed. No speedup ratio against that job is justified because the versions and budgets differ.

Logs, stack samples, output and a copy of the tested executable are under the sibling `.agent-work/ovasm-boundary-fixes-20261103/` directory. Unit-test log: `.agent-work/solver-unit-tests.log`. Earlier debug builds also timed out after 45 seconds and are diagnostic runs, not successful validations.

## Final default-budget verification (2026-11-04 Asia/Shanghai)

The earlier default-budget timeout above is superseded by this verified run. Further profiling located candidate-report DNA spelling/reverse-complement allocation and whole-walk allocation inside unconstrained DFS. Commit `39e1a19` deduplicates covers before describing them and computes lengths from indexed overlaps; commit `4b6499b` skips full walk construction when no phase constraints are active. The default search budget remains **5,000,000**, unchanged. These are equivalent-work optimizations, not a reduction of the default budget.

Added tests compare indexed lengths with actual spelling for linear/circular walks, reverse orientation, duplicate edges with different overlaps and self-loops; another test compares complete candidate order and truncation flags between fast and forced full-check paths across 0/10/10000 budgets and 1/2/3-molecule covers. All **32** linearization tests passed. The final frozen working-tree binary passed **353 tests, 0 failures, 6 ignored**; full-source/Cargo before/after hashes matched. Concurrent WIP is recorded in the build manifest; the executable is not represented as a clean release tag.

Frozen binary SHA-256: `19e78482fef473a265ca4cd248267ab460cd672849fbbde767f084a1dca66138`.

On the same 41,747-segment real graph/evidence, `solve linearize` with default arguments (no `--max-steps` override) exited naturally with code **0**, in **15.73 seconds**, peak RSS **213,616 KiB**. The external 300-second deadline did not fire. Run resources were nice 19 and CPUs 24–27,80–83. The existing old job was untouched.

The result truthfully reports **43 components, 13 truncated components, 48 molecules, 7,734 partial paths**, `search_truncated=true`, `decisive=false`, and `cover_size_limited=true`. This verifies bounded practical completion of the pathological default-parameter case; it does **not** certify a complete mitochondrial structure. No speedup ratio against the old job is claimed because its source version differs.

Final artifacts: sibling `.agent-work/ovasm-five-fixes-final-v2-20261103/solver-default.{json,fa,log,time,exit}` and `solver-default-compact.json`, with the binary and source manifest in the same directory. The intermediate attempted build in `.agent-work/solver-fast-phase-20261103/` failed its source-stability check and was not used for this verification.
