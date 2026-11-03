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
