# Native graph solver validation, 2026-11-04

The implemented change preserves read-supported alternative configurations and
reports ambiguity instead of requiring a single circular representation. It does
not establish globally optimal physical genome reconstruction. The principal
measured improvement is explanatory output coverage, not improved accuracy of
the single representative.

## Implementation

- Original-record deduplication, graph-bound strict continuous paths, and local
  adjacency evidence kept separate from complete-path evidence.
- Topology-neutral native ranking; all closing choices; bounded candidate and
  cover-size search with explicit truncation and fallback flags.
- A minimum explanatory family representing supported paths and local links;
  equivalent minimum-family candidates retained in the GFA ambiguity envelope.
- Missing supported connections and unobserved entrance/exit pairings at doubly
  branching segments prevent decisive output, even if depth missed a repeat.
- Representative FASTA, configuration-family FASTA and GFA alternatives have
  distinct meanings. Unresolved GFA paths retain their provenance for pan-graph
  consumers. Configuration files follow pipeline retry snapshot/restore rules.

Core commits: `c0e38d4`, `efd0ced`, `f4f2024`; final pipeline artifact integration:
`4095974`.

## Real within-library checks

Six pre-recruited plant mitochondrial datasets were split by original record ID
using SHA-256 modulo five. Assembly used training records only. Legacy and native
solvers received the same resulting graph; the legacy ablation removed only the
strict-record metadata. No sample-specific thresholds were adjusted. The two
prospective input hashes were frozen before their sequences/results were examined;
the tested V4 model and parameter manifest were fixed before running them.

| Dataset | Legacy representative | Native representative | Native explanatory family |
| --- | ---: | ---: | ---: |
| SRR31693266 | 58/109 | 58/109 | 109/109 (100%) |
| SRR31693270 | 44/77 | 44/77 | 77/77 (100%) |
| Buxus sempervirens | 115/180 | 105/180 | 179/180 (99.44%) |
| Adonis annua | 18/43 | 18/43 | 43/43 (100%) |
| Alnus glutinosa | 92/98 | 92/98 | 98/98 (100%) |
| Erodium maritimum | 111/167 | 111/167 | 167/167 (100%) |

Numerators count explained eligible held-out strict full-record paths. Denominators
exclude records that do not form an unambiguous continuous path across at least
three graph segments. These values are neither whole-library mapping recall nor
biological structure accuracy. Recruitment was already performed on the source
libraries; recruitment itself was not independently held out.

The native representative is unchanged on five datasets and worse on Buxus.
Reporting family recall as if it were improved single-genome reconstruction would
be incorrect. Buxus retains a cover-size limit and one unexplained held-out path;
the identical training path has one witness, below the unchanged minimum support
of two. That observation remains reported as low support. The threshold was not
lowered to improve the held-out result.

All real checks remain non-decisive. No physical topology, chromosome count or
isoform abundance is identified by these tests. Animal, fungal and microbial
performance is not established by this plant-only cohort.

## Independent optimum checks

A separate Python implementation reconstructed constraints from graph walks and
brute-forced candidate subsets. It completed on all ten cases: four clean controls
and six real datasets. Native minimum cardinalities, equivalent-family counts and
candidate membership matched its exhaustive results. For example, the Alnus
component has 59 candidates, minimum family size three and 216 equivalent minimum
families; selecting one arbitrary minimum family would hide uncertainty.

This proves minimum explanatory cardinality only within the enumerated candidates
and representable supported constraints. It cannot prove that candidate enumeration
is biologically complete. Buxus explicitly retains its search-size limitation.

## Synthetic controls and retained failure

The four clean-sequence controls recover 100% of truth canonical 31-mers, including
graph-link boundaries. All representative and configuration FASTAs spell valid
source-graph walks.

- Single circle: zero eligible full-chain held-out records, so no full-chain recall
  is measurable.
- Two genomes sharing a repeat: 29/29 eligible held-out paths represented.
- A 12 kb repeat with 6.5 kb reads: no spanning evidence; result remains
  non-decisive, all four supported local connections represented, and both
  equivalent minimum families retained. The V3 false-decisive result exposed by
  this control was fixed before V4 acceptance.
- A mixture of original and recombinant configurations: representative 39/55,
  explanatory family 55/55; result remains non-decisive.

The original mixed/noisy end-to-end stress test FAILED and remains failed. With
0.1% independent substitution errors, the upstream graph recovered 88.6648% of truth
31-mers and supplied no eligible full-chain observations. The original full-recall
assertion and failed artifacts are retained; its absence of eligible paths is not
a successful structure result.

On a separate fixed graph assembled from clean training records, the same archived
noisy records gave representative 21/27 and explanatory family 27/27, non-decisive.
This isolates evidence/solver robustness. It does not resolve or erase the noisy
end-to-end assembly failure.

## Formal-source and performance receipts

All tested Rust code was built in the formal server repository:
`/home/jiazc/data24t/work-live/PhytoOrga_AI_c/ovasm`. Frozen binaries have source and
binary hashes. Formal source manifests were checked before and after accepted
builds/tests and end-to-end validations. Editing copies and historic worktrees
were not used as tested source.

V4 release suite: **415 passed, 0 failed, 6 ignored**. Final V7 release suite:
**416 passed, 0 failed, 6 ignored**. Its release binary passed all four clean
end-to-end controls plus the Buxus integration replay, including the summary file
index, FASTA spelling, unique GFA path names and unresolved provenance. Only
`src/pipeline.rs` differs between the V4 model-validation manifest
and the final integration candidate; native evidence/ranking/set-cover code is
identical.

The 41,747-node historical graph is an input-data regression using legacy evidence,
so it tests bounded legacy-mode performance rather than strict-record ranking.
V3 completed in 24.77 seconds with 214,436 KiB peak RSS, 43 components, 48
representative molecules and 7,734 partial paths. It retained `search_truncated`,
`cover_size_limited` and `decisive=false`. Final V7 repeated the result in
23.82 seconds (wrapper elapsed; timed child 23.68 seconds), with 214,568 KiB peak
RSS, the same 43 components/48 representatives/7,734 partial paths, and unchanged
uncertainty flags. This is a regression receipt, not a performance improvement
claim against the older differently timed run.

Final binary: `.agent-work/graph-solver-optimal-20261104/ovasm-model-v7`;
SHA-256 `8a35cf619a4deda3a1c38da47b0ae72380912a97063c76582f808ec5a3ab477b`.
`model-v7-success.json` records stable sources and all final validation stages.

Receipts, commands, source manifests, frozen binaries and unabridged failed logs:
`.agent-work/graph-solver-optimal-20261104/`. Main accepted model cohort:
`e2e-v4-positive/`; exhaustive checks: `independent-v4-positive/`; failed noisy
assembly: `e2e-v4/synthetic-mixed-noisy/`; fixed-graph noise isolation:
`fixed-graph-noisy-v4/`. Final integration artifacts: `e2e-v7/`, `large-v7/`.

Windows/browser execution was not validated. No additional physical-topology truth
set or cross-taxon performance claim is supplied by this work.
