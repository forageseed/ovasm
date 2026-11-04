# Read-backed graph configurations

New `check evidence` reports include graph-bound, deduplicated original-record
evidence. `solve linearize` and the native pipeline select the new model
automatically when that metadata is present. Older evidence reports retain the
legacy bridge/copy model. A graph fingerprint mismatch is an error: evidence
cannot be silently reused after changing sequences, links or overlaps.

## Evidence and ranking

An original record contributes at most once to each canonical adjacency and once
to its complete accepted graph walk. Identical duplicate records do not increase
support. Conflicting informative long-read sequences under one record ID fail.
Short-read mates retain distinct sequence identities. Disconnected chains,
unresolved assignments and non-ACGT bases inside a proposed spanning interval do
not become a single phased walk; their valid local adjacencies remain available.

The native objective is lexicographic: supported joint paths spanning at least
two complete repeat visits, supported complete record paths, then supported local
adjacency occurrences, followed by junction and copy/representation tie breaks.
Circularity has no evidence bonus. The adjacency count is correlated: one record
can support multiple adjacencies, so this is not an independent sample size, a
likelihood, an abundance estimate or a posterior probability. `margin_ln` is null
for this model.

All admissible closing bridges and open alternatives enter the bounded search.
Small components search all cover sizes automatically; larger components retain
explicit limits. Exhausting the read-scoring budget discards every partial native
score and falls back to a labelled, non-decisive legacy ranking.

A candidate cannot be decisive if it omits supported local connections or accepted
full-record paths. At a segment with multiple supported entrances and exits,
separate entrance and exit coverage cannot identify their pairing: an original
record must span that pairing with sufficient support. This check also applies
when the earlier depth heuristic has not labelled the segment a repeat.

## Three outputs with different meanings

| Output | Meaning |
| --- | --- |
| `molecules.fasta` | One ranked representative configuration. |
| `molecules.configurations.fasta` | One compact explanatory family of configurations representing supported record paths and local adjacencies. Headers identify configuration membership. |
| `organelle.gfa` | Original graph plus representative and alternative configuration paths, including the union of equivalent minimum-family candidates when enumerable. |

The explanatory family is a set-cover representation of observed constraints,
not an inferred number of physical molecules, isoforms or chromosomes. It may
contain overlapping or incompatible configurations. Do not concatenate records
across configurations into a proposed genome.

The family is chosen from the full enumerated candidate list before the ordinary
alternative display cap. A bounded exact search can prove its minimum cardinality
when the supported constraint universe fits the exact solver; otherwise the
valid greedy family remains available without a minimum proof. A second bounded
enumeration retains equivalent minimum families. Its completeness has a separate
flag. Unresolved exported GFA paths carry
`st:Z:unresolved_configuration`, so downstream pan-graph code does not promote
them to trusted priors.

## Reading certainty and limits

Inspect `linearization.json`, especially these fields:

- `ranking_model`, `read_model_budget_exhausted`: objective used and any fallback.
- `search_truncated`, `cover_size_limited`: candidate enumeration limits.
- `read_support.unexplained_thread_records`: accepted full-record evidence absent
  from that candidate, including below-threshold observations.
- `read_support.unexplained_adjacency_record_occurrences` and
  `unphased_branch_visits`: missing supported connections and unobserved pairings.
- `read_path_cover`: supported/represented/unrepresentable record counts and local
  connection counts; low-support records are reported separately.
- `minimum_family_size_proven`: optimum only within the enumerated candidates and
  supported, representable constraints.
- `equivalent_families_complete`: whether equivalent minimum families were fully
  enumerated within that same finite candidate space.

`decisive` describes the bounded ranking result under these criteria.
`physical_topology_identified` remains false. Even complete minimum-family
enumeration cannot establish physical topology from an observationally ambiguous
graph, recover a missing source-graph connection, or prove that an upstream
depth/copy classification is correct.

## Reproducing validation

`scripts/validate_native_read_solver.py` requires a binary built from the formal
source directory and a manifest of the exact Rust sources and Cargo files. It
verifies source stability, hashes the binary, analysis scripts and input reads,
and checks FASTA spelling against graph paths. Synthetic truth never enters the
assembler or solver. Real checks split original record IDs into training and
held-out sets, build the graph from training only, and compare the native model
with its legacy ablation on the same graph.

`scripts/verify_native_cover.py` independently brute-forces candidate subsets from
a full report (`solve linearize --max-alternatives 10000`), reconstructing path
and adjacency constraints from evidence. It compares minimum family cardinality
and equivalent-family membership without importing the native set-cover code.
It reports its own combination limit.

Within-library held-out linkage is a predictive check, not independent biological
topology truth. Synthetic sequence recall, primary linkage recall, family linkage
recall and ambiguity-envelope recall must be reported separately. A simple graph
with zero eligible full-chain held-out records has no measurable full-chain recall.
These checks do not establish animal/fungal performance or solve read recruitment
across distant taxa; those require separate data and controls.

The separate `synthetic-mixed-noisy` end-to-end assembly stress test is retained
with its original full-sequence-recall assertion. The November 4 development run
failed that assertion: the upstream graph retained 88.66% of truth 31-mers and had
no eligible three-node held-out paths. This is an unresolved assembly limitation,
not a successful graph-resolution test. Explicit `--synthetic-case` arguments can
select the four clean-sequence graph controls for the positive validation suite;
the selected protocol is recorded in provenance. Testing the same noisy records
on a graph built from separate clean training reads measures evidence/solver
robustness only and cannot turn the noisy end-to-end failure into a pass.

For example, from the formal repository after verifying the recorded source
manifest still matches:

```sh
python3 -I scripts/validate_native_read_solver.py \
  --binary .agent-work/graph-solver-optimal-20261104/ovasm-model-v7 \
  --source-manifest .agent-work/graph-solver-optimal-20261104/model-v7-before.json \
  --out .agent-work/graph-solver-optimal-20261104/recheck-clean-controls \
  --synthetic-case synthetic-single \
  --synthetic-case synthetic-shared-repeat \
  --synthetic-case synthetic-long-repeat \
  --synthetic-case synthetic-mixed
```

Add `--sample NAME=/absolute/recruited.fastq` for real within-library checks. Select
`--synthetic-case synthetic-mixed-noisy` separately to reproduce the retained
upstream-assembly stress failure; a nonzero exit must not be reported as success.
