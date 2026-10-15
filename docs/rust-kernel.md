# Graph implementations and source attribution

OVASM has its own Rust sparse de Bruijn graph in `src/assemble.rs`. The normal
`ovasm run` pipeline and `ovasm assemble` use it. The dense short-read graph is
also Rust. This does not mean every optional backend in a Linux binary is Rust.

| Entry | Implementation | Coordinates | Scope |
|---|---|---|---|
| `run`; `assemble` without a backend flag | OVASM Rust (`src/assemble.rs`) | Bases | Existing sparse graph, cleaning, gap rescue and read-supported junction restoration; short reads use dense mode |
| `assemble --multichrom-graph` | Attributed Syncasm C from Oatk 1.0, with documented OVASM patches | Homopolymer runs, restored to bases | Run-coordinate correction, repeat unzipping and shared-coordinate consensus; Linux only |
| `assemble --rust-hpc` | OVASM Rust (`src/hpc.rs` + existing Rust graph) | Homopolymer runs, restored to bases | Experimental graph and run consensus; **no read correction or repeat unzipping** |

`--rust-hpc` and the Cargo feature selection described here are source development
after the `v0.2.0` tag. They are **not present in the v0.2.0 precompiled packages**.
Existing defaults and the published tag are unchanged.

## Method versus implementation

Sparse de Bruijn graphs and syncmer sampling are scientific methods, not an
exclusive algorithm invented by Syncasm. Syncmers were introduced by Robert
Edgar ([PeerJ 2021](https://peerj.com/articles/10805/)). Oatk provides a specific
assembler, Syncasm, combining sampling, correction, graph cleaning and
disentangling ([Oatk paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12329965/)).

The embedded C implementation really is third-party code, not an independent
Rust implementation. Its MIT license, Chenxi Zhou attribution and OVASM patches
remain in `vendor/syncasm`. Calling it through Rust FFI does not change its origin.
The new Rust HPC backend reuses OVASM's existing Rust graph implementation and
does not call or translate the C kernel.

## Independent Rust HPC development

The Rust backend compresses each original ACGT homopolymer to one run. Ambiguous
bases break sampling. It constructs the graph using the existing closed-syncmer
sampler, counts read adjacencies and compacts unitigs. During this stage, unitig
overlaps are retained; base-coordinate gap filling and boundary trimming are
disabled until run lengths have been restored.

Exact, oriented HPC overlaps identify shared graph positions using a union-find.
Exact sampled read words anchor each original read to the graph; hash hits must
pass a full sequence comparison. Collinear anchors extend only through exact
sequence. Each distinct `(input record, compressed read position, shared graph
position)` contributes once. Terminal runs are excluded because the read can end
inside a homopolymer. Multiple repeat placements remain multiple placements;
they do not establish locus, organelle origin or transfer direction.

The base length is the rounded mean of these original run observations, pooled
on shared graph coordinates. Unsupported positions fail without emitting a base
graph. Restored overlaps must match exactly in both orientations; sequence is
not clipped to conceal inconsistent overlaps. Intermediates and input statistics
are retained beside the requested graph. Reports explicitly mark the backend as
experimental, `uses_syncasm: false`, `read_error_correction: false`,
`read_threaded_unzipping: false` and `physical_topology_identified: false`.

This is not yet a validated replacement for the corrected C multichrom backend.
It still inherits the Rust graph's existing filtering and gap-rescue rules.
Those rules have not been newly certified as safe for every chromosome, rare
branch or nuclear-context alternative. Neither a graph self-loop nor recovery of
synthetic circles establishes physical circular chromosomes in a real sample.

Build without compiling or linking the Syncasm C kernel:

```sh
cargo build --release --no-default-features --locked
```

For development testing, use already recruited HiFi reads and an explicit,
read-derived cutoff. No reference or expected molecule count is used to build
the graph:

```sh
ovasm assemble --rust-hpc --reads recruited.hifi.fastq \
  --k 1001 --s 31 --min-count 30 --threads 8 \
  --out candidate.gfa --out-json candidate.json
```

`30` is an example argument, not an automatic or universally recommended cutoff.
The experimental command requires a new output path and does not silently reuse
intermediates from a failed run. It is not automatically selected by `run` or by
`--multichrom-graph`.

The finite real-input development check is reproducible through
`scripts/assembly_benchmarks/validate_rust_hpc.py`. It compares the existing ZM4
recruited inputs at their saved k/s/cutoff, records original input counts,
raw-read link evidence, representative sequences, wall time and GNU peak RSS.
The earlier OVASM outputs are development controls, **not independent truth**.
Report outcomes, including failures and limitations, before considering a default
backend change. Corrected read walks, repeat phasing, source-context decisions,
minor branches and independently evaluated chromosomes remain separate gates.

### Existing ZM4 development input (15 October 2026)

The C-disabled build completed both recruited-input cases at their saved settings:

| Input | Complete records / bases | k / s / cutoff | Graph segments / links | Primary records | Graph wall time | GNU maximum RSS |
|---|---|---|---|---|---|---|
| MT | 29,287 / 516,099,590 | 1001 / 31 / 53 | 5 / 8 | 307,026 bp | 30.65 s | 371.85 MiB |
| PT | 95,680 / 1,674,746,317 | 1001 / 31 / 323 | 1 / 1 | 125,509 bp | 95.94 s | 576.09 MiB |

These times include compression, Rust graph construction and run-length
restoration, using eight threads on a shared server. Recruitment is excluded;
these are not timings for scanning the original 7.35-million-read WGS input.
Graph segment lengths include overlap sequence, while primary lengths exclude
walk overlaps. They must not be summed or compared as though they were the same
quantity. Peak RSS is GNU `time -v` for the assembly command, not a process-tree
RSS sum or the maximum for the earlier full recruitment pipeline.

All compressed graph positions had uncensored original observations, all
restored-base overlaps matched, and subsequent evidence scans found raw-read
support for 8/8 MT and 1/1 PT links. This is local junction support, not proof that
every complete path is supported or that the graph is biologically complete.

The PT sequence matches the previous 125,509-bp representative allowing reverse
complement and rotation. MT exports a 307,026-bp representative instead of the
previous two 160,953/146,073-bp representatives; those two records do not match the
new record individually. This comparison therefore **does not demonstrate
equivalence of the MT path representation**. No independent truth, physical
chromosome count or certified organelle origin is inferred from this check.

The initial shared-overlap unit test exposed unobserved end positions: exact
sampled anchors did not necessarily reach an arbitrary segment boundary. The
fix extends a matched diagonal through exact sequence to its actual boundary;
it does not assign a default run length to unknown positions. The original
failing logs are retained alongside the passing rerun. Five new HPC/CLI checks
and 16 existing Rust assembly checks passed on each of Windows and Linux; these
small checks and the two development inputs are separate evidence scopes.

As a separate default-route regression, the new C-disabled executable and the
published v0.2.0 executable assembled both saved MT/PT inputs without `--rust-hpc`.
Both output GFAs were byte-identical and all the science fields enumerated in
`tests/rust_hpc_default_regression.py` agreed. This checks preservation of the
existing default on these two inputs; it is not a full regression panel or a
proof of equivalence between the experimental HPC and default backends.
