# External pan graph builders

`pan build --backend pggb|minigraph --assembly SAMPLE=FASTA ... --out DIR`
accepts assemblies, not raw WGS. Each input file represents one haploid sample assembly.
The first input is the minigraph reference backbone. Sample labels must be unique and
use letters, digits, underscore, dot or hyphen. Record IDs are generated as
`SAMPLE#1#contig_N`; original headers are retained in `input-names.tsv`.

An optional `--executable PATH` selects an existing environment. Its sibling directory
is added only to the backend child process PATH so PGGB helpers can resolve. No shell
is used for dispatch and no dependency is installed automatically. Existing output
directories are refused. The prepared FASTA uses fixed 80-base lines and an FAI index.

PGGB is called with `-i all.fa -o external -t THREADS -n SAMPLE_COUNT`;
minigraph with `-cxggs -tTHREADS reference.fa sample1.fa ...`.
Raw outputs, exact arguments, backend version and logs are retained. A single
`.smooth.final.gfa` is preferred for PGGB; multiple ambiguous graphs cause an error.
The import adapter writes `normalized/normalized.gfa` and its interpretation report.
Minigraph ordinarily lacks non-reference sample walks: inferred rank-zero backbone
paths must not be presented as all sample configurations. External builders do not
establish biological truth or short-read completeness.

Official command contracts:
- https://github.com/pangenome/pggb
- https://github.com/lh3/minigraph

2026-11-03 server validation: PGGB 0.7.4 on deterministic A=150,000 bp and B=152,000 bp
(5 kb insertion, 3 kb deletion and SNPs) completed through the actual build CLI and
normalizer. Both normalized P paths spelled their complete inputs exactly, 223 nodes.
This is synthetic integration validation, not a biological benchmark.
Minigraph command dispatch/error paths were unit-tested with a fake executable;
real minigraph was not available in checked tool/environment locations. Official
GitHub clone/codeload/API downloads failed or stalled; no real minigraph execution
is claimed. Task receipts: `.agent-work/pan-backend-20261103/`.

FASTA evidence annotations are retained across build: `reads_decide=false`,
`chosen_by=fewest_circular_molecules`, `unplaced=true`, and `length_uncertain=true`
mark the corresponding normalized P path `st:Z:unresolved_configuration`.
Nonzero `reference_fill_bp` marks `rf:Z:input_reference_fill`. The input-to-path
provenance map is also saved. If an external backend loses exact path correspondence,
all emitted paths conservatively inherit the uncertainty; conflicting existing status
annotations cause an error rather than overwriting evidence or silently trusting it.
