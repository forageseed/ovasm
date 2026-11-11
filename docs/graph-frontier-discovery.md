# Discover graph-anchored raw sequence differences

`scripts/discover_graph_frontier.py` nominates original-read intervals containing exact words that the fixed graph cannot spell. The report uses graph and read coordinates and carries no biological source labels.

## Inputs and output

Use a JSON manifest with schema `ovasm.graph-frontier-discovery.v1`, a nonempty `scope`, a `graph` descriptor and a nonempty `raw_reads` array. Each descriptor supplies `path` and lowercase `sha256`. Paths resolve relative to the manifest. Provide either `recruit`, the descriptor of an OVASM `find recruit` report, or an explicit `scan_k` from 1 to 31. The scan length refers to the recruitment word length. It is independent of the assembly k.

~~~bash
python3 -B -I scripts/discover_graph_frontier.py manifest.json \
  --out discovery.json --out-evidence raw-evidence.jsonl.gz
~~~

The graph must use GFA1 inline nucleotide sequences and links with exact `nM` overlaps. Each overlap must agree with both oriented segment sequences. The scanner includes reciprocal links and traverses paths across any number of short segments. Graph headers and path annotations are ignored. Unknown graph records, sequence-less segments, and non-exact overlaps are rejected.

The raw files can contain FASTA or four-line FASTQ, optionally gzip compressed. Evidence records include the input and record ordinals, the full original header encoded as hex, the SHA-256 of the literal sequence and the SHA-256 of FASTQ quality when available. Coordinates are zero-based and half-open. Graph anchors include both oriented and forward-segment offsets.

Sequence matching is case-insensitive. FASTA sequence lines are trimmed and joined, preserving the case of the resulting sequence bytes for the digest. The report separates records with no represented word among observable positions from records with no observable word at all. It also counts positions containing any N or other non-ACGT symbol as unobservable. Unobservable positions break candidate runs. `--include-no-graph-words` emits individual rows for both types of record. Counts are FASTX records. Paired read ends, repeated exports and independent original molecules require separate identity handling.

## Interpretation and bounds

The scanner treats an exact word as a unique node anchor only when it appears once in the oriented node sequences and never across a graph edge. The cross-edge rule is conservative for words with multiple overlap representations. A candidate is a run of missing raw words with at least one directly adjacent unique node anchor. The affected raw interval is the union of those words and can overlap an anchor. The report records two-sided anchoring, distinct graph trajectories and terminal or branch geometry separately.

The graph traversal budget includes reached continuation states and defaults to 1,000,000. Set a different budget with `--max-traversal-states`. If the graph index is truncated, the report still binds all inputs but leaves raw counts unavailable and writes no evidence rows. A partial index cannot establish word absence or anchor uniqueness.

A complete scan finds differences relative to the supplied graph. These differences can arise from variants, sequencing errors, incomplete assembly or sequences from another source. The output does not assign origin or authorize recruitment, and reads disconnected from the graph cannot become graph-anchored candidates.
