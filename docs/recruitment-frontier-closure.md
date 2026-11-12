# Fixed recruitment frontier replay

Run the checker on a directed graph of proposed physical connections and evidence prepared under one fixed eligibility policy.

~~~sh
python3 -B -I scripts/recruit_frontier_closure.py manifest.json \
  --out frontier-report.json --round-budget 6 --edge-budget 500
~~~

The report describes closure within the declared graph and read pool. The checker verifies hashes, receipt binding, and reachability. Source, molecule independence, and statistical calibration come from the supplying evidence. This command is separate from default recruitment. An optional export writes a new bridge FASTA from verified external receipts.

## Inputs

A file reference is an object with `path` and lowercase `sha256` of the exact file bytes. Relative paths resolve against the file containing that reference. Use immutable evidence artifacts. The read-pool reference ordinarily points to an upstream manifest of original read-file hashes. The checker locks that manifest and does not fetch or realign the original reads.

| Manifest field | Value |
| --- | --- |
| `schema` | `ovasm.wp4.frontier-closure.v1` |
| `scope` | Stable local scope identifier |
| `graph` | Frozen graph reference |
| `eligibility` | Frozen policy reference |
| `read_pool` | Frozen original read-pool receipt reference |
| `scope_receipt` | Coverage attestation reference |
| `seeds` | Nonempty list of the original seed node IDs |
| `edge_evidence` | List of edge evidence references, which may be empty |

The frozen graph uses schema `ovasm.wp4.frontier-graph.v1`, a `nodes` list of identifiers, and `edges` entries with string `id`, `from`, and `to` fields. Edge IDs are unique even when endpoints coincide. Normalize biological orientation and physical edge aliases upstream so the graph binds the actual directional claims.

Both scope and edge receipts contain the following binding fields.

| Field | Value |
| --- | --- |
| `scope` | Manifest scope identifier |
| `graph_sha256` | Frozen graph file digest |
| `eligibility_sha256` | Fixed policy file digest |
| `read_pool_sha256` | Frozen read-pool receipt digest |
| `seed_set_sha256` | Digest of the JSON encoding of sorted distinct seed IDs, with no spaces and Python default ASCII escaping |

Use schema `ovasm.wp4.frontier-scope.v1` for the coverage attestation. It also contains `coverage` of `complete` or `partial`, a nonempty `explanation`, and a nonempty `provenance` list of file references. A complete attestation means the candidate graph covers the stated local domain. A partial attestation prevents convergence even if the declared graph reaches a fixed point.

Use schema `ovasm.wp4.frontier-edge.v1` for each edge receipt. Include the binding fields, `edge_id`, `decision`, `basis`, `full_edge`, `opportunity_fragment_ids`, `support_fragment_ids`, a nonempty `explanation`, and a nonempty `provenance` list. Use upstream identifiers for independent original fragments. Both lists are deduplicated. Supporting fragments must occur in the opportunity set. Support refers to the stated basis, for example the competing source for `competing_context`.

| Decision | Basis | Additional requirements |
| --- | --- | --- |
| `admitted` | `target_exact_bridge` | A full-edge certificate, nonempty opportunity and support sets, and external target-context evidence |
| `excluded` | `competing_context` or `calibrated_absence` | A full-edge certificate and nonempty opportunity set. Competing-context exclusion also requires supporting fragments. |
| `unresolved` | `opportunity_only`, `no_opportunity`, `unmeasured`, or `ambiguous` | Opportunity-only receipts identify fragments. No-opportunity and unmeasured receipts leave opportunity and support sets empty. |

An absent receipt leaves that edge unresolved if its origin becomes reachable. Distinct receipts with conflicting decisions also leave it unresolved. Copying a receipt or repeating a fragment ID contributes no further decision. A stale hash or incompatible binding rejects the manifest.

## Reading the result

Each complete round evaluates all currently reachable frontier edges together and adds only destinations of admitted edges. Cycles and duplicate inputs do not change the final reachability. If the remaining edge budget cannot cover the next full wave, the checker stops before that wave. Both budgets may be omitted. The budgets govern graph evaluation after verifying supplied evidence.

- `converged` means the scoped reachable frontier was exhausted, the scope attestation was complete, and no evaluated edge remained unresolved.
- `unresolved` means the frontier was exhausted for the supplied admissions, with a held edge or incomplete scope.
- `truncated` preserves evaluated decisions and names the next unevaluated frontier. An edge beyond that frontier is `not_reached` in the current result.

The output records seed and reachable sets, separate admitted and held edges, round transitions, and per-edge evidence with counts and provenance. The decision digest excludes file paths and receipt order. It is a reproducibility aid alongside the exact input file digests.

## Exporting admitted bridge seeds

Add `--out-seeds admitted.fa` to export the exact sequences of the edges the replay has evaluated and admitted.

~~~sh
python3 -B -I scripts/recruit_frontier_closure.py manifest.json \
  --out frontier-report.json --out-seeds admitted.fa --round-budget 6
~~~

Each certificate for an exported edge must supply a `bridge_sequence` object:

~~~json
{
  "bridge_sequence": {
    "path": "exact-bridges.fa",
    "sha256": "<SHA-256 of the complete FASTA file>",
    "record": "bridge-17",
    "sequence_sha256": "<SHA-256 of the normalized record sequence>"
  }
}
~~~

The path is relative to the edge receipt. The record is the first token of a FASTA header and must occur exactly once in that file. Normalize the selected sequence by removing whitespace and uppercasing ASCII bases, then hash the resulting A, C, G and T bytes without a header or newline. The same edge's admitted certificates must agree on the oriented sequence. Identical normalized sequences from multiple edges appear once in the FASTA.

The report's `seed_export` gives the output FASTA digest, the record-to-edge mapping, and the source file and receipt for each record. `budget_truncated` indicates that a later wave was left unevaluated. The FASTA contains the verified admissions already evaluated in that budget. It is empty when there are no admissions. Held, excluded, unevaluated, and unreachable edges supply no sequences. The replay's ordinary decision digest is independent of this opt-in export.

Both outputs must be new, distinct paths. The checker validates the complete export before publishing either file and places the report last as the receipt for the FASTA. A handled publication failure rolls back the links installed by that invocation. The report preserves the replay's scope and status. The exported sequences come from external certificates, and the checker does not rerun their read or source analyses. Adding the FASTA to a k-mer bait can recruit ambiguous or competing reads that share its k-mers. Use the FASTA as a diagnostic input. It does not automatically admit reads or establish their physical source.
