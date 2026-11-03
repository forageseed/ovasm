# Pan-assisted WGS and external builders: verification, 2026-11-04

Actual source: `/home/jiazc/data24t/work-live/PhytoOrga_AI_c/ovasm`.
All new product runs were compiled there. Old worktree directories supplied historical inputs
and comparisons only. Concurrent depth/ONT changes remained untouched; these are working-tree
builds with source/binary SHA manifests, not claims about a clean tagged release. Server CLI
verification only; Windows and browser were not verified.

## Correctness fixes

- Unique reference-prior probe IDs despite repeated flank sequences.
- Preserve P-path provenance through pan merge; propagate reference-filled S nodes to paths.
- Exclude unresolved/reference-filled/interval/backbone-only paths from complete-sample priors.
- Reject conflicting reference fills rather than selecting by global similarity.
- Independent cumulative bait-growth guard in pan recruitment; growing is non-converged.
- Frontier depth uses supported reference anchors and rejects isolated sparse support.
- Paired-read name normalization, shared anchors, circular-origin invariant repeat auditing.
- Closed graphs with unresolved repeat lengths report closed_length_uncertain, and FASTA
  headers retain length_uncertain=true. Uncertain outputs are not automatically certified
  and written back as independent reference evidence.

## Tests and genuine runs

Relevant final pan/recruit/frontier/uncertainty tests: 72 passed before external integration.
Final external/reference/evidence/closure/uncertainty test set: 69 passed; external backend
cwd regression set also passed (7). Earlier broad regression: 326 passed, 5 ignored, before
later additions. Counts are separate test runs, not summed coverage.

Final source hashes and release binary: `.agent-work/pan-evidence-fixes-20261103/backend-cwd-final/`.
A build interrupted by concurrent source edits was explicitly rejected; a stable source-hash
build was used for final verification. The root target/release binary is built from the formal
source, and immutable tested copies are retained alongside their manifests.

Direct external GFA import followed by actual target reconstruction reproduced synthetic
8,000 bp MT and 6,000 bp PT circles exactly (rotation/reverse-complement equivalent).
Evidence: `.agent-work/pan-external-direct-smoke-20261103/`.

Real PGGB 0.7.4 build/import: 223 nodes, two complete paths (150,000 and 152,000 bp), both
exactly spell the inputs including a 5 kb insertion, 3 kb deletion and SNPs. Final replay
also preserved uncertain/reference-fill FASTA provenance on normalized paths. PGGB requires
FAI; it is now generated. seqwish/SDSL attempted cwd caches in the broken home directory;
backend cwd isolation fixed this, confirmed using identical FASTA/PAF and a final actual run.
Evidence: `.agent-work/pan-backend-20261103/` and backend-cwd-final/verified.json.

Minigraph dispatch/rGFA import have fixture and direct-parser coverage, but the actual optional
executable was unavailable. Official GitHub and Bioconda retrieval failed on the server. No
mock result is counted as real minigraph execution. Existing installations can be specified
with --executable. Outputs are normalized GFA suitable as reconstruction references; rank-zero
backbones do not imply all sample paths are known.

## Biological outputs and costs (single shared-server measurements)

Final same-recruited-input ZM4 reconstruction retains MT 307,026 bp + 453 bp and PT
125,243 bp + 92 bp. MT configuration is unresolved; PT has seven length-uncertain regions.
No reference bases were forced in to reproduce the different HiFi specimen's lengths.
The initial fixed target-stage run took 165.94 s; this excludes WGS recruitment.

Full WGS fixed-binary pan reconstruction: 1,304.85 s (21m45s), MT+PT.
GetOrganelle 1.7.7.1 full-WGS PT: 1,454.27 s (24m14s), peak sampled RSS sum 4.469 GiB,
exit 0 but an explicitly incomplete 124,952 bp scaffold. It auto-sampled 5,215,048 pairs.
These different output scopes do not justify a cross-tool full-workflow speedup ratio.

Matched PT assembly-only input (856,699 identical records): OVASM 121.38 s vs SPAdes+
GetOrg extraction 160.78 s; 1.325x in this one trial, 24.5% less time. Peak sampled RSS sum
1.034 vs 1.957 GiB. Input preparation 7.45 s was reported separately. OV gives a 125,243 bp
circle + 92 bp; GetOrg an incomplete 124,952 bp scaffold. Repeats remain ambiguous in both.
A launcher-specific failed SPAdes attempt is preserved and not called an algorithm failure.

Optional 20% pair-name sampling experiment: all 77,777,446 pairs checked; 15,556,683 retained.
Sampling 290.28 s + pan reconstruction 263.83 s = 554.11 s (9m14s), versus 1,304.85 s full
input using the same frozen binary. MT main sequence and all 15,135 bp historical missing
intervals were retained, but PT called sites changed 232->231 and MT minor-allele/low-coverage
presence calls changed. This is a candidate fast assembly workflow, not a replacement for
full-input variant analysis or a universal performance guarantee. See pan-fast-assembly.md.

GetOrg comparison evidence: `/home/jiazc/data10t/ovasm-dev-work-20261103/getorg-comparison/`.
Full and target pan runs: sibling `pan-evidence-validation/`; sampling: sibling
`pan-wgs-pair20-experiment/`. All original failed attempts and alternative outputs retained.

## Remaining limits

The complete original WGS anchor scan finished: 155,554,892 records, 2,558.74 s. At the four
PT difference intervals (20-187, 412-496, 107981-108116, 108541-108800), observed gap lengths
and read-record support were respectively {53:18,110:2}, {2:7,21:13,42:8,63:4,84:6}, {61:6},
and none. This reveals some evidence missed by the retained cache, but does not determine the
dominant long allele: exact anchors miss errors, records are not independent molecules,
shorter alleles are preferentially observable in 150 bp reads, and origin is not certified.
Do not force the HiFi lengths or replace the assembled allele from these biased sparse counts.
The old JN-L growth-guard baseline recruited 1,300,083 MT reads and reached saturated without
triggering the growth guard. This outcome is not causal proof of the guard's effectiveness;
its assembly and the separate frontier-guard run must retain their own version/status records.
