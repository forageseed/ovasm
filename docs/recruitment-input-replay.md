# Recruitment input replay

`ovasm run` and `ovasm recruit` accept `--recruit-cache NEW_DIRECTORY`.
This is an opt-in execution optimization. The scratch directory must not exist
and must have space for the original sequences and identifiers. Uppercase ACGT
is stored in two bits per base; other sequences are stored verbatim. No records
are sampled, shortened or dropped, and the original batch boundaries are retained.

The first recruitment round reads the original input and caches its sequences.
Later rounds read the cache. Each bait mutation records every canonical key whose
target mask changed, including an existing key gaining a competing target bit.
A read without any such key has precisely the same hits, compatible blocks,
scores and ties as in the preceding round, so its classification is reused.
Every other read is classified by the original classifier. Biological iterations,
extension, stopping, depth gating and deterministic subsampling are unchanged.

Intermediate private pass files need only sequence and identifier. After the
last biological round, one final original-input pass supplies the actual quality
strings. Identifiers and sequences are compared directly, record by record,
against the cache; input metadata and final record/base counts are checked.
For a multi-round run this reduces original-input reads from N classification rounds
to **two** input passes, not one. A one-round run already has original qualities
and omits the extraction pass. Cache passes still traverse complete sequences to test changed
keys. Optional bootstrap/depth-gate scans are separate and are not removed by
this optimization.

`recruit.json.input_replay` reports the source classification pass, replay passes,
cache bytes, final quality extraction pass and number of reads actually
reclassified in each biological round. Existing `reads_scanned` counts continue
to mean all logical records considered in that round. On successful completion
the invocation removes only its own scratch file/directory; a failed cache is
retained for inspection and never reused automatically.

Large-input speed and memory must be measured on the complete original input.
Toy equivalence checks or prefixes do not establish the full-run speedup. Compare
final target FASTQ and merged candidate FASTQ byte for byte, the two organelle
FASTA files, biological recruitment statistics and scientific reports; execution
timings, paths, new evidence exports and replay profiling are reported separately.

```sh
ovasm run --reads full-hifi.fastq.gz --organelle both --read-type hifi \
  --read-set whole-genome --recruit seeds --mode fast --threads 8 \
  --out sample --recruit-cache /scratch/new-sample-recruit-cache
```
