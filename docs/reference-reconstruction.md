# Short-read reconstruction from species reference graphs

`reconstruct` accepts mitochondrial and plastid OV-GFA graphs produced by `ovasm pan`.
It uses every source configuration path for recruitment and for read-to-graph projection.
Both reference graphs are required so that shared organelle sequence competes under the
existing Illumina recruitment rules.

```bash
ovasm reconstruct \
  --reference mitochondrion=species-mt.gfa \
  --reference plastid=species-pt.gfa \
  --reads individual_R1.fastq.gz --reads individual_R2.fastq.gz \
  --sample individual --threads 8 --target-depth 150 --out individual-reconstruction
```

The optional depth argument controls resource use. Its denominator is the total stored
node sequence of each reference graph, with each node counted once. This is explicitly
reported as a sampling basis, not physical genome size. Summing the lengths of all
source genomes and alternative configurations would inflate that denominator. Omit
the argument to retain every recruited read. The existing recruitment classification,
extension, saturation and short-read de novo assembly rules remain in use. Pairing is
preserved in the input files but insert-size constraints and mate phasing are not used.

The workflow proceeds as follows:

1. Seed recruitment with all reference P walks (including circular-boundary words) and
   all graph node sequences. Both targets compete in the same recruitment call.
2. Match complete retained reads, in either orientation, exactly against every source
   P walk. A sparse 21-mer index proposes positions; full sequence comparisons verify
   them. Reads span multiple graph nodes, including one-base variant nodes. Equivalent
   placements in different source paths count once per read record and graph coordinate.
3. Retain reference intervals covered in every compatible graph placement of a read.
   Ambiguous alternatives contribute separately to possible coverage. Keep a reference
   link only if a read crosses it in every compatible placement. This is local adjacency
   evidence, never long-range repeat phase. No coverage or origin threshold is inferred.
4. Assemble the retained target reads de novo using the existing `sr` pipeline. Merge
   that graph with the supported reference intervals using the existing exact pan merge,
   retaining every input segment, link and path. This retains novel target sequence;
   unobserved reference sequence is not used to fill gaps.
5. Write `individual.gfa` and `individual.fasta` for each organelle. FASTA records are
   maximal nonbranching graph walks. Branches remain in the GFA, and a graph cycle does
   not certify a physical circular chromosome.

Important outputs for each organelle are:

| File | Meaning |
| --- | --- |
| `individual.gfa` | Integrated individual graph, including unresolved branches |
| `individual.fasta` | Nonbranching contigs from that graph |
| `observed-reference.gfa` | Only reference intervals/links supported by exact target reads |
| `reference-intervals.tsv` | Every reference interval: certain exact, ambiguous only or unsupported |
| `node-coverage.tsv` | Certain and possible coverage per reference node |
| `edge-support.tsv` | Record counts crossing reference links and whether each link was emitted |
| `projection.json` | Actual scanned records/bases, exact/ambiguous/unmatched counts and limitations |
| `de-novo/` | Independent target de novo graph, FASTA, evidence and acceptance report |
| `integration.json` | Preservation checks for the two input graphs; input provenance is not sample count |

`reconstruction.json` reports both organelles and actual recruitment counts. A de novo
error is retained explicitly even if a reference projection is available. `complete_assembly`
remains false: local read support alone does not identify complete chromosome paths.
Strict and relaxed scores stay null until an independent compatible truth is supplied.
Physical molecule count stays null and physical topology stays unidentified.

An explicit target-input stage can reuse previously recruited files without scanning the
whole WGS library again:

```bash
ovasm reconstruct \
  --reference mitochondrion=species-mt.gfa --reference plastid=species-pt.gfa \
  --target-reads mitochondrion=earlier/recruit/mitochondrion.fastq \
  --target-reads plastid=earlier/recruit/plastid.fastq \
  --sample individual --threads 8 --out target-stage
```

This mode reports `explicit_target_reads` and does not run recruitment or depth
sampling. Its upstream recruitment is not certified by the command itself; retain
the original input/recruitment manifests when chaining stages. It cannot be combined
with raw `--reads` or `--target-depth`. A de novo failure writes its diagnostic report
and partial graph evidence, then exits with an error rather than declaring the whole
workflow completed.

The exact projection has a deliberately bounded interpretation. Sequencing errors and
alleles absent from all source paths fail the complete-read exact test; the de novo
graph is their separate sequence channel. A new recombinant target walk may not occur
in a source P walk. Compatible nuclear insertions and sequence shared between organelles
can also match perfectly. Neither exact coverage nor a HiFi reference certifies nuclear
versus organelle origin. Short reads cannot select a whole configuration through long
identical repeats. Those uncertainties must be retained when comparing materials or
building a species pan graph from the resulting individual graphs.

An eight-run reference graph is a usable panel, not a universal minimum sample count
or proof of eight independent biological individuals. Reference and target identities
must be checked separately; a target's held-out HiFi graph must not be added to its
reference or used to select paths during reconstruction.
