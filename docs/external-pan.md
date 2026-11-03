# External pangenome graph import

`pan import` converts an existing graph to the strict GFA1 subset consumed by reconstruction. The CLI connects to `pan_external::Config { input, extra_paths, out, sample }`. Output is a new directory containing `normalized.gfa`, `import.json`, and `path-names.tsv`; existing outputs are never overwritten.

Supported inputs: GFA1 S/L/P/W with complete IUPAC DNA sequences, explicit 0M links, and existing links for every oriented path adjacency. Gzip files ending in `.gz` are accepted. Nonzero or unknown link overlaps, gap/jump records, missing sequences, duplicate paths, incompatible W lengths, and absent path links fail explicitly rather than silently discarding information. P path names already valid PanSN retain identity; others are encoded with an original-to-normalized mapping. P/S/L optional tags, including uncertainty and reference-fill provenance, are retained. W walks remain independently named intervals; their coordinates are never interpreted as proof of a complete chromosome.

For minigraph rGFA without P/W records, rank-0 (`SR:i:0`) nodes are grouped by `SN:Z`, ordered by `SO:i`, and joined only if coordinates are contiguous from zero and each edge exists in the graph. Such paths are tagged `st:Z:reference_backbone_only`; novel nodes remain in the graph. The report explicitly marks reference-backbone-only coverage and does not claim all samples have paths. Additional P/W records can be supplied in a separate file; partial GAF alignments are not accepted or converted into whole-sample paths.

Import validation preserves segment sequences and oriented path steps, checks every link and the summed path spelling length, and checks W coordinate lengths. This is structural validity, not biological correctness or read-resolved topology. PGGB/minigraph execution and benchmarking are separate from importer tests.

Primary specifications consulted:
- https://github.com/pangenome/pggb (FASTA/PanSN input and final GFA output)
- https://github.com/lh3/minigraph (rGFA graph output)
- https://gfa-spec.github.io/GFA-spec/GFA1.html (S/L/P/W semantics)


## Build with an installed backend

```
ovasm pan build --backend pggb --assembly A=A.fa --assembly B=B.fa --sample panel --threads 8 --out built-pggb
ovasm pan build --backend minigraph --assembly A=A.fa --assembly B=B.fa --sample panel --threads 8 --out built-minigraph
```

Each labeled FASTA is one haploid sample assembly. The first assembly is the minigraph backbone.
Use `--executable /path/to/backend` for a dedicated environment; sibling executables are resolved
within that child process without changing global PATH. Input sequence names are mapped to PanSN
in `input-names.tsv`. PGGB's FASTA index is generated without an extra samtools dependency.
The original backend output, command, version, exit status and logs are retained. Failure or
ambiguous final graph selection is an error, never a successful empty reference.

Build results are ready for reconstruction at `built-pggb/normalized/normalized.gfa` (the same
layout for minigraph). Import-only usage is:

```
ovasm pan import --input existing.gfa.gz --sample panel --out imported
ovasm pan reconstruct --reference mitochondrion=mt/normalized.gfa --reference plastid=pt/normalized.gfa --reads R1.fastq.gz --reads R2.fastq.gz --sample individual --out individual
```

Build or import MT and PT references separately. Reconstruction requires both for competitive
recruitment. A minigraph rank-zero path is usable as a reference backbone, but it does not supply
all sample haplotypes. W intervals and reference-backbone-only paths are excluded from complete
sample-configuration priors; source-filled segments also make their traversing paths ineligible.
An external graph lacking a circular edge does not become circular simply because its original
FASTA was described as a circle. External assembly engines remain optional dependencies.
