# Assembly assessment (`assessment.json`)

`ovasm run` ends every organelle with `ovasm check assess`, which writes `assessment.json` and a
summary of it under `"assessment"` in `summary.json`. Everything in it comes from files the run
already wrote (the assembly graph, `evidence.json`, `linearization.json`, `chromosomes.json`,
`molecules.fasta`) and the built-in protein panel. No reads are re-read and no other program is
called, so it costs well under a second. `ovasm check assess` also runs on its own on earlier
results.

| Field | What it says | How it is derived |
|---|---|---|
| `depth.single_copy_depth` | read depth of one copy of the genome | length-weighted median depth of the segments the molecules walk once |
| `copy_number` | whether the molecules walk each segment as often as its depth says | segment depth / single-copy depth against the number of walks, within 0.5; `agreement` is the share of walked bases that agree, `disagreements` the segments that do not |
| `shallow` | sequence far below one copy (NUMT, contamination, errors) | segments under 0.25 of the single-copy depth |
| `coverage_of_graph` | graph bases the molecules leave out | segments no molecule walks |
| `junctions` | reads behind every link the molecules use | link reads (and short-read adjacency reads) of `evidence.json`, both strands; `without_reads`, `min_reads`, `median_reads`, the five weakest |
| `complexity` | graph shape, as assembly graphs are classed in comparative studies, and what the reads add | `simple` (one segment), `complex` (more), `multipartite` (more than one connected component); segments, links, components, dead ends, double forks (two or more links at both ends); repeat calls of `chromosomes.json` (`interconvertible`: reads for both pairings, i.e. recombining, with the pairing counts; `not_crossed`: no read spans the repeat) and the chromosome range |
| `genes` | panel genes on each molecule and its depth | the target organelle's and the other organelle's genes; `call` is the target when it holds at least 5 target genes (plant mitochondria hold dozens of plastid genes), the other organelle only with 5 or more of its genes and twice the target's, else `few_genes`; `depth_copies` is the molecule's depth in single-copy units |
| `plastid_structure` | LSC, SSC and inverted repeat of the longest circular plastid molecule | steps whose segment the molecule walks once in each orientation are the repeat; two equal repeat runs separating two single-copy runs are `quadripartite` |

What it cannot say: per-base errors and coverage within a segment, which need the reads mapped
back (OV's `ov.qc.assembly` does that with minimap2, samtools and meryl).

Examples (HiFi unless noted):

| Assembly | Copy-number agreement | Junctions (min reads) | Complexity | Notes |
|---|---|---|---|---|
| Carex echinata plastid | 1.000 | 4 (170) | complex, 1 double fork | quadripartite: LSC 102,453, SSC 9,949, IR 33,342 |
| Lophophora williamsii plastid | 1.000 | 4 (125) | complex | quadripartite with a 1,334 bp IR (cactus plastids have lost most of theirs) |
| Carex echinata mitochondrion | 1.000 | 28 (167) | complex, 7 double forks, 6 recombining repeats | one molecule, 33 mitochondrial genes |
| Lophophora williamsii mitochondrion | 0.982 | 24 (53) | complex, 5 double forks, 3 recombining repeats | three repeats at 1.4-1.7 copies by depth: substoichiometric conformations |
| WL354 mitochondrion, short reads | 1.000 | 4 (47) | multipartite, 83 components, 170 dead ends | fragmented: 49 pieces with few genes each |
