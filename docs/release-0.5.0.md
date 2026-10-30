# OVASM 0.5.0

`ovasm reconstruct` turns another individual's short WGS reads and a species reference graph
into that individual's organelle genomes, as circles where the reads allow, with its genotypes
against the graph, and writes it back into the graph.

## Reconstruction (`ovasm reconstruct`)

- Frontier walk: one pass counts the library's 31-mers (count-min, `--frontier-memory-gib`,
  default 4 GiB) and walks from every place where the reference stops being present, through
  k-mers at 0.4-10x the target's depth, until it has stayed on the reference for 300 k-mers;
  the walked pieces join the bait (`--no-frontier` to skip). ZM4's mitochondrion holds 15 kb no
  reference run has; recruitment alone stopped at its edges, and extending to a fixed point ran
  into the nucleus.
- Short reads are assembled at several k, as GetOrganelle does: the read-length k, then 99, 77
  and 55 until one molecule; the fewest pieces of at least 1 kb, then the most bases, wins.
  This also applies to `ovasm run` with short reads.
- A single linear molecule is closed with read pairs: pairs across its ends give their distance,
  and the reads anchored there with their mates spell the bases between. When they cannot, a
  reference walk of agreeing length, from the run most like the individual, may fill the join,
  labelled as reference bases. Recruitment for reconstruction keeps mates together.
- When the reads leave repeat pairings open, every configuration the graph allows is written
  (`<organelle>.configurations.fasta`, P walks of `<organelle>.gfa`), ranked by how often the
  reference runs hold each pairing; the circular one with fewest molecules is the representative
  (`circular_configuration_chosen`).
- Genotypes at every simple bubble of the reference graph (`genotypes.tsv`, `presence.tsv`),
  a Li-Stephens mosaic over the reference runs that imputes uncalled sites (`mosaic.tsv`), a VCF
  on the backbone walk (`genotypes.vcf`), and the individual merged back into the graph
  (`<organelle>.pan.gfa`).
- Stretches whose length the reads do not fix (tandem arrays, microsatellites longer than a
  read) are listed in `uncertain_regions.tsv`; lengths are not corrected.
- Results at the top of the output: `<organelle>.fasta`, `<organelle>.gfa`, and a `final`
  record per organelle in `reconstruction.json` (status, reasons, closure, configurations,
  genotypes, write-back).

## Assembly

- A dead-end branch as shallow as a NUMT's reads is cleaned (f6d1745).
- `chromosomes`: a component far from the main one's depth is `source_pending` (6f8a965).

## Validation

ZM4 (`SRR13333747`, 155,554,892 records, 23.3 Gb) against the eight-run HiFi panel, end to end
with `--target-depth 1000`: 1:07:10, 4.73 GB peak (31-mer counting 2,734 s on a slow disk). The
held-out ZM4 HiFi assembly is used only to check.

| | Mitochondrion | Plastid |
|---|---|---|
| Representative | one 307,026 bp circle (configuration 3 of 32) + 453 bp unplaced | one 125,243 bp circle |
| Against the HiFi assembly | both HiFi molecules joined through the 8.3 kb repeat, every base identical; configuration 25 is the two HiFi molecules exactly | 266 bp shorter, at four AT-rich arrays, all marked uncertain |
| Genotypes: major allele = HiFi | 569 of 570 | 230 of 230 |
| Written back into the graph | +23,779 bp | +5,307 bp |

Details: [reference reconstruction](reference-reconstruction.md),
[validation](reference-reconstruction-validation.md). 0.4.0 full runs of the ten assembly
datasets were unchanged by f6d1745's cleaning (molecules identical).
