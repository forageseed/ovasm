# Short-read reconstruction from species reference graphs

`pan reconstruct` accepts mitochondrial and plastid OV-GFA graphs produced by `ovasm pan merge`.
It uses every source configuration path for recruitment and for read-to-graph projection.
Both reference graphs are required so that shared organelle sequence competes under the
existing Illumina recruitment rules.

```bash
ovasm pan reconstruct \
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

0. Walk the individual's sequence that no reference path holds (frontier walk; skipped
   with `--no-frontier` and in explicit target-input mode). One pass counts every 31-mer
   of the WGS reads in a count-min table (`--frontier-memory-gib`, default 4 GiB). A
   target's depth is the median count of its reference 31-mers. From every place where a
   reference k-mer at that depth is followed by one that is not, the walk takes next
   k-mers from 0.4 of the target's depth to ten times it, follows both of two only when
   the lesser is at least half the greater, and stops after 300 reference k-mers in a row
   (a short stretch the reference holds elsewhere does not end it), where no next k-mer
   is at depth, or at its budget (50,000 k-mers per place, twice the stored reference
   bases per target). Walked pieces shorter than 150 bp that lead nowhere at one end are
   dropped. Nuclear sequence is far shallower than an organelle and a plastid far
   deeper than a mitochondrion, so the walk does not step into either. The walked pieces
   are written to `reference/<organelle>.frontier.fasta` and join the bait only: reads
   still enter the target by the recruitment rules below, and no walked base is emitted.
   Recruitment alone reaches only about a fifth of a read past the reference in each pass,
   and its saturation rule stops after the second pass; extending to a fixed point instead
   walked into the nucleus (see the validation record).
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
   One-base alleles cannot match a long node-internal k-mer. After the exact merge,
   parallel nodes are joined only if their sequence and both oriented neighbour
   sets are identical. This removes duplicate representations of the same allele
   while preserving different alleles, isolated components and tandem self-loops.
   Every affected segment, link and P-walk spelling is checked; the direct node
   mapping is retained in `individual.parallel-nodes.tsv`.
5. Write `individual.gfa` and `individual.fasta` for each organelle. FASTA records are
   maximal nonbranching graph walks. Branches remain in the GFA, and a graph cycle does
   not certify a physical circular chromosome.

The de novo step assembles short reads at several k, as GetOrganelle does: the k the read
length gives (127 for 150 bp reads), then 99, 77 and 55 until one gives a single molecule;
the attempt with the fewest molecules of at least 1 kb, then the most molecule bases, then
the larger k is kept (`k_tried` in its summary).

The result of each organelle is at the top of the output directory, as `ovasm run` gives it:
`<organelle>.gfa` (the de novo organelle graph with its molecule paths) and
`<organelle>.fasta` (its molecules). A single linear molecule of at least 1 kb is first
closed with read pairs, as GetOrganelle's assemblies close on paired reads: pairs with one
mate forward in its last 600 bp and the other reverse in its first 600 bp join the ends, and
their distance follows from the insert size of the pairs within one end. With at least five
such pairs, the reads anchored at either end and their mates spell the bases between, one at a
time by majority (at least three reads and 80% of them, after the longest context of 90, 70, 50
or 31 bases they hold, so an inverted repeat around the join does not mislead them), until the
molecule's first 31 bases are reached. The spelled bases become a segment of the GFA
(`<sample>.read_pairs_join`) linked to the walk's last and first steps, and the walk is marked
circular (`cj:Z:read_pairs`). Recruitment keeps mates together when it downsamples for
reconstruction, or few pairs would remain at a deep plastid.

When the reads do not reach the start but the pairs join the ends, a reference walk may supply
the bases between: exact 31 bp anchors at the molecule's end and start (moved inwards up to
200 bp, with the walk spelling the trimmed bases as the molecule does) at most 5 kb apart, and
a length within 100 bp or 30% of the pairs' distance. Of the walks that agree, the one most like
the individual (the share of the molecule's sampled 31-mers it holds) is taken; the record
names it, its similarity and the most similar walk with a different fill. These bases are the
reference's: the segment carries `rf:Z:<walk>`, the molecule's header `reference_fill_bp`, and
the status is `complete_with_reference_fill`. No bases are filled otherwise, and N is never
used.

`reconstruction.json` gives each organelle a `final` record: `complete` for one circular
molecule of at least 1 kb from a graph without loose segment ends; otherwise `incomplete`, with
every reason. Whether reads decided every repeat pairing is reported beside it (`decisive`).
Short reads cannot decide a repeat longer than their fragments: the pairings of its copies are
then a configuration the graph allows, not one the reads choose, and a plant mitochondrion with
such a repeat stays linear here however much of it the graph holds.

When an organelle is incomplete and its reads leave repeat pairings open, every
configuration its graph allows is listed: each cover of the anchor graph `linearize` finds
(every pairing of every repeat's copies the copy numbers allow), as the molecules it makes,
duplicates removed (a circular walk from any step in either direction, a linear one either
way). They go to `<organelle>.configurations.fasta` (`configuration=<N>`,
`reads_decide=false`) and to `<organelle>.gfa` as P walks
`<sample>_configuration<N>#1#<mt|pt>.<M>` tagged `st:Z:unresolved_configuration`, spelled from
the organelle graph's segments. `reconstruction.json` lists each with its score and how many
are as well supported as the best. The representative molecule in `<organelle>.fasta` is one
of them, not a choice the reads made.

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

`reconstruction.json` reports both organelles, the frontier walk of each (depth, ends,
new k-mers, joins, branches, budget) and actual recruitment counts. A de novo
error is retained explicitly even if a reference projection is available. `complete_assembly`
remains false: local read support alone does not identify complete chromosome paths.
Strict and relaxed scores stay null until an independent compatible truth is supplied.
Physical molecule count stays null and physical topology stays unidentified.

An explicit target-input stage can reuse previously recruited files without scanning the
whole WGS library again:

```bash
ovasm pan reconstruct \
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

The default representative-path search includes the necessary lower bound from
path ends. Its existing four-cover computational envelope is raised only when that
bound exceeds four; otherwise the search interval would be empty without examining
a feasible size. `solve linearize --max-molecules N` with positive N retains an explicit
limit. Omitted larger cover sizes and expansion-budget exhaustion remain reported
and continue to prevent a decisive result. No support, optional-anchor or scoring
threshold is loosened.

An eight-run reference graph is a usable panel, not a universal minimum sample count
or proof of eight independent biological individuals. Reference and target identities
must be checked separately; a target's held-out HiFi graph must not be added to its
reference or used to select paths during reconstruction.

See [the eight-run validation record](reference-reconstruction-validation.md) for
actual input scope, measured stages, first failures, graph results and remaining limits.

The configurations are ranked by the reference runs: for each repeat pairing a configuration
makes (the step before a repeat's copy, the copy, the step after), the reference walks of each
run are searched for the entry's last 31 bases followed, at the distance the copy makes
(within 5%, at least 20 bp), by the exit's first 31 bases, on either strand. A configuration's
`reference_prior_ln` sums ln((runs holding the pairing + 0.5) / (runs holding its entry + 1))
over its pairings; configuration 1 is the one most like the reference runs. Each pairing's runs
are listed, and the FASTA headers and P walks (`rp:f`) carry the score. This is a prior from
other individuals: a repeat that recombines differs between individuals, and the reads of this
one still decide nothing between the configurations.

## Genotypes against the reference graph

Each organelle's reads are also genotyped at every known variant site of its reference graph
(`<organelle>/genotypes.tsv`, `<organelle>/presence.tsv`, `genotypes` in
`reconstruction.json`). A site is a simple bubble: two oriented segments with two or more ways
between them, each through one segment or straight across (an allele lacking that sequence).
Each allele is spelled with 30 bases of either flank, and its depth is the median read count
of its 31-mers that occur nowhere else in the graph. The allele with the most depth is
`major`; another with at least 10% of the site's depth is `minor` (heteroplasmy); a site whose
alleles together have less than a fifth of the organelle's depth (the median of its segments'
own 31-mer depths) is `no_call`. Each allele lists the reference runs whose walks take it.
A segment is `present` from a quarter of the organelle's depth and `absent` below 2%; short
variant segments without 31-mers of their own are judged by their sites instead.

The individual is then laid over the reference runs as a Li-Stephens mosaic
(`<organelle>/mosaic.tsv`, `genotypes.mosaic` in `reconstruction.json`): the runs are the
states, each run's allele at a site is the one its walks take, a run's allele differs from the
call with probability 0.01 and the closest run changes between neighbouring sites with
probability 0.001 (an organelle genome is inherited whole). Sites are ordered along the run
walk that passes most of them; sites off it are predicted from the runs' overall weights. The
record gives each run's weight, the stretches of the backbone each run is closest over
(Viterbi), and every site predicted from the others without its own call: their agreement with
the called sites measures how far the reference runs predict this individual. Sites the reads
could not call are imputed with the predicted allele and its probability, marked `imputed`.
Imputed alleles are the runs', not the reads'.

When the reads leave an organelle's repeat pairings open and some configuration's molecules
of at least 1 kb are all circular, the one with the fewest of them (a master circle when there
is one), the best ranked of those, is the representative: `<organelle>.fasta` holds its
molecules (`chosen_by=fewest_circular_molecules reads_decide=false`; pieces under 1 kb,
the graph's short dead-end tips, are marked `unplaced=true`), and the status is
`circular_configuration_chosen`, with the choice in `representative`. The configuration walks
of `<organelle>.gfa` follow the organelle graph's pieces exactly as `unify` makes them blunt
(a repeat's bases shared with an overlapping neighbour are the neighbour's piece).

The genotypes are also written as `<organelle>/genotypes.vcf` on the backbone walk, haploid:
each site where the individual's allele (called, else imputed: `IMP`, `PP`) or a minor allele
differs from the walk's, with the alleles' read shares (`AF`).

Finally the individual is written back into the species graph: `<organelle>.pan.gfa` is the
reference graph and the individual's representative molecules (their walks renamed to the
sample, `<sample>#1#<mt|pt>.<n>`) merged by `ovasm pan merge` (k 63). Sequence the reference holds is
shared; the rest becomes new segments (`pan.individual_new_bases`).

Stretches of the final molecules whose length the reads do not fix are listed in
`<organelle>/uncertain_regions.tsv` (`uncertain_regions` in the final record): local repeats
(20-mers repeated within 500 bp: tandem arrays, microsatellites) between their nearest unique
20-mers, marked when no read can hold both anchors or under a tenth of the reads that would
at the assembled length do. Each row also gives the read pairs across it and how much longer
they say it is against pairs 1.5 kb away (an estimate only: inserts vary with the place by
about 20 bp). Lengths are not corrected.


## Evidence safeguards (2026-11-03)

A circular graph with repeat stretches whose lengths are not established by reads is
reported as `closed_length_uncertain`, rather than `complete`. It is not repaired by
copying repeat lengths from a potentially different individual's HiFi assembly.
Conflicting reference gap fills are withheld even when one source is globally more
similar. Identical fill candidates can still be used, retaining the reference-fill label.

Automatic pan write-back now requires `complete` and `decisive=true`; otherwise its
report is `withheld_uncertain` and individual outputs remain available. Existing
reference graphs are not modified in place. Pan merge preserves P-walk evidence tags
and propagates reference-filled segment membership to a path-level `rf` marker.
Tagged unresolved/reference-filled paths remain recruitment candidates but are excluded
from reference-prior ranking and gap filling. Legacy graphs whose tags were already
lost cannot retroactively acquire trusted provenance and require a provenance audit.
