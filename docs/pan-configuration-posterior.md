# The configuration a species panel supports (`pan reconstruct`)

Short reads do not span a repeat longer than their fragments. The graph they give is then
right, link for link, but every pairing of a repeat's copies is a walk the graph allows, and
the reads rank all of them alike. `pan reconstruct` used to emit the first cover its search
found. This page describes how the reference graph's other individuals now choose among those
covers, and what the choice is worth.

## What went wrong before

Leave-one-out on 48 *Arabidopsis thaliana* accessions (Zou 2025; each has HiFi and Illumina
reads; the pan graph of a fold is merged from the other 47 HiFi assemblies):

* Every walk of every fold's pan graph carried `st:Z:unresolved_configuration`, because each
  HiFi assembly is open at its two large recombining repeats. `Reference::read` therefore
  took none of them as a prior ("0 reference path runs" in the log), although each of them
  spells how every shorter repeat is joined.
* A closed de novo graph kept its molecule ("A closed but unresolved result is not replaced"):
  the first of, for example, 327,859 covers searched, all at 48 ln 1/2.
* `reference_prior` took its 31-base flanks at the segment ends, which lie inside a link's
  overlap where the assembly graph overlaps, so that all entries of a repeat look alike.

The held-out HiFi reads were then laid on the same graphs
(`scripts/eval_conformations.py`): 2,920 of 2,993 links had at least two reads and 5 were
contradicted, so the graphs are right. The molecules were not: where one was emitted, 38% of
its junctions were not the individual's.

## Objects

`linearize` calls the single-copy segments anchors and the others repeats. An anchor has two
ends. A bridge is a walk from one anchor end through repeats to another anchor end. A
configuration takes one bridge at every anchor end that a bridge leaves, so it is a matching
of the anchor ends: perfect on a closed graph, and with an open end wherever the graph stops
(that end ends a linear molecule). The matching fixes the molecules: leave an anchor by one
end, follow the bridge, enter the next anchor, leave it by its other end.

## What a panel run tells

A run (another individual's assembly in the reference graph) holds a bridge when one of its
walks holds the first anchor's last 31 own bases followed, at the distance the repeats between
make (within 5%, at least 20 bases), by the last anchor's first 31 own bases. "Own" means not
shared with the neighbour through a link overlap.

A walk is evidence only where its own reads decided it. `pan merge` now writes, for every walk
tagged unresolved, the offsets at which the sample's own graph branches (`uj:Z:0,300,400`;
`uj:Z:-` for none): a link that is not the only one at both of its ends. A bridge that spans
such an offset is not held. So a HiFi assembly that is open at two 4-7 kb repeats still
counts for the 100-600 bp repeats inside its unitigs, and a short-read assembly counts for
nothing its own reads could not span. Walks with reference fill or external intervals stay
out, and so does an unresolved walk of a graph merged before `uj` existed (merge it again).
The individual's own walks are never used.

## The model

Nothing below is fitted to results.

1. **Closest relative.** Organelles are inherited whole, so the individual's junctions are
   those of its closest relative in the panel unless they changed since. The posterior weight
   of each run being that relative comes from the individual's major alleles at the graph's
   variant sites (`genotypes.tsv`). A site differs between relatives at an unknown rate below
   one half; with Jeffreys' prior on that rate, d mismatches and m matches have probability
   (4/pi) * integral over t in (0, pi/4) of sin(t)^(2d) cos(t)^(2m). A site the run does not
   pass says nothing (1 / alleles).
2. **Votes and frequency.** At an anchor end with K bridges, `relatives_share` is the weighted
   share of the runs holding each bridge (a run holding several counts for each in equal
   parts) and `species_frequency` the share among all runs alike, with half a count added per
   bridge (Krichevsky-Trofimov).
3. **Conservation.** A junction that changed since the closest relative is drawn from the
   species' frequency: P(bridge) = rho * relatives_share + (1 - rho) * species_frequency. rho
   is the panel's own: each run is predicted from the others in turn, and rho is the fixed
   point of rho = (sum of rho a / (rho a + (1 - rho) b) + 1/2) / (n + 1) over the junctions
   the runs decide (a: the others' weighted share for what the run holds, b: its frequency).
4. **Choice.** The configuration taken maximizes the expected number of correctly joined
   anchor ends: the sum of P over its bridges, at both ends of each. Expectation is linear, so
   this needs no independence between junctions. It is a maximum-weight matching, found
   exactly: the end graph falls into small components (the entries and exits of a repeat, and
   bridges that share a repeat), and each component's perfect matchings are enumerated
   (2,000,000 steps at most; beyond that the panel's cover is not built and the report says
   so).

The individual's reads come first throughout. A matching is ranked as `linearize` ranks a
cover on reads: bridges that reads measured before the graph's pairings, then the shares of
reads less the penalty per repeat copy laid down beyond or short of its depth. The panel
orders only what those leave tied, and among matchings the panel ranks alike the one with the
fewest molecules is given.

## What comes out

* `<organelle>.fasta`: the panel's configuration when the panel tells bridges apart at one
  open end or more, `chosen_by=pan_posterior reads_decide=false
  expected_junction_accuracy=<share> open_junction_ends=<n> panel_decided_ends=<n>`. Status
  `pan_configuration_chosen` (`..._length_uncertain` with stretches of uncertain length). Where
  the graph stops at two ends the molecule is linear between them, and the join that read
  pairs gave the de novo molecule closes it. A reference-filled closure is kept as it was.
* `reconstruction.json`, `final.configurations.pan_posterior`: panel runs, variant sites, the
  five closest relatives with their weights, rho and the junctions it rests on, the cover's
  expected share of correct ends, how many configurations rank as well, and every bridge of
  an open end with its probability, relatives' share, species frequency, runs holding it and
  whether it was taken.
* The configuration is number 1 of `<organelle>.configurations.fasta`. All walks stay tagged
  unresolved (`st:Z:unresolved_configuration`, `up:Z:unresolved`) and are not written back
  into the species graph: a panel's choice is not an observation of the individual.
* `ovasm solve pan-configuration --gfa de-novo/assembly.gfa --evidence de-novo/evidence.json
  --reference species.gfa --genotypes genotypes.tsv --sample NAME --out-fasta F --out-json J`
  gives the same choice from the files a reconstruction left, without the reads.

`expected_junction_accuracy` is the model's own statement: the mean of P over the anchor ends
of the configuration. An end where the panel is silent counts 1/K.

## Validation

The 48 folds above, mitochondrion. The new choice was made from each fold's own files
(`ovasm solve pan-configuration` on the de novo graph, evidence and genotypes that ps-2107407
left, against the fold's pan graph merged again so that its walks carry `uj`), so the graphs
are the same and only the walk through them differs. Two answers that took no part in either:
the held-out accession's own HiFi assembly (the partner of every anchor end in its walks; an
end where the HiFi graph is open has several right partners), and its HiFi reads on the same
graph (`eval_conformations.py`: a traversal of a repeat with under 5% of the reads that span
the repeat is a minor form).

| | ps-2107407 molecule | panel's configuration |
| --- | ---: | ---: |
| molecules, all folds | 443 | 116 |
| junction ends with a HiFi-assembly answer | 2,276 | 3,240 |
| of them wrong | 866 (38.0%) | 78 (2.4%) |
| repeat traversals the HiFi reads were counted on | 742 | 1,011 |
| of them minor forms | 309 (41.6%) | 24 (2.4%) |
| folds without a wrong end | 2 | 27 |
| truth covered, median (at least 95% in) | 96.2% (26) | 100.0% (41) |
| identity over primary alignments, median (at least 99% in) | 84.5% (20) | 100.0% (42) |
| one molecule covers at least 95% of the truth | 14 | 34 |

The last three rows are the scores of the first report (minimap2 asm5 against the HiFi
assembly's representative, blocks of at least 1 kb; an alignment that runs across a
rearrangement loses identity). Against the best of the configurations the HiFi graph itself
lists, identity is at least 99% in 48 of 48 (before: 20).

The model's own statement was cautious: it expected 90.8% of the open ends right (ends where
the panel is silent count one half; it told bridges apart at 2,603 of 2,925 open ends), and
97.6% of those with an answer were. rho was 0.78 (0.65 to 0.84 over the folds), from 1,358 to
5,079 junctions of the panel each. No fold needed more than 284 matchings to be enumerated.

What is still wrong, by cause: an arrangement only the individual has (SAMEA113427122, 8
ends); a closest relative that is one variant site closer than five others and joined
differently (run ERR11437347 for four accessions, 4 ends each); a closest relative with a
change of its own; and, in the graphs that lack a link, the place where the molecule is left
open.

Not covered by this validation: other species, panels of fewer runs, and panels whose runs
are short-read assemblies (they hold no junction across a repeat their own reads did not
span, so they add nothing and take nothing away).

## Limits

* An arrangement no panel run has cannot be predicted (SAMEA113427122 differs from all 47
  others at repeats of 181, 364 and 508 bp; its HiFi reads confirm it). The report shows such
  junctions as confidently wrong only after the fact; before it, `species_frequency` near 1
  is all that can be said.
* The closest relative by alleles can carry a change of its own. rho is one number for the
  whole graph, estimated from how often that happens within the panel.
* Where every panel run is open (the large recombining repeats) the panel says nothing and
  the cover with the fewest molecules is given, as before. Both pairings exist in the plant.
* An optional (shallow) anchor is left out unless an end beside it would stay open, as the
  search does. Where the bridged ends of a repeat cannot all be joined (a graph that lacks a
  link), the fewest are left open, two at most per repeat; a graph that needs more gets no
  panel cover and the searched covers stand as before.
* The distance tolerance (5%, 20 bases) is `reference_prior`'s; a repeat whose copies differ
  more in length between individuals is not recognized.
