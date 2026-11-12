# OVASM: mathematical and algorithmic principles

Scope: `ovasm run hifi` in full detail; `ovasm run short` and `ovasm run ont|hybrid` briefly. Source read: branch
`agent/protein-seeds`, `src/*.rs`. Each section names the file and function that holds the rule. A threshold is marked
**derived** when the code or its docs give an argument from the data model, **empirical** when it was set on named
samples, and **unstated** when no justification is written down. Lengths are in bp. "Depth" has three meanings, always
qualified: *read depth* (coverage), *sketch depth* (count of a read's sampled 21-mers over the whole input) and
*node depth* (count of a graph node, `assemble.rs`).

## 0. Pipeline in run order

```
reads -> depth gate -> protein seeds or seed library -> recruitment (passes, bait growth) -> downsampling
      -> sparse graph + k ladder -> component filter -> read evidence -> linearize (molecules)
      -> unified graph -> chromosomes -> assess -> result.json
```

`main.rs::run` builds the seeds first (section 2) and then calls `pipeline::run`. That validates the config, runs
`recruit_all` once for all requested organelles, then `assemble_organelle` per organelle (output in
`<out>/<organelle>/` when both are requested). A failure of one organelle does not stop the other; `summary.json` records
both and the run returns an error if any failed.

## 1. Inputs and read sets

Source: `main.rs` (`RunArgs`), `pipeline.rs` (`ReadType`, `ReadSet`, `RecruitBy`), `recruit.rs::read_batches`.

- Files: FASTQ or FASTA, optionally gzipped (needletail). Several `--reads` are one library. Batches are cut by bases,
  not reads (`BATCH_BASES` = 8 Mi), so the work per batch is the same for 150 bp and 15 kb reads (2,048 short reads
  were too little work for 32 threads; empirical).
- `ReadType`: `hifi` (sparse closed-syncmer graph), `sr` (dense graph, every k-mer a node); `ont` and `hybrid` are flows
  that end in the HiFi stages (section 7).
- `ReadSet`: `whole-genome` (the organelle's reads must be recruited) or `target` (reads already are the organelle's;
  recruitment and protein seeds are skipped and the reads go to assembly as given).
- `RecruitBy`: `seeds` (default), `discover` (no seeds; depth clusters called by the protein panel, HiFi only), `both`
  (seed reads plus clusters called as the target).
- `--lineage landplants|metazoa|fungi` selects the seed library (`seeddb.rs`), the protein panel and genetic code
  (`panel.rs`) and several defaults below. Libraries: land plants 16 mitochondria and 44 plastids, metazoa 32
  mitochondria, fungi 19 mitochondria. Animals and fungi have no plastid, so the organelles sought are filtered by
  `lineage.fasta(o)`.
- The depth-gate options (`--genome-size`, `--nuclear-depth`, `--extend-gate-fold`, `--no-extend-gate`) are accepted for
  HiFi whole-genome runs only (`main.rs::run`).
- Fewest recruited bases before assembly is attempted (`pipeline::min_recruited_bases`): 100,000 for land plants
  (empirical: Luzula sylvatica recruited 4 reads, 4 kb, and failed later), 30,000 for animals and fungi (empirical: Canis
  lupus at 10 Mb of reads recruited 11 reads, 75 kb, about four mitogenomes).

## 2. Depth gate and protein seeds

Protein seeds replace the seed library when the sample may be far from every library genome. They are the default for
every lineage on HiFi whole-genome reads (`--no-protein-seeds` turns them off), and never used with `--recruit discover` or
`--read-set target` (`main.rs::run`). An organelle without a called cluster falls back to the library
(`seeddb::fill_missing`).

### 2.1 Depth gate (`seedfree.rs::DepthGate`, shared with `find discover` and recruitment)

Quantity: which reads are organelle-like, i.e. come from sequence present at many times the nuclear copy number.

1. Sketch. Canonical 21-mers are kept when `mix64(kmer) < u64::MAX/64` (FracMinHash, scale 64). The whole input is read
   (`sketch_bases = 0`), because a prefix sketch failed on a run whose plastid reads all came after the first 3 Gb
   (empirical, doc on `DiscoverParams::sketch_bases`). `count(kmer)` is its number of occurrences.
2. Read features (`features_of`). For a read with sampled k-mers $K$, keep depths $d_i = \mathrm{count}(K_i)$ with
   $d_i \ge 2$ (a k-mer seen once carries the read's own error). Let $q_1$, med, $q_3$ be quantiles of the kept depths,
   distinct $= |\mathrm{set}(K)|/|K|$, high $=$ share of kept $d_i \ge T$, spread $= (q_3 - q_1)/\mathrm{med}$.
3. Organelle-like (`passes`): $|K| \ge 20$, distinct $\ge 0.8$, high $\ge 0.9$, spread $\le 1.5$, med $> 0$. Intent: deep
   along the whole read and uniform, which excludes NUMT/NUPT reads (single-copy flank) and tandem arrays (repeat their
   k-mers in the read). The four numbers are `long_read_gate_params`; no evidence is given for 20, 0.8, 0.9 or 1.5
   (**unstated**; they are `find discover`'s defaults).
4. Threshold $T = f \cdot D_n$ with fold $f = 3$ (**unstated** beyond the example in the module doc: Nipponbare
   mitochondria ~433x, plastids ~2,665x, nucleus ~80x) and $D_n$ the nuclear depth.
5. Nuclear depth $D_n$ (`nuclear_read_depth`). The median sketch depth of every judged read ($|K| \ge 20$, median
   $\ge 3$ = `MIN_JUDGED_DEPTH`) is binned on $\log_2$ with width 0.25 (`DEPTH_BIN`, 19% of the depth, so a
   Poisson-spread peak holds the same share at any depth: **derived**). Bins are smoothed with their two neighbours and
   candidate bins must be local maxima of (smoothed, raw). Take the fullest candidate whose reads within $\pm 2$ bins
   stand for at least `MIN_NUCLEAR_GENOME` = 20 Mb, where the sequence a window stands for is
   $\sum \mathrm{bases} / \mathrm{depth}$ and depth is the window's read median. The 20 Mb floor is **derived** from
   genome sizes (smallest land-plant nuclear genomes ~60 Mb, mitochondria at most ~12 Mb); no animal or fungal check is
   cited. Reads are used rather than k-mers because a few percent of shallow reads (contaminants, 1-3x sequence) bring more
   distinct k-mers at count 2 than any count of the nuclear peak holds; the k-mer estimate gave 2x for Carex echinata,
   Medicago polymorpha and Lophophora (empirical, doc on `nuclear_read_depth`).
6. Fallbacks: a separate low-depth integer peak (`low_depth_read_peak`: count at depth $d$ exceeds both neighbours by more
   than $2\sqrt{n_d + n_{d\pm1}}$, $3 \le d < 10$, at least 20 Mb of sequence); the mode of the read-context k-mer histogram when
   no read peak stands for a nucleus (`nuclear_source = "kmers"`). A warning is recorded when the read peak is a flat
   shoulder (`read_peak_shape_supported`: its mass must exceed the masses four bins above and below).
7. Overrides: `--nuclear-depth D` sets $D_n = D$; `--genome-size G` sets $D_n = \mathrm{round}(\text{sketch bases}/G)$.
   Needed where read medians stand at two copies (Lespedeza davurica: read peak 56x, k-mer peak 28x, the 72x
   mitochondrion under the 3x threshold; empirical).

### 2.2 Depth clusters (`seedfree::depth_peaks`, used by `protein_seed::run`)

Quantity: the separate copy-number populations among the organelle-like reads.

Histogram of $\log_2(\text{read median depth})$ with bin 0.25. A bin $b$ is a peak when its smoothed count
$s(b) = n_{b-1} + n_b + n_{b+1} \ge 20$ (`min_cluster_reads`, **unstated**) and $(s, n_b)$ is $\ge$ the left and $>$ the
right neighbour's (the raw count breaks single-bin plateaus). Peaks closer than 1 $\log_2$ unit are merged by averaging
("one genome, not two"; **unstated** beyond that). Each organelle-like read joins the peak nearest its median in $\log_2$.
Peaks are sorted deepest first (usually the plastid).

Cluster span $= \sum \text{read bases} / 2^{\text{peak}}$, the sequence the cluster stands for. A cluster seeds only if span
$\ge$ `MIN_CLUSTER_SPAN` = 4,000 (`protein_seed.rs`). Evidence: a deep nuclear repeat holding organelle-derived genes stood
for about 36 bp (Orcinus orca, 101,070x, 356 reads); the smallest mitogenomes are ~6 kb (empirical, one negative case and a
genome-size fact).

### 2.3 Six-frame translation and gene stretches (`panel.rs`)

Quantity: which panel gene, if any, a stop-free peptide stretch of a read belongs to.

- Genetic code (`genetic_code(table)`): 1 for plants, 5 for metazoa, 4 for fungi. TGA is Trp in 4 and 5, so the standard
  table would cut animal and fungal genes into pieces shorter than `min_orf`. Table 2 is available but no lineage default
  uses it, so vertebrate AGA/AGG (stops in table 2) are read as Ser.
- Each read is translated in six frames. A stretch is a maximal stop-free run of at least `min_orf` = 80 residues. About 2% of
  stretches of random sequence reach that length (**derived**, doc in `panel.rs`).
- Lookup: amino-acid k-mers with $k = 6$ (5 bits per residue); the panel index maps a 6-mer to the genes holding it. The
  panel's 6-mers are about 1% of all 6-mers, so chance hits stay far below `min_hits` (**derived**, same doc). For a stretch,
  $h_g = \#\{i : \text{6-mer}_i \in \mathrm{genes}(g)\}$; the best gene is $g^* = \arg\max h_g$ (ties to the lower gene id, so the
  call does not depend on hash order) and the stretch is assigned only if $h_{g^*} \ge$ `min_hits`.
- `min_hits` (`Panel::lineage_min_hits`): 5 for plants and metazoa, 6 for fungi. Rule: the lowest value at which random DNA
  (six frames, 10 Mb, GC 0.30 and 0.45) gives at most 0.2 gene stretches per Mb while keeping the most held-out genes
  (`scripts/panel_tradeoff.py`, 2026-11-05). Empirical. The land-plant panel gives 0.4 per Mb at 5, i.e. it exceeds the stated
  0.2 criterion.
- Significance (metazoa and fungi only; `Panel::max_evalue`, `panel::significant_bits`). Shared 6-mers pick the gene but
  do not by themselves show homology: hydrophobic, low-complexity peptides from AT-rich nuclear or bacterial DNA share many
  6-mers with the short membrane proteins (atp8, nad4l, nad6, nad3). A stretch of $m$ residues keeps $g^*$ only if its best
  ungapped BLOSUM62 alignment to a protein of $g^*$ (extended from every diagonal of a shared 6-mer, X-drop 16;
  $\lambda = 0.3176$, $K = 0.134$), less what its composition alone scores, reaches E $\le 10^{-6}$ against the panel's $n$
  residues: $S_\mathrm{fwd} - \max(0, S_\mathrm{rev} - \log_2 mn) \ge \log_2(mn / 10^{-6})$ bits, where $S_\mathrm{rev}$ is the
  same score of the stretch reversed (composition kept, homology gone). Empirical check: of 84 assigned stretches of seven
  reference mitochondria (two sponges, a sea anemone, two bryozoans, a moth, a sawfly) 83 are kept (60-885 bits forward, at
  most 35 reversed; the one lost, a sawfly nad6, scores 20 bits). Of 403 assigned stretches of the whole-read seeds of two
  PMAT2 read sets without the mitochondrion (Calliactis parasitica, Myriapora truncata; DIAMOND blastx found homology in none
  of 31 and 3 of 46 seed reads) only the 12 of the real homologs are kept: 348 fall short forward, 43 reach the threshold
  forward (up to 120 bits) but score as much reversed. The land-plant panel has no such check (`max_evalue` None).
- Metazoa and fungi panels are small covers chosen by leave-genus-out greedy selection (`scripts/polish_panel.py`; no genus of
  the PMAT2 benchmark): 551 metazoan and 696 fungal proteins. Larger panels need higher `min_hits` and find slightly fewer
  held-out genes (metazoa 551 proteins at 5 hits 98.96%, 13,868 at 15 hits 98.08%; fungi 696 at 6 hits 99.34%, 8,650 at 10 hits
  99.10%).
- Land-plant panel (`data/organelle_proteins.faa`): mitochondrial proteins and plastome CDS translations, headers
  `class|gene|species`.

### 2.4 Read class, cluster call, seeds (`protein_seed.rs::run`, `panel::call_from`)

Four passes over the input: depth sketch; nuclear peak; organelle-like reads with their medians; scan of the sampled reads.
At most `sample` = 3,000 reads per cluster are scanned, every `step`-th in input order, $\mathrm{step} = \max(1, \lfloor n_c/3000
\rfloor)$ (**unstated**; a cost cap).

1. Read class (`read_class`). With $m$ stretches of mitochondrial genes and $p$ of plastid genes, the read is mito if
   $m \ge$ `min_read_stretches` (2) and $m > p$; plastid symmetrically. Evidence for 2: with one stretch, 84 of 1,452 Polytrichum
   commune "mitochondrial" seeds were on no organelle (empirical).
2. Gene tallies per cluster count only the stretches of the read's own class, so a homologous gene of the other organelle or an
   insert does not add to its class.
3. Per class and cluster: genes $G$ = distinct genes with at least `min_gene_segments` = 2 stretches; density $\rho$ = stretches of
   those genes per scanned read. The class counts when $G \ge$ `min_genes` = 5 and $\rho \ge$ `min_density` = 0.3. Density evidence:
   organelle reads carry plastid 2-9 and mitochondrion 0.7-3 stretches per read, organelle-derived nuclear sequence 0.03-0.16
   (empirical, v3 clusters). `min_genes` = 5 has no stated evidence.
4. `panel::call_from` returns plastid, mitochondrion, mixed or none: a counting class is the call if the other does not count or if it
   has at least twice the other's genes; both counting without that margin is mixed. `ClusterCall.seeded` does not use this call:
   each counting class seeds independently (Ulota crispa: part of the mitochondrion at the plastid's 1,878x).
5. The reversed control scans the same read with its sequence reversed. A gene contributes its forward count minus its
   reversed count when the forward count is at least twice the reversed count. Otherwise it contributes zero. A class that
   meets the gene and density thresholds in the reversed scan does not seed from the cluster.
6. Seed = a read cropped to the span of its class's stretches plus `flank` = 2,000 on each side (**unstated**). Seeds go to
   `protein_seeds.<organelle>.fasta`.
7. Depth reference: for each organelle, the seeded cluster with the most genes (ties: the deeper) writes
   `protein_seeds.<o>.depth.fasta`, used for assembly node depth (section 4.3).

The `protein_seeds.json` report includes the loaded panel families per organelle in `panel_gene_sets`. Gene labels
are the exact panel labels. Each depth cluster reports per-family stretch and FASTX record counts before control
subtraction and minimum-stretch filtering. It also counts records with the same family detected in both scans.
These are input records. Physical molecule identity and source require separate evidence. The whole-read
fallback is outside this paired report because it scans forward only.

#### Export candidates from paired gene sets

`ovasm find protein-seeds --paired-gene-candidates` exports nominated depth-cluster
crops in `protein_candidates.<organelle>.fasta` under a new output directory.
Each panel family needs more
forward than reversed input records and at least `min_gene_segments` net forward
stretches. A cluster nominates the requested organelle when its supported family
count reaches `--min-genes`. The report `protein_candidates.json` includes the
effective family denominator, the paired evidence, and the legacy cluster calls
for comparison on the same sampled records. Density and estimated cluster span
are recorded for comparison.

These crops have pending source qualification. Use the original reads and the
registered source windows with `audit_seed_source.py` and
`--retention-policy target-supported-only` before passing the qualified FASTA and
its audit into recruitment. The ordinary seed path includes whole-read fallback.

### 2.5 Whole-read scan (`protein_seed::whole_read_seeds`)

For organelles no cluster seeded, every read is scanned with $\text{min\_hits}' = \max(\text{min\_hits}, 12)$
(`WHOLE_READ_MIN_HITS`). A read is a seed when its class is that organelle and it holds at least `WHOLE_READ_MIN_GENES` = 3 distinct
genes of the class. The seed is cropped to the gene span with the configured flanks, up to `WHOLE_READ_MAX_SEEDS` = 3,000 seeds. Evidence for 12: chance stretches of nuclear
reads reach 5-8 hits, mitochondrial genes 13-140 on the best stretch (PMAT2 reference genomes); at 5 hits and 2 genes, 1,514 of
74,216 Gibbula magus reads passed against 19 mitochondrial (empirical). Cases: organelle about as deep as the nucleus (Orcinus orca,
19 mitochondrial reads in 1 Gb; Sterechinus neumayeri, 9 in 2.5 Gb). The value 3 rests on a rule of thumb in the constant's doc (a NUMT
is a few kb of one or two genes; a 10 kb animal HiFi read holds five or more); no distribution is cited. If any whole-read seed is found
(`ProteinSeedReport::shallow`), recruitment runs without the depth gate and with one pass.

### 2.6 Audit seed source with original read coordinates

`scripts/audit_seed_source.py MANIFEST --out AUDIT.json` checks native protein seed crops
against predeclared source windows. Each seed header supplies an original read identifier
and a `span=START-END` interval. The tool requires a unique original record and exact crop
sequence at those coordinates. It also requires an exact, unambiguous source window
inside the original read and the crop.

The manifest uses `ovasm.seed-source-audit.v1` and names the comparison scope. Its `seeds`
and `raw_reads` entries lock the seed FASTA and each original FASTX file by path and SHA-256.
A `raw_reads` entry can select a one-based `id_field` when an archive prefixes its own ID.
The `references` entries add an ID and either a `target` or `competing` role to a locked FASTA.
Every registered `windows` entry provides an ID, reference ID, record ID, half-open
reference interval, interior boundary coordinate, and a locked `provenance` file.
Select intervals that cross source-informative boundaries before evaluating crops.

`--out-fasta RETAINED.fa` writes a separate FASTA that omits crops with a supported competing
context. It preserves target-supported and unresolved records byte for byte. The audit
records exact placements in the named panel, both strands, and uncertainty from repeated
IDs, repeated windows, mixed source observations, or missing crop evidence. Source claims
are scoped to the supplied references and provenance. Use an independent evidence set
when validating the windows or measuring the source filter.

For positive-only release, add `--retention-policy target-supported-only`.
The producer retains only crops that spell a unique target window completely
within their audited original span. If original IDs repeat, strict mode
requires the ID, span and literal crop to match exactly one original record.
The audit records its file index, record ordinal and sequence digest. Missing,
ambiguous, mixed or competing evidence leaves the crop out of the strict view.
An empty strict view is reportable. OVASM requires at least one qualifying crop
for that organelle before starting recruitment.

For an audited view prepared with a manifest target of `mitochondrion` or `plastid`, pass
the retained FASTA to normal seed recruitment. Supply the audit with the same organelle.

~~~sh
ovasm run hifi --reads sample.fastq.gz --out assembly \
  --seeds mt=RETAINED.fa --seed-source-audit mt=AUDIT.json
~~~

OVASM checks the digest, target, and complete ordered list of retained records.
For a strict view it also checks original-coordinate spans and the literal
target window spelling within each crop. This check uses the source evidence
in the named audit panel and applies to the organelle named on the command. It uses
a copy of the verified bytes in the run directory and records the audit in
`summary.json` under `seed_database.seed_source_audit`. The audit describes the source
comparison that produced the seed view. The run may use a different reads file.

#### Recover original records from candidate crops

Run `python3 scripts/extract_seed_crop_originals.py PLAN.json --out AUDIT.json --out-reads ORIGINALS.fastq`
before repeated source comparisons on the same crops. Use schema
`ovasm.seed-crop-originals-plan.v1` with `scope`, a `candidate_sets` list of
FASTA file bindings with an `id` per set, and a `raw_reads` list of original
plain or gzip FASTQ file bindings. Each binding has `path` and `sha256`.
An original input can also set its one-based `id_field`.

The extractor verifies the archive SHA-256 while scanning the entire file
and joins each crop by identifier, exact interval and literal sequence bytes.
It records missing, mismatched and ambiguous matches. The output includes each
uniquely bound original once, with its full FASTQ title, sequence and qualities.
The report ties each crop to an input record ordinal and an output record index.
Use these bindings when assigning mapper names. Source qualification requires
separate reference and read evidence.

### 2.7 Audit protein gene sets

Use `ovasm find protein-seeds --gene-set-audit` or add `--protein-gene-set-audit`
to an active protein seed scan in `ovasm run hifi`. The report records gene sets for
each sampled depth cluster. If the whole-read fallback runs, it also scans every
fallback input in reverse sequence order, using the same translation, gene
assignment, read qualification, and crop rules in both orientations. The fallback
counts qualifying records before the seed output cap and depth check.

Each gene set lists observed forward families and families with both more forward
than reversed records and at least `min_gene_segments` excess forward stretches.
The shared record count is provided for each family. `support_fraction` is the
supported family count divided by the indexed panel families for that organelle,
or null if none are indexed. A family present in the loaded panel with no amino
acid k-mer is reported in `unindexed_families`.

The depth-cluster and fallback scans qualify reads differently. Interpret their
fractions within each scan. Counts refer to input records, and multiple families
can share a record. The reversal control measures directional signal under these
rules. True homologues in nuclear and plastid copies can also produce that signal,
so assigning physical source requires independent original-read evidence. This
audit leaves seed selection and FASTA output unchanged.


## 3. Seed recruitment

Source: `recruit.rs`, `bait.rs`, `seeddb.rs`, `main.rs::run`, `pipeline.rs::recruit_all`.

### 3.1 Bait set (`bait.rs`)

- Bait = map from canonical $k$-mer to a bitmask of the targets (organelles) it came from, at most 8 targets. $k = 21$ for HiFi and short
  reads, 15 for the `ont` preset (`run ont` sets 21; section 7).
- Seed k-mers: all canonical k-mers of the seed sequences. The target's `seed_length` is the summed seed length and serves as the
  genome-size estimate until the depth sketch gives a better one.
- `add_seed_reads` (reads of the same sample as bait, e.g. corrected long reads): solid k-mers are those seen in at least `min_count` = 3
  reads (counted once per read) and at least `SEED_PEAK_FRACTION` = 0.25 of the count peak (`solid_kmers`; the peak is the maximum after
  the valley of the count histogram, bins to 4,097 with the last bin holding the rest so a ~10,000x target is not misread). 0.25 is
  empirical: at 0.1, NUMT flanks recruited 16,079 nuclear short reads (Nipponbare ONT) and gold 21-mer recall moved from 99.65% to 99.64%
  when going from 0.1 to 0.25.

### 3.2 Read classification (`recruit::classify`)

For a read and target $t$, with bait hits at positions of the read's 21-mers:

- Hits closer than `block_gap` are bridged into a block. $C_t$ = total block length, $H_t$ = number of hits, $L_t$ = longest block,
  $N$ = k-mers in the read.
- Best target $= \arg\max_t (C_t, H_t)$ lexicographically. If several tie on both, the read goes to each (identical MTPT/NUPT copies).
- cover $= C_{best}/|\mathrm{read}|$, hit share $= H_{best}/N$.
- Call `Organelle` if cover $\ge$ `min_cover` and hit share $\ge$ `min_hit_frac`, unless the read is a tandem array; else `NumtCandidate` if
  $L_{best} \ge$ `numt_min_block`; else `Other`.

Defaults (`main.rs::Preset::defaults`, with `min_cover` = 0.8 for all):

| preset | k | block_gap | min_hit_frac | numt_min_block |
|---|---|---|---|---|
| hifi | 21 | 50 | 0.2 | 1000 |
| sr | 21 | 40 | 0.2 | none |
| ont | 15 | 200 | 0.05 | 2000 |

Stated rationale: blocks bridge gaps so sequencing errors and seed/sample SNVs do not break coverage while NUMT flanks (kilobases of
nuclear DNA) do. The numbers themselves have no stated evidence.

Tandem arrays (`HIFI_MIN_DISTINCT` = 0.5): over an eighth of the read's k-mers (chosen by hash) the share of distinct ones must be at least 0.5
for HiFi; 0 (off) for short and ont reads. Evidence: Medicago truncatula's second pass recruited 115,110 reads, 90% of them a nuclear
satellite with 6-16% distinct 21-mers; organelle reads are near 100%; a read across two tandem copies of an organelle segment is still
about two thirds distinct (empirical).

### 3.3 Passes and bait extension

Each pass scans the whole input with the current bait. After it, a target's bait grows by solid k-mers of the reads new to it
(seen in at least `extend_min_count` = 3 reads and at least 0.1 of the count peak; the "noisy reads run away" argument is in `bait.rs`).
Controls:

- Subsampling. New reads are thinned to `extend_depth` = 200 over the genome they stand for. With $d$ = bases/genome and $n$ = reads,
  the kept fraction is $\phi = \min(1, \max(200/d, 5000/n))$ (`extension_fraction`). The 5,000-read floor (`MIN_EXTEND_READS`) is
  empirical: the Schoenoplectus lacustris plastid, a distant seed, stopped at 9 kb when 9 reads extended it. The genome is the size the
  target's reads stand for when the depth sketch exists (`implied_size`), otherwise the seed length.
- Extension gate (HiFi, `--extend-gate`; on in `run` unless `--no-extend-gate` or the whole-read scan seeded). A read extends a target
  only if it is organelle-like by the depth gate with deep threshold $T_e$ below, has sketch median at least `EXTEND_CORE_FLOOR` = 0.25 of
  the target's core depth (median read depth in the first pass with at least 20 reads, `CORE_MIN_READS`) and at most `HIGH_COPY_FOLD` = 10
  times it. With core depth $c$, nuclear depth $n$ (`extension_nuclear`) and gate organelle threshold $3n$:
  $$T_e = \min\bigl(3n,\ \sqrt{nc},\ \max(0.75\,c,\ n)\bigr).$$
  The geometric mean lowers the threshold for targets less than nine times as deep as the nucleus (Lophophora williamsii mitochondrion:
  core 104x, nucleus 22x, 66x threshold stopped it); the $0.75c$ cap keeps nine tenths of a read's k-mers above it (97% of Poisson counts
  at 72x are at least 54). The geometric mean and the 0.25 and 10 are empirical (named samples); 0.75 has the Poisson argument.
  A read that fails as a whole is retried on only the k-mers the bait lacks (`extends`), so reads across a target's own deep tandem
  arrays pass (Schoenoplectus plastid inverted repeat, 30-40x deeper than the rest).
- `extension_nuclear`: the k-mer peak replaces the read peak when the ratio read peak / k-mer peak lies in [1.7, 2.5] and the k-mer peak is at
  least 10 (`NUCLEAR_DOUBT_RATIO`, `TWO_COPY_RATIO`, `MIN_KMER_PEAK`; a read median at two copies, Lespedeza davurica). The three cut-offs
  rest on that one species.
- Core depth (`target_depth`): the target's median read sketch depth, unless that exceeds `HIGH_COPY_FOLD` times the span-weighted median
  of its organelle-deep reads (each weighted by length/depth), in which case the latter. Evidence: Lophophora plastid first pass, a repeat
  at ~700,000x beside plastid reads at ~440x gave a 24 bp "genome".
- Stop rule (`settled`, `genome_settled`). A target stops when its reads grew by less than `saturation` = 0.05 of what it had, unless the
  gate is on and the growth is still at least `STEADY_GROWTH` = 0.5 of the previous pass's growth; it also needs the implied genome to
  have grown by less than the mean read length. Rationale: a large mitochondrion gains about as many reads every pass (Lophophora, 2.2 Mb,
  250-330 reads per pass, stopped at 418 kb by the 5% rule alone). 0.05 and 0.5 are empirical.
- Frontier extension within a pass (`FRONTIER_STEPS` = 10, at most `FRONTIER_READS` = 10,000 per target): partial reads that the grown bait now
  covers by `min_cover` extend it again in the same pass. Motivation: one pass moves an edge by about a tenth of a read length. Both numbers
  are unstated beyond that.
- Size cap (`size_cap`): a plastid target stops when its reads stand for more than 1,000,000 bp, a mitochondrion more than 20,000,000 bp
  (four times the largest plastome, 242 kb; nearly twice the largest mitogenome, 11.3 Mb): derived from genome sizes. Stop reason `size_cap`
  counts as not converged.
- Bootstrap (`--bootstrap`, HiFi). When a seed is too distant, deep partial reads (NUMT candidates that are organelle-like) seed it once
  every other target has stopped. A candidate must be no shallower than 0.5 of the target's own reads (`BOOTSTRAP_DEPTH_FLOOR`, needs 10
  reads) and covered no more than 0.5 by another target's bait. Leftover bootstrap (mitochondrion only): at least 5 reads at 0.5-2x the
  target's depth (`LEFTOVER_DEPTH`) to find the other circles of a multipartite genome (Buxus 1.05, Bellardia 0.76, Euphorbia 0.71-0.90 of
  the published circles' depth). The plastid is excluded because what is left at its depth is nuclear rDNA or satellites (Pseudognaphalium ran to
  the 1 Mb cap). Empirical.
- Short reads (`--kmer-gate`). The library's k-mers are counted once (count-min, $2^{31}$ 16-bit counters, 4 GiB). A read extends a target only
  if the median count of its k-mers that the bait lacks is within [0.25, 10] times the target's k-mer depth (`KMER_GATE_FLOOR`,
  `KMER_GATE_CEILING`); fewer than 5 new k-mers extends. Evidence: ungated WL354 grew from 0.72 M to 54 M reads in 20 passes.

### 3.4 Passes per lineage (`main.rs::protein_seed_passes`)

| case | `--recruit-iterations` |
|---|---|
| land plants with protein seeds | 60 (`GATED_RECRUIT_PASSES`), stopping once the bait stops growing |
| metazoa, fungi with protein seeds | up to 3 (`PROTEIN_SEED_PASSES_SMALL_MITO`); fewer once no read is added |
| seeds from the whole-read scan (`shallow`) | 1, any lineage |
| seeds from the library | 60 (`GATED_RECRUIT_PASSES`) for gated HiFi, else 20 |

Empirical evidence. Plants: a land-plant mitochondrion of 200-500 kb is mostly intergenic, and seeds are genes $\pm$ a read. One pass
left it partial on 9 of 10 DToL angiosperms (Oatk molecule covered: Juncus bufonius 0.27, Euphorbia lathyris 0.13), where the usual gated
passes reach 0.999 and 0.886. Mosses (the earlier case for one pass, recruitment precision after 1, 2, 60 passes Polytrichum commune
0.92/0.59/0.59) assemble the same after 60 (Polytrichum 0.9996): the extra reads are not kept in the assembly. Animals and fungi: seeds are genes $\pm$ 2 kb and a
mitochondrion whose genes sit apart misses part of the molecule after one pass (Rhizopus microsporus 186 reads for 29 of 43 kb, Ganoderma
resinaceum 60 of 87 kb, Bibio marci 13 of 16 kb). Whole-read seeds are capped at one pass because later passes enter the nucleus (Orcinus
orca: 237 reads after one pass, 7,237 reads/660 kb after three, killed at 30.6 GB). The 60-pass cap is a bound, not derived: a read joins only
when the bait covers 80% of it, so an edge advances about a tenth of a read per pass (Carex echinata: 27 kb per pass, 830 kb after 20).
When protein-seed passes exceed one, the mitochondrial component filter also runs in gene-free mode (section 5).

### 3.5 NUMT flags and competing targets

- `NumtCandidate` describes a partial read (cover below `min_cover`) with a contiguous organelle-like block of at least `numt_min_block`.
  The unexplained sequence may reflect a nuclear insertion or sequence missing from the target bait. OVASM writes these reads to
  `numt_candidates.fastq` and uses them in bootstrap and frontier extension.
- Both organelles are seeded unless `--single-seed`: a mitochondrion seeded alone extended into the plastid genome through plastid-derived
  stretches (Ajuga and Alisma "mitochondria" were plastids; plastid reads are 10-30x deeper).

### 3.6 Downsampling (`recruit.rs`, final output block)

A read is kept when $u(\mathrm{name}, \mathrm{salt}) < \phi$, with $u$ the splitmix64 finaliser of the FxHash of the name mapped to $[0,1)$
(deterministic; with `mates_together` the first word of the name is hashed so mates stay together). $\phi$ is the smallest of
`target_depth/depth` (only with `--recruit-depth`, when the seed length is a genome size), `max_read_depth/d_read` (HiFi; `run` default 250;
$d_{read}$ = median sketch depth of the target's reads) and `max_bases/bases` (`run` default $10^9$, a memory bound: 4.3 Gb of Medicago
polymorpha plastid reads took 2.5 GB, 8.5 Gb of Carex echinata "plastid" reads 309 GB). The 250x and $10^9$ are empirical bounds. The input
cache (`recruit/input-cache`, a quarter byte per base plus names) is used when the output disk has room for twice the need
(`default_recruit_cache`).

### Recruitment status in the report

Each target in `recruit.json` has a `recruitment_state`.

- `converged` requires an evaluated pass with no bait gain and no unresolved evidence within the configured search. The check covers recruited reads and the partial candidates assigned to that target by the classifier.
- `unresolved` covers heuristic saturation, seeds supplied as reads, or evidence withheld during extension. For example, a read can contain a new k-mer that lacks independent support. The run retains this uncertainty when a later pass has no new evidence for that read.
- `truncated` covers stops due to the scan budget, the target size cap, or explicit bootstrap limits. A single scan classifies reads without extending the bait. The last allowed scan also skips extension, so a target can have the historical `bait_complete` reason and still be truncated.

The `stop_reason` field retains the historical stop label. The `converged` boolean is true only for the `converged` state. The state describes the observed bait search under the run's rules. Final FASTQ sampling is reported separately, and the state cannot establish physical genome completeness.

For a diagnostic run, add `--frontier-audit` to `ovasm find recruit`. The report links
to `source_candidates/frontier-evidence.jsonl`, with up to eight read examples per
target and reason. The report counts events, including repeated observations of the
same record. It reports some stages only as counts. Example coordinates refer to
the final bait and matching uses the read identifier and sequence. A match does not
establish the read's physical source. This option uses the existing scans and keeps
the recruitment and FASTQ rules.

### 3.7 Admit original reads with external source evidence

Use `--read-admission-manifest MANIFEST.json` with `ovasm find recruit` or with the
seed recruitment arm of `ovasm run hifi` and `ovasm run short`. Prepare the manifest
from independent source evidence for the input sample. The manifest has schema
`ovasm.read-admission.v1`, a nonempty `scope` and a `records` array. Each record has
the following fields.

| Field | Value |
| --- | --- |
| `target` | Exact seed target name, such as `mitochondrion` |
| `read_id_hex` | Lowercase hexadecimal bytes of the complete FASTA or FASTQ header after its initial marker |
| `sequence_sha256` | Lowercase SHA-256 of the parsed sequence bytes, preserving case and removing line breaks |
| `quality_sha256` | Lowercase SHA-256 of the FASTQ quality bytes, or JSON null for FASTA |
| `provenance` | Object with `path` and lowercase SHA-256 `sha256` for the evidence file. Relative paths resolve beside the canonical manifest. |

Each declared identifier must occur once in the manifest and once across the original
inputs. OVASM validates the evidence file and exact original record digests. The source
assessment itself comes from the evidence author. Record bytes come from the original
input stream. Existing passes establish the final bait call. New admissions are
appended in manifest order after ordinary downsampling, without extending the bait
or adding the other mate. An ordinary call for the requested target keeps its usual
sampling behavior. An ordinary call exclusively for another target takes precedence.
A partial record keeps its diagnostic entry even when it is explicitly admitted.
The same sampling hash applies to added records, and an addition must fit within
the remaining configured base cap.

The optional `read_admission` section of `recruit.json` records manifest and evidence
digests, the baseline classification and each record's outcome. `already_recruited`
describes the baseline call before sampling. `admitted` confirms the record was
written. The remaining outcomes are `competing_target`, `downsampled` and `base_cap`.
Final output counts include admitted records. Iteration counts, bait closure and
sketch-based depth describe the ordinary bait search.

### Exclude exact source-bound originals during seed recruitment

`--read-source-exclusions MANIFEST.json` applies to `find recruit` and
whole-genome `run` with `--recruit seeds`. The
`ovasm.read-source-exclusions.v1` manifest uses the same original identity
and provenance fields as read admission. Its target names the prohibited
destination. One original can exclude several targets. A positive admission
for the same original and destination conflicts with the exclusion.

The exclusion removes that target assignment before pass files, sampling
and read-derived bait growth. Partial reads can still be inspected in the
diagnostic files. A prohibited original cannot feed within-pass extension
or either leftover bootstrap for that target. Replay validates the same
identities and restores and checks original qualities for the final outputs.
Other target assignments remain eligible. No alternate target is inferred
from exclusion alone.

The starting seed files are supplied independently. The whole-library
k-mer walk cannot track an original's source, so it is skipped for protected
targets while read-level depth gates remain available. The recruitment
report records match passes, bindings and targets whose library walk was
skipped. Exclusion claims have the scope of the supplied evidence.

Generate this manifest from registered source junctions with:

```sh
python3 scripts/prepare_read_source_exclusions.py PLAN.json --out AUDIT.json \
  --out-exclusions EXCLUSIONS.json
```

The plan schema is `ovasm.read-source-exclusion-plan.v1`. Set `scope` and
`target`, list the original FASTA or four-line FASTQ inputs in `raw_reads`,
and provide `references` with an `id` and a physical `role` for each FASTA.
Inputs and provenance objects carry `path` and `sha256`. Each registered
window in `windows` gives its `id`, panel `reference` ID, FASTA `record`
ID, zero-based half-open `interval`, interior `boundary`, and `provenance`.
For a competing role, `exterior_side` marks the `left` or `right` side
of the boundary that identifies the competing physical source. Target windows
are veto contexts.

The exporter requires the competing window to occur exactly once in the
supplied linear reference records and its exterior sequence to be absent
from all other supplied roles. In the original reads, it accepts either
orientation. It abstains for mixed roles, ambiguous competing windows,
repeated matches of a competing window, and duplicate full read headers.
The audit records the checks and positions. The manifest binds every
exported record to that audit and preserves the original sequence and
quality hashes. A run with zero qualifying originals produces only the
audit. Establish the physical role and boundary independently of the
recruitment outcome. Use the same original inputs for recruitment.

### 3.8 Screen original reads at registered junctions

Run `python3 scripts/prepare_read_admission.py PLAN.json --out SCREEN_DIR` to prepare
a native read admission manifest. Use the same original four-line FASTQ inputs for
screening and recruitment. Plain and gzip inputs are supported. The plan uses schema
`ovasm.read-source-screen.v1`, a nonempty `scope`, and a `target` of `mitochondrion`
or `plastid`.

- `reads` lists original FASTQs. `panel` gives one FASTA for each of
  `mitochondrion`, `plastid`, and `nuclear`. A FASTA can contain several reference records.
  Each file uses an object with `path` and lowercase `sha256`.
- `windows` lists registered target junctions. Each object has an `id`, target FASTA
  `reference` identifier, and integers `start0`, `boundary0`, and `end0`. Coordinates
  are zero-based and half-open. Supply at least one base on either side of the boundary.
- Optional `qualification` binds an independent comparison file with schema
  `ovasm.read-source-qualification.v1`. That file declares `scope`, `target`,
  `panel_sha256` keyed by the same three roles, a `provenance` file binding, and
  `records`. Each record has `read_id_hex`, `sequence_sha256`, `quality_sha256`, the
  complete registered `window` object, and a `verdict` of `target_supported`,
  `competing`, or `unresolved`. The evidence author supplies the verdict using
  independent evidence for the same original record and interval.

Relative paths resolve beside the JSON file that declares them. The screen searches
all input reads for registered windows in either orientation, then compares the
complete exact read end and exact window against the named panel. Only a unique
target placement for both, at the registered coordinates and with an independently
qualified original record, is eligible. Exact comparison applies to literal panel
bases. Assembly gaps, unrepresented sequence and approximate matches remain outside
that comparison. A matching competing or unresolved declaration also withholds the
record. The report keeps up to two placements per role, enough to establish
ambiguity. Multiple windows on the same record yield at most one native declaration.

`source-screen.json` records the frozen inputs, exact placements and individual
decisions. If screening finishes with eligible records, `read-admission.json` is
ready for `--read-admission-manifest`. Reaching the candidate cap before all reads
are checked produces a truncated report and no admission manifest. This route
currently handles exact complete reads at registered windows.

For an additional short-read competition check, set `--phred-error-quantile Q`
with a probability strictly between zero and one. The script calculates each
read end's error-count quantile from its Phred+33 qualities, conditional on
independent calibrated errors. It searches both strands for equal-length placements
within that number of substitutions. Splitting the read into one more block than
the allowed substitutions guarantees that every such placement has an exact seed.
The search confirms each candidate by comparing the full read end.

A competing role, another position or orientation in the target, or an exhausted
search budget withholds an otherwise eligible record. The report records up to
two placements per role. `saturated` means at least two were found.
`--max-substitution-starts` sets the work budget per original read and role, with
500,000 starts by default. The independent qualification remains part of the
admission rule. This comparison covers literal equal-length sequence. Indels,
assembly omissions, sequence changes between individuals, and correlated quality
errors are outside its model.

## 4. Assembly (`assemble.rs`, `pipeline.rs::GraphAttempts::attempt`, `k_ladder`)

### 4.1 Closed-syncmer sparse graph

- Window $w = k - s + 1$ with $s = 31$ (default `--s`). Hash every s-mer (canonical). A k-mer starting at $i$ is a node iff the minimum-hash
  s-mer inside it is at its first or last position. Every window of $w$ k-mer starts holds at least one such k-mer, so consecutive sampled k-mers
  of a read are at most $k - s$ apart and overlap by at least $s$; density is about $2/w$ ($k = 1001$: $w = 971$, about 0.2% of positions).
- Nodes are keyed by the hash of the whole k-mer, so an error anywhere in it creates a different, rarely seen node. Edges join consecutive
  sampled k-mers of a read with their distance $d < k$; an edge and its reverse complement are one canonical edge. The share of error-free
  k-mers is $(1-e)^k$; at $e = 0.001$, $k = 1001$ it is 0.37, so node depth is about a third of read depth (doc, `run_auto_k`).
- Compacting non-branching runs gives unitigs; the GFA carries exact `(k - d)M` overlaps. With k shorter than the organelle's repeats each repeat
  collapses into a branching unitig, which evidence and linearization resolve with reads.
- Short reads use dense mode ($s = 0$, every k-mer a node) with Bloom filters so the first sighting of a k-mer or edge is not counted (an error makes
  $k$ new k-mers); gap rescue is off there (memory).

### 4.2 Node depth: count-weighted median

Over nodes with count $c \ge 3$, with $n_c$ nodes of count $c$:
$$D = \min\Bigl\{ c : \sum_{c' \le c} c'\,n_{c'} \ge \tfrac12 \sum_{c'} c'\,n_{c'} \Bigr\}.$$
This is the count at which half of all k-mer occurrences lie, the depth of the read mass. Excluding $c < 3$ removes the many error nodes.

### 4.3 Seed-anchored node depth (`depth_reference`, `assemble::run_engine`)

Problem: a deep repeat among the reads makes $D$ the repeat's depth. Eilema depressum: 3,861 for a repeat at ~140x read depth, threshold 386 left no
unitig. Platycheirus albimanus: 375 for a 14 kb element with 4 gene fragments, threshold 37 against a 13-gene mitochondrion at ~35x.

Rule. Sample the closed-syncmer k-mers (same $k$, $s$) of the sequences in `protein_seeds.<o>.depth.fasta`, the seeds of the protein-seed cluster with
the most genes. Let $R$ be the distinct nodes among them present in the graph with count $\ge 3$. If $|R| \ge$ `DEPTH_REFERENCE_MIN_NODES` = 50,
$D := \mathrm{median}\{\mathrm{count}(v) : v \in R\}$; otherwise $D$ stays the count-weighted median. The count $\ge 3$ restriction is needed because at
$k = 1001$ most of a read's k-mers hold an error ($0.999^{1001} = 0.37$ error-free) and each is a count-1 node: the median of all seed nodes was 1 for
Trichoderma longibrachiatum (threshold 1, 127 unitigs; Acleris sparsana killed at 30.6 GB). The minimum of 50 nodes is unstated. The reference is supplied
only with protein seeds.

### 4.4 Automatic minimum count

$$m = \max\bigl(\lfloor D/10 \rfloor,\ f(D)\bigr),\qquad
f(D) = \begin{cases} 1 & D < 6\ (\texttt{VERY\_LOW\_DEPTH\_NODES}) \\ 2 & 6 \le D < 8\ (\texttt{LOW\_DEPTH\_NODES}) \\ 3 & D \ge 8 .\end{cases}$$

`--min-count` overrides; `--min-count-fraction` $\varphi$ replaces $1/10$ ($\max(\lfloor \varphi D + 10^{-9} \rfloor, f)$). Nodes with count $\ge m$ are kept.
The fraction 1/10 has no stated derivation (`docs/rust-kernel.md`: "suits recruited organelle reads", not whole-genome input). The floors are empirical: Canis
lupus mitochondrion from 11 reads (~5x) gave one circular 16,720 bp unitig with threshold 1, two to four pieces of 9.8-14.9 kb with 2, and 14.2 of 16.7 kb with 3.

Edges (`build_graph`): an edge between kept nodes is kept if its count is $\ge m$, or $\ge$ `ARC_MIN_READS` = 3 and $\ge$ `ARC_FRACTION` $\cdot \min(n_a, n_b)$
with `ARC_FRACTION` = 0.35 (the arc rule of Syncasm in Oatk 1.0; Schoenoplectus repeat exits with 61-76 reads between nodes of ~130 were lost at an absolute
cutoff of 80).

### 4.5 Automatic k and the open-end rule (`run_auto_k_cleaning`, `ladder_choice`)

`K_LADDER` = 1001, 701, 501, 351, 251, 201, 151, 101. Assembly takes the largest k whose node depth reaches `AUTO_K_NODE_DEPTH` = 22, going down only as far as depth
demands. Evidence (empirical, Zou cohort, 12 accessions thinned to 15-30x): at 30x one circle in 2/12 at k = 1001 and 12/12 at 251; across 15-30x accessions came out
as one circle once node depth reached about 22; too short a k merges short repeats (k = 101 at 30x: 10/12). A shorter k sees more error-free k-mers but resolves fewer
repeats alone (Nipponbare at 501: 10-14 of 16 junctions; at 1001: 16).

Each k tried is scored by the tuple (open ends, unitigs, $-$total length, $-k$), smaller is better. The ladder is done when node depth reached the target and the graph
is closed or has more than `OPEN_END_MAX_UNITIGS` = 4 unitigs. A graph of at most 4 unitigs with an open end is read as a broken circle and the ladder continues. If no k
reached the target, or the chosen graph is small and open, the k with the minimal tuple is taken (`ladder_choice`). Examples: Rhizosphaera kalkhoffii (604 reads, ~106x):
k = 1001 lost a 10.5 kb stretch under the threshold and gave one open 33.8 kb unitig, k = 201 one circular 44,388 bp unitig; Platycnemis pennipes (~18x, no k at the
target): k = 101 three unitigs, k = 1001 a closed two-unitig graph of 17.2 kb.

Second ladder. If the best result is at most 4 unitigs and still open and `--min-count` was not given, the ladder is run again with node threshold
`OPEN_END_MIN_COUNT` = 3 (a stretch present in a fraction of copies, such as an optional intron, stays under $D/10$ at every k) and kept only if it gives fewer open ends.
Evidence is one species (Rhizosphaera).

Correction fallback. If even k = 101 stays below 22 (about 20x and less), reads are self-corrected (`assemble correct`, k 17 then 25, self rounds 15, 21, 25) and the ladder
is climbed again; the corrected graph is used only if it keeps at least half the bases of the uncorrected one (Canis lupus, 11 reads, 75 kb: 5 reads, 1.4 kb after
correction). Thinned Zou accessions: one circle within 10 bp in 9/12 at 15x and 12/12 at 20x with correction against 8 and 11 without; at 30x, 12/12 without and 11 with.

### 4.6 Starting k in `run` and the pipeline-level ladder (`pipeline.rs::start_k`, `k_ladder`)

If the recruited reads' kept depth (their sketch median times the kept share) is at least `START_K_READ_DEPTH` = 3 x 22 = 66, the first attempt is made at k = 1001 directly
instead of by node depth (Medicago polymorpha, ~250x kept: node depth 16 at k = 1001 because of tandem-repeat nodes, k fell to 251 and the graph into 1,974 unitigs; at 1001
one circle of 287,741 bp).

This is a second, independent ladder over whole attempts. `GraphAttempts::attempt(k)` runs assemble, the optional component filter, evidence and linearize and returns an
`Outcome`. An attempt is accepted when it has at least one molecule and
$$\mathrm{explained} = \frac{\sum \text{molecule bases}}{\text{graph bases}} \ge \texttt{MIN\_EXPLAINED} = 0.5 .$$
If not accepted, the next smaller k of `K_LADDER` is tried; the first accepted attempt is the result; if none, the attempt with the largest explained share (earliest on ties) is kept
and `k_accepted = false`. Evidence: Carex laevigata mitochondrion (2.77 Mb): k = 1001 to 351 gave no molecule, k = 251 two molecules of 11,560 and 349 bp sharing nothing with the
published genome. The 0.5 rests on that case. If the accepted attempt came from a given start k and is a small open graph (at most 4 unitigs, an end without a link), one more
attempt at the automatic k is made and kept if accepted with fewer open ends (on small animal/fungal mitochondria a fixed 1001 left broken circles: Rhizosphaera 0.76 of
reference at 1001, one 44,388 bp circle at 201).

### 4.7 Cleaning (`assemble::clean`, up to 50 rounds)

Rules act on unitigs of depth $dp$, with graph depth $D$. An error touches about $2k$ k-mers and syncmer sampling adds up to $k - s$ per side, so error shapes span up to ~$4k$.

- Reference depth floor: $0.5D$.
- Tips: length $< 4k$ and $dp < 0.5\times$ the strongest way on from its junction, unless it is the strongest at some end. The deepest way is never removed (ties: more
  k-mer observations, then lower index), because at systematic-error sites (ONT homopolymers) the true branch can hold barely half the reads.
- Bubbles: a detour rejoining within length $+\,2k$ is removed when $dp < 0.5\times$ the surrounding depth, or when it is shorter than $2k$ (an arm that short differs by few bases,
  e.g. a homopolymer that splits reads half and half) and is not the deepest way somewhere.
- Isolated unitigs: $dp < 0.25D$, or length $< 2k$ (those are kept in the first, deferred round so gap rescue sees them).
- Components whose deepest unitig is below $0.3D$ are removed; so are components with no unitig of at least $\min(4k, 1000)$ and no circle (tangles of nuclear repeat; Nipponbare
  short reads: 25 such components, none mitochondrial).
- Near copies: an arm shallower than half the deeper way beside it is removed, whatever its length, when at least `VARIANT_SHARE` = 0.8 of its 15-mers (`VARIANT_W`) occur in the deeper
  way (one difference per 100 bp keeps ~85%; Col-0 chromosome 2 NUMT, ~99.9% identical).
- NUMT flanks: with `nuclear_share` $= D_n/\text{core depth}$ from recruitment, a dead-end branch off a junction where a deeper way continues, no deeper than
  `NUMT_BRANCH_FOLD` (1.5) $\times$ nuclear_share $\times D$, is removed whatever its length (Medicago truncatula A17: 623x, nucleus 84x, a 4,732 bp dead end at 12% depth).
- Isolated circles shorter than k are tandem arrays (a 62 bp satellite in Altai-5 and Est-1) and are removed.

The constants 0.5, 0.25, 0.3, $2k$ and $4k$ have stated purposes but no stated validation; 0.8/15 and 1.5 rest on one or two named samples.

### 4.8 Gap rescue and fill (`rescue_gaps`, `fill_gaps`)

- Rescue: for a dead exit that some read leaves by and that the read later reaches a dead entry from (within `MAX_SPAN` = 30 kb, not passing through the live graph), the most frequent
  k-mer run between them is added if at least two reads show it exactly. Long reads: a deferred first round plus two more; dense mode: two rounds. Aim: systematic errors (a homopolymer
  called differently by most reads) that keep a k-mer under the count threshold while a group of reads still agrees. "At least two" is unstated.
- Fill (`fill_gaps`): a local assembly at `GAP_K` = 31 from reads selected by the 150 bp (`GAP_WINDOW`) at each end; at most `GAP_MAX` = 600 bp in dense mode, `LONG_GAP_MAX` = 10,000 bp
  between long-read dead ends (within one read's span). None of these numbers has stated evidence.
- Then boundaries are normalised and links blunted.

## 5. Component filter (`components.rs`, `pipeline.rs::ComponentFilter`)

Applied when reads came partly from depth clusters (`--recruit discover|both`) and, in gene-free mode, to the mitochondrion when protein-seed passes exceeded one.

For each connected component: length, length-weighted mean depth $d$, closed (every segment end has a link), and distinct panel genes of the target organelle $G_t$ and the other $G_o$
(assembled sequence is scanned like reads but with `min_gene_segments = 1`: one stretch is enough). The reference component has the most target genes (with $G_t > G_o$ and a depth), the
deeper on ties. With reference depth $D_r$, `decide` tests in order:

1. reference: kept; no depth: kept.
2. `reject_other_only` (multichrom primary profile only): no target gene and some other gene: removed.
3. `reject_beside_closed_reference`: reference closed and all of the component's target genes (by name) already on the reference, or none: removed. Enabled only for the mitochondrion
   when protein-seed passes > 1. Evidence (PMAT2 benchmark, 47 graphs): Aporia crataegi, a 98 kb component at 25x with nad2/nad4l beside a 15.1 kb circle at 34x; Nemoura dubitans, a 51 kb
   gene-free component at 100-167x beside a 15.7 kb circle at 70x; none of the removed sequence on the published genome. A component with a gene the reference lacks stays (circles of a
   multi-chromosome mitochondrion carry different genes); so does everything beside an open reference.
4. other organelle: $G_o \ge$ `min_other_genes` (5) and $G_o \ge 2G_t$: removed.
5. $G_t > 0$: kept unless $D_r > 2.5\,d$ (`gene_depth_factor`; organelle-derived nuclear sequence sat at 0.10-0.28 of the reference depth in Climacium and Luzula, published components at 0.51 and above).
6. $G_t = 0$: kept if $d$ is within a factor 2 (`depth_factor`) of $D_r$ either way. Gene-free components absent from the published genome sat at 0.09-0.45 and 2.5-3.3 times the reference;
   no published component was gene-free in these graphs.

The module doc states that the factors were set on graphs whose published genomes had been seen and are not validated on unseen data. Gene calling uses $k = 6$, `min_orf` 80 and `min_hits` 5
(6 for fungi in gene-free mode only; the discover-path filter always uses 5). Optional original-read context can exclude an unmarked component when every associated read favours nuclear
context and covers its full sequence.

## 6. Evidence, linearization, chromosomes, assessment

### 6.1 Read evidence on junctions (`evidence.rs`)

Each GFA link becomes a junction sequence: `flank` = 2,000 bp from the end of the left segment plus the start of the right. A read supports a link when its exact 21-mer hits on that junction
are collinear (diagonals within `diag_tol` = 30) and anchor at least `min_anchor` = 500 bp on both sides of the joint (capped at 90% of a short segment), with at least `min_side_hits` = 10
hits per side. Short reads: anchor 40, side hits 5, and exact overlap-edge words counted separately (`short_read_adjacency`). Competitive assignment: a read counts for the junction it fits
best, with ties settled by the read's next or previous crossing (`assign`); a witness needs `JUNCTION_MIN_SHARE` = 0.5 (Schoenoplectus plastid: 21 sub-k inverted repeats gained 8-512
alternative reads each; one case). The preset numbers cite `_EVIDENCE_PRESETS` of a Python file that is not in this tree.

Repeat segments (`repeat_segments`): a segment is a repeat when at least two distinct read-supported links leave an end and, if depth is recorded, its median depth (`dm`) is at least `REPEAT_FOLD`
= 1.35 times the graph's length-weighted median (Nipponbare gold graph: repeats 1.56-2.35x, unique at most 1.19x). Further rules: crossroads (two ways at each end, the second-strongest holding at least
`CROSSROAD_SHARE` = 1/3 of the strongest's reads) are repeats unless depth-single-copy (at most 1.2x and at least 5,000 bp: Medicago arabica 69.7 kb circle at 0.92x); ambiguous junctions at
$\ge$ 1.1x; hidden repeats crossed in two contexts of at least 3 reads each between segments at $\ge$ 0.5 of the median. The constants 1/3, 1.1, 0.8 (`HIDDEN_DOMINANT`) and 0.5 (`MAJOR_DEPTH`)
carry no stated evidence.

Traversals: from a read's supported links, consecutive joints spaced by a segment's length (tolerance 5% of the length, floor 50 bp) form X $\to$ S $\to$ Y traversals; at a repeat with two entries
and two exits they give the two pairings. Minority share $\hat p = m/(m+M)$ with a percentile bootstrap over 1,000 resamples (`minority_share`: draw $n = m+M$ reads each minority with probability
$\hat p$, take the 2.5% and 97.5% quantiles). Every read is also recorded as a single-record thread (`thread_evidence`), used by the linearization read model.

### 6.2 Linearization (`linearize.rs`, `read_model.rs`)

- Anchors are the non-repeat segments. Two anchors are bridged by a read-supported phased path through repeats (reads $\ge$ `min_reads` = 2, set by `run`), by a supported direct link, or
  topologically (`TOPO_READS` = 1, at most `TOPO_MAX_REPEATS` = 8 repeats) at an anchor end that fewer than `ABSENCE_READS` = 10 reads cross. Topological and weak bridges make a molecule non-decisive.
- Connected components of the anchor graph are solved separately. Within one, a cover is a set of molecules visiting each anchor once through bridges (one molecule is a Hamiltonian path in the small
  anchor graph); repeats are laid down as often as the chosen bridges cross them; a molecule is circular when a bridge returns from its last anchor to its first. Covers of every size compete together.
- Once the required-anchor covers have been enumerated, the remaining search budget can extend an open cover to an observed record end. The entire inbound thread must match the existing oriented
  path. Its tail must follow existing bridges to unused optional anchors and respect phase constraints. These additional open candidates carry the same representation and topology limits as the original
  covers. On exhaustion, the original covers remain available and the search is marked incomplete.
- Ranking when thread evidence is present (always, since `check evidence` writes it): `read_model::ReadFit::order`, lexicographic and not a likelihood: more joint records explained, then more
  thread records, then more adjacency-record occurrences, then fewer unobserved junctions, then fewer weak junctions; then fewer copy-number deviations, then fewer molecules (a representation
  tie-break). Search work is bounded by `WORK_BUDGET` = $10^7$; on exhaustion `ranking_model = "legacy_after_read_budget_exhaustion"`.
- Legacy ranking (the module doc's description; `ranking_model = "legacy_bridge_copy"`): fewer linear molecules, then complete junction evidence, then the higher
  $$\mathrm{score}_{\mathrm{ranked}} = \sum_i \ln \max(f_i, 10^{-9}) \;-\; 2\,\Delta_{\mathrm{copy}},$$
  where $f_i$ is the fraction of reads leaving an anchor end that take bridge $i$ (geometric mean of the two directions) and
  $\Delta_{\mathrm{copy}} = \sum_{\mathrm{repeats}} |\mathrm{passes} - \max(1, \mathrm{round}(\mathrm{depth}/\mathrm{median}))|$ over repeats the cover reaches (`COPY_PENALTY_LN` = 2, a factor $e^{-2}$ per copy; unstated justification).
- Decisive: the best path is at least $\ln 2$ more log-likely than the runner-up (`DECISIVE_LN`).
- Tails: a best path whose end leads back into its own interior is cut into a circle and a tail, only if the closing bridge is fully read-supported and the tail is at most `TAIL_MAX_FRACTION` =
  0.1 of the circle's anchor sequence (Nipponbare's 376 kb path would fall into 190 + 186 kb from ONT reads and 226 + 150 kb from short reads, so it stays whole).
- Unbridged anchors of at least `MIN_UNBRIDGED_BP` = 1,000 are linear molecules; shorter ones and components without a cover go to `partial_paths` and the component stays unsolved. Shallow anchors
  (judged by depth when recorded) are optional low-frequency configurations.

### 6.3 Unified graph (`unify.rs`)

`assemble unify` writes `organelle.gfa`: assembly segments plus up to `--max-alternatives 3` read-supported alternative connections, and the replicate tags of `--mode standard`. It is a re-encoding, not a separate estimator.

### 6.4 Chromosomes (`chromosomes.rs`)

`linearize` writes the fewest molecules that explain the graph, which is not a chromosome count. Rules in order:

1. Sequences joined by a link with at least `MIN_LINK_READS` = 2 reads are one chromosome.
2. Molecules of one chromosome are its conformations; each has a share (the reads of the pairing it takes at the least supported repeat it crosses).
3. Two molecules are different chromosomes when no supported link joins them, or when every read across a shared repeat takes the separating pairing and the missing pairing's share is below
   `MAX_MISSING_SHARE` = 0.05 at 95%: with $n$ reads all of one pairing the 95% upper bound of the missing share is $1 - 0.05^{1/n}$, below 0.05 from $n = 59$ (**derived**).
4. A repeat no read crosses, or too few, leaves a range of counts.
5. A component below `SHALLOW` = 0.5 of the main one's depth, with fewer than `MIN_SEED_BP` = 300 bp of the organelle's seed 21-mer blocks, is possibly nuclear or contamination and is not counted.
6. Any other unjoined component at `SOURCE_FOLD` = 2 times the main depth or more, or at its reciprocal or less, is "source pending": not counted, only in the range's maximum.

#### Audit read length for a repeat route

`python3 scripts/audit_repeat_span_opportunity.py PLAN.json --out AUDIT.json`
counts original reads long enough to contain the literal spelling of a
declared repeat route. Use schema `ovasm.repeat-span-opportunity-plan.v1`
with `scope`, one `graph`, a list of `raw_reads`, and a list of `routes`.
The graph and each FASTA or FASTQ input use `path` and `sha256`. Each route
has an `id` and a `path` of oriented GFA segment names, with a flanking
segment at each end and at least one interior segment.

The minimum length is the spelled interior plus one base from each
exterior overhang. Internal overlaps are counted once and repeated visits
are spelled separately. The graph must supply unambiguous exact overlaps.
The report gives the number of original records meeting the minimum per
route. This length condition is necessary for a literal, gap-free span.
Read alignment and physical pairing require separate evidence.

### 6.5 Assessment (`assess.rs`)

Derived from the run's own files; the reads are not re-read.

- Single-copy depth: weighted median over the bases of segments the molecules walk exactly once; segment depth is expressed in these units.
- Copy number: the walk count and the depth in single-copy units agree when within `COPY_TOLERANCE` = 0.5; the report gives the share of bases where they agree.
- Shallow sequence: depth below `SHALLOW` = 0.25 of single-copy (a different constant from the chromosome rule's 0.5 and the cleaning rule's 0.3).
- Junctions: reads behind every link used, the five weakest, any with none.
- Complexity: simple (one segment), complex (more), multipartite (more than one component); dead ends and double forks (segments with at least two links at both ends).
- Genes: panel genes on each molecule and the call of `call_from`.
- Plastid: quadripartite structure when the circular molecule walks one segment twice in opposite orientations in two runs of equal length separating two single-copy runs.

## 7. Short reads and noisy long reads (brief)

Short reads (`run short`; `assemble::dense_k`, `pipeline::sr_ladder`, `close.rs`). Recruitment uses the k-mer gate (section 3.3) and keeps mates together. The dense graph uses
$k = $ median read length $- 23$ (first 10,000 reads of the first file), made odd and at least 21, so that each read gives 24 k-mers (Nipponbare 150 bp: 127, Salvia 100 bp: 77).
`sr_ladder` tries `SR_K_LADDER` = 99, 77, 55 below that until one attempt gives a single molecule; the result is the accepted attempt with the fewest pieces (molecules of at least `MIN_PIECE`
= 1,000), then the most molecule bases, then the larger k. `close.rs` closes a linear molecule into a circle with read pairs: at least `MIN_JOINING_PAIRS` = 5 pairs joining end to start; extension by
majority with at least 3 local reads and 80% dominance (contexts 90, 70, 50, 31), up to `MAX_FILL` = 2,000 bp; if reads do not reach, a reference walk may fill up to 5,000 bp only when its length
agrees with the pair distance within max(100 bp, 30%), and those bases are labelled the reference's. Evidence is one sample (ZM4 plastid, hairpin region).

Noisy long reads (`run ont|hybrid`, `longread.rs`). Noisy reads cannot extend a bait (error k-mers are short random words that the nuclear genome also holds), so: (1) one recruitment pass with
the seeds at k = 21 (`RECRUIT_K`; at 15 the library covered every read by chance, 61 Gb of "plastid" reads out of 131 Gb for ZM4); (2) self-correction of each organelle's reads (k 17 then 25);
(3) further rounds in which the round's reads and the partly covered ones (a 2 kb bait block) join the bait, a round adding less than `ROUND_GROWTH` = 5% of the bases ends the rounds, partly
covered reads join only at depth $\ge$ `PARTIAL_DEPTH_SHARE` = 0.25 of the target's recruited reads (or the 0.9 quantile of their own depths, `SELF_DEPTH_QUANTILE`, when the target recruited nothing),
reads deeper than `CORRECT_MAX_DEPTH` = 300 are not corrected, microsatellite runs of at least 21 bases with period at most 6 are masked (`TANDEM_*`); (4) hybrid only: the corrected long reads seed one
pass over the short reads, then the long reads are corrected again with them; (5) the corrected reads enter the HiFi stages as target reads at k 201 unless the graph at k 101 holds at least 1.10 times
its sequence (`ASSEMBLE_K_*`; Nipponbare ONT 273 kb at 201 and 405 kb at 101 of 376 kb). These numbers rest on two ONT libraries (Nipponbare, ZM4); CLR is untested.

## 8. Outputs

Per organelle directory (or `<out>` for one organelle): `molecules.fasta` (headers carry `length=`, `circular=`, `path=`), `organelle.gfa`, `assembly.gfa` (with a filter also
`assembly.filtered.gfa` and `components.json`), `assembly.json`, `evidence.json`, `linearization.json`, `chromosomes.json`, `assessment.json`, `summary.json` (k, `k_accepted`, the attempts of the
ladder, recruitment and seed records, warnings, replicate record), `recruit/recruit.json`, and `protein-seeds/protein_seeds.json` when protein seeds were used.

`--mode standard` adds downsampled replicates. The replicate depth is half the full depth, at most `REPLICATE_DEPTH_CAP` = 100x and at least `REPLICATE_MIN_DEPTH` = 30x (so full depth 60x; automatic
k gave one circle in 12/12 thinned Zou accessions from 30x), and a replicate holds at most `REPLICATE_MAX_FRACTION` = 0.5 of the reads. Three replicates are drawn, then more until the junction set
was the same in three in a row or `--max-replicates`. Each junction of the full-depth graph gets the share of replicates that reproduce it (`rs:f`, `rn:i` on `organelle.gfa` links,
`replicates/junctions.tsv`). Replicates assemble at the same k as the representative.

### 8.1 `result.json` (`result.rs`, schema `ovasm.result.v1`)

Written whether the run worked or not. Fields: `schema`, `version`, `command` (argv), `status`, `error`, `organelles[]`, `seconds`, `peak_memory_mb`. Each organelle entry has `organelle`, `status`,
`error`, `molecules[]` (`name`, `length`, `circular`, parsed from the `molecules.fasta` headers), `total_length`, `k`, `gene_calls` (from `assessment`), `seconds`, `files` (paths relative to `<out>`;
only those that exist among `molecules.fasta`, `organelle.gfa`, `assembly.gfa`, `summary.json`, `assessment.json`, `evidence.json`) and `warnings` (summary warnings, recruitment warnings, and "no k
explained half of the graph" when `k_accepted` is false).

| organelle status | rule |
|---|---|
| `failed` | the organelle's step stopped with an error |
| `complete` | molecules exist and all are circular, or exactly one molecule and the linearization is `decisive` |
| `partial` | molecules exist and some are open |
| `no_molecule` | a graph but no molecule |

The run status is the worst organelle status in the order complete < partial < no_molecule < failed, and `failed` when the run itself errored. `complete` states topology and decisiveness only: it does
not encode `k_accepted`, gene completeness or the chromosome call (these are in `warnings`, `gene_calls` and `chromosomes.json`), so a circular result below the 0.5 explained-share rule is still
`complete` with a warning.

## 9. Provenance of the numbers

- **Derived:** `min_orf` stretch rate (~2%), panel 6-mer share (~1%), log-bin width 0.25 (19%), the 20 Mb minimum nuclear genome and the size caps (genome-size facts), closed-syncmer density
  $2/(k-s+1)$, error-free fraction $(1-e)^k$, `MAX_MISSING_SHARE` (binomial bound), the Poisson argument for $0.75c$, the dense-mode k from read length.
- **Empirical, set on named samples:** every other threshold that carries an evidence note above. Most are fitted on one to four species. The component filter's module doc states that its factors
  are not validated on unseen data; the noisy-read constants rest on two ONT libraries.
- **Unstated:** the depth-gate fold and uniformity numbers, `min_cluster_reads`, `min_genes`, flank, sample size, `min_cover` and the other preset values, `extend_depth`, `extend_min_count`, `DEPTH_REFERENCE_MIN_NODES`,
  the node-count fraction 1/10, most cleaning fractions, gap constants, evidence anchor/flank numbers, `COPY_PENALTY_LN`, `MAJOR_DEPTH`, `CROSSROAD_SHARE`, `AMBIGUOUS_FOLD`, `HIDDEN_DOMINANT`, `MIN_SEED_BP`.
