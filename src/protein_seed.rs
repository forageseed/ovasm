//! Protein seeds: the sample's own organelle reads as recruitment seeds, found by depth and
//! conserved proteins instead of nucleotide similarity to a seed library (`ovasm find
//! protein-seeds`, `ovasm run hifi --protein-seeds`).
//!
//! Nucleotide seeds start recruitment only where the sample shares 21-mers with a library genome.
//! With the land-plant library's two bryophytes left out, its seeds recruited none of the
//! mitochondrial reads of three of five mosses even after 60 passes, while the same mosses'
//! organelle proteins were found at once (leave-out test, 2026-11-05). Protein sequence stays
//! recognisable far beyond nucleotide k-mers, so the seeds are taken from the reads themselves:
//!
//! 1. Depth: only organelle-like reads (`DepthGate`: deep and uniformly deep along the read, the
//!    test of `find discover`) are looked at, so nuclear copies of organelle genes (NUMTs and
//!    NUPTs, at nuclear depth) never start a seed.
//! 2. Clusters: those reads are split by their median depth into depth peaks (`find discover`'s
//!    clusters; a plastid and a mitochondrion form separate peaks).
//! 3. Proteins: up to `sample` reads of each cluster are translated in six frames and their
//!    stop-free stretches assigned to the protein panel's genes (`ovasm find classify`).
//! 4. Reads: a read is of a class when at least `min_read_stretches` of its gene stretches are of
//!    that class and more than of the other; a single stretch is too weak a sign.
//! 5. Call: a class counts in a cluster only when the cluster's reads of that class hold many
//!    distinct genes at many stretches per scanned read (`min_genes`, `min_density`, as in
//!    `panel::call_from`); a deep repeat cluster with a few gene-bearing reads seeds nothing
//!    (Polytrichum commune: a 6,317x cluster with 24 mitochondrial hits in 3,000 reads). Both
//!    classes can count in one cluster where parts of both genomes share a depth (Ulota crispa:
//!    much of the mitochondrion at the plastid's 1,878x).
//! 6. Seeds: the counting classes' reads, cropped to the span of their class's stretches plus
//!    `flank` bases on each side.
//!
//! The same reads reversed (not complemented) are scanned as a control of chance hits: a gene
//! counts with its stretches in excess of the control's (`excess_over_control`), and a class that
//! counts in the control seeds nothing there (`class_seeds`). An organelle
//! without a called cluster gets no protein seeds; `run` then seeds it from the library.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::Result;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::Serialize;

use crate::panel::{self, ClassReport, Panel, PanelParams};
use crate::recruit::{peak_rss_mb, read_batches};
use crate::seedfree::{depth_peaks, DepthGate, GateReport};

pub const ORGANELLES: [&str; 2] = ["mitochondrion", "plastid"];

#[derive(Debug, Clone, Serialize)]
pub struct ProteinSeedParams {
    /// Reads scanned per depth cluster at most (evenly spread over the cluster).
    pub sample: usize,
    /// Bases kept on each side of a seed read's gene stretches.
    pub flank: usize,
    /// Reads a depth peak needs to become a cluster.
    pub min_cluster_reads: usize,
    /// Gene stretches of one class a read needs to be of that class (and more than of the other):
    /// a single stretch is too weak (Polytrichum commune: 84 of 1,452 such mitochondrial seeds were
    /// on no organelle).
    pub min_read_stretches: usize,
    /// The depth gate's nuclear peak, given instead of judged (`find recruit --nuclear-depth`).
    pub nuclear_depth: Option<u32>,
    /// The nuclear genome's collapsed size: the nuclear peak is the sketched bases over it
    /// (`find recruit --genome-size`), when no `nuclear_depth` is given.
    pub genome_size: Option<u64>,
    /// The depth gate's organelle threshold in multiples of the nuclear peak (default 3).
    pub depth_fold: Option<f64>,
    pub panel: PanelParams,
}

impl Default for ProteinSeedParams {
    fn default() -> Self {
        Self {
            sample: 3000,
            flank: 2000,
            min_cluster_reads: 20,
            min_read_stretches: 2,
            nuclear_depth: None,
            genome_size: None,
            depth_fold: None,
            panel: PanelParams::default(),
        }
    }
}

#[derive(Serialize)]
pub struct PanelGeneSet {
    pub organelle: &'static str,
    pub family_count: usize,
    pub families: Vec<String>,
}

#[derive(Serialize)]
pub struct ClusterGeneCount {
    pub organelle: &'static str,
    pub family: String,
    pub forward_stretches: u64,
    pub reversed_stretches: u64,
    pub forward_records: u64,
    pub reversed_records: u64,
    /// Input records contributing this family in both orientations.
    pub both_orientation_records: u64,
}

#[derive(Serialize)]
pub struct ClusterCall {
    /// Centre of the depth peak (sketch depth).
    pub depth: f64,
    /// Organelle-like reads nearest this peak.
    pub reads: u64,
    /// The sequence those reads stand for: their bases over the depth (`MIN_CLUSTER_SPAN`).
    pub span_bp: u64,
    pub scanned: u64,
    /// Each class's genes by their stretches in excess of the reversed control's
    /// (`excess_over_control`); the densities and the call are of these.
    pub mitochondrion: ClassReport,
    pub plastid: ClassReport,
    pub mitochondrion_density: f64,
    pub plastid_density: f64,
    /// "mitochondrion", "plastid", "mixed" or "none" (`panel::call_from`).
    pub call: &'static str,
    /// The classes that count in this cluster and seed from its reads of that class.
    pub seeded: Vec<&'static str>,
    /// The call of the same reads reversed, not complemented: "none" unless chance or
    /// low-complexity sequence fakes genes.
    pub reversed_control_call: &'static str,
    pub reversed_control_genes: (usize, usize),
    /// The reversed control's stretches per gene (both classes, most first).
    pub reversed_control_segments: Vec<(String, u64)>,
    /// Per-family counts after assigning each record a class, before subtraction or filtering.
    pub gene_evidence: Vec<ClusterGeneCount>,
    /// Seed reads this cluster gave.
    pub seed_reads: usize,
}

#[derive(Serialize)]
pub struct ProteinSeedReport {
    pub tool: &'static str,
    pub version: &'static str,
    pub params: ProteinSeedParams,
    pub panel_proteins: usize,
    pub panel_genes: usize,
    /// Effective family identities after panel exclusions, separated by organelle.
    pub panel_gene_sets: Vec<PanelGeneSet>,
    /// Records are FASTX input records, with no claim of original-molecule deduplication.
    pub gene_evidence_record_unit: &'static str,
    /// The audit covers the existing depth-cluster scan, not the whole-read fallback.
    pub gene_evidence_scope: &'static str,
    pub reverse_control: &'static str,
    /// The alignment E-value a gene stretch needs (`panel::Panel::max_evalue`); null: shared
    /// k-mers alone decide.
    pub panel_max_evalue: Option<f64>,
    pub gate: GateReport,
    pub reads: u64,
    pub organelle_like_reads: u64,
    pub clusters: Vec<ClusterCall>,
    /// Per organelle: seed reads, their bases and the FASTA (none without seeds).
    pub seed_reads: FxHashMap<&'static str, usize>,
    pub seed_bases: FxHashMap<&'static str, u64>,
    pub seeds: FxHashMap<&'static str, PathBuf>,
    /// The scan of every read for the organelles no depth cluster seeded (`whole_read_seeds`).
    pub whole_read_scan: Option<WholeReadScan>,
    pub elapsed_seconds: f64,
    pub peak_rss_mb: Option<f64>,
}

/// Seeds found without the depth gate: every read scanned, a read a seed when it holds at least
/// `min_read_stretches` stretches of at least `WHOLE_READ_MIN_GENES` distinct genes of the class.
#[derive(Serialize)]
pub struct WholeReadScan {
    /// The organelles it looked for (those without depth-cluster seeds).
    pub organelles: Vec<&'static str>,
    pub reads: u64,
    pub bases: u64,
    /// Seed reads per organelle (cropped to their genes and flanks, as cluster seeds).
    pub seed_reads: FxHashMap<&'static str, usize>,
    /// Of those, the reads the depth gate calls organelle-like (sketch depth above its threshold).
    /// The scan runs when no depth cluster was seeded, which happens for two reasons: the organelle
    /// is no deeper than the nucleus, or it is deep but its clusters' sampled reads called no class
    /// (Acanthosoma haemorrhoidale: mitochondrion at 312x over a 4x nucleus, clusters called none).
    /// Only the first is shallow.
    pub seed_reads_above_gate: FxHashMap<&'static str, usize>,
    /// The depth gate's nuclear peak stands for less sequence than a nuclear genome
    /// (`seedfree::MIN_NUCLEAR_GENOME`): the nucleus does not show in the sketch, and the gate's
    /// threshold, a multiple of that peak, sets nothing apart.
    pub gate_without_nucleus: bool,
    pub seconds: f64,
}

impl WholeReadScan {
    /// Some organelle has seeds, and at most half of them pass the depth gate or the gate has no
    /// nucleus to judge them against. Gliophorus irrigatus (organelle ~10x; nuclear peak 2x
    /// standing for 0.8 Mb): 13 of 21 cropped seeds passed a 6x threshold, recruitment was
    /// gated over the usual passes and took in nuclear repeats above it (59 molecules, 927 kb,
    /// none of the 46.5 kb mitochondrion); judged shallow, 27 kb of it at 99.7% identity.
    fn shallow(&self) -> bool {
        self.seed_reads.iter().any(|(o, &n)| {
            n > 0
                && (self.gate_without_nucleus
                    || 2 * self.seed_reads_above_gate.get(o).copied().unwrap_or(0) <= n)
        })
    }
}

/// Whether a whole-read seed, the span `lo..hi` of `seq` (`crop`), is organelle-like by the depth
/// gate: judged on the seed itself, since a read that is partly nuclear is not organelle-like
/// along its whole length however deep its organelle part (Gliophorus irrigatus: 6 of 21 whole
/// reads above the gate, the organelle 10x over a 2x nucleus; called shallow, recruited ungated
/// in one pass, 7 reads from the cropped seeds).
fn seed_above_gate(gate: &DepthGate, seq: &[u8], lo: usize, hi: usize) -> bool {
    gate.profile(&seq[lo..hi]).is_some_and(|d| d.organelle_like)
}

/// Distinct genes of its class a read needs in the whole-read scan. Without the depth gate a
/// nuclear copy of organelle DNA (NUMT, NUPT) is not set apart by depth, and most are a few kb of
/// one or two genes; an animal mitochondrial HiFi read of 10 kb holds five or more.
pub const WHOLE_READ_MIN_GENES: usize = 3;
/// Shared k-mers each counted stretch needs with its gene in the whole-read scan. A chance stretch
/// of a nuclear read reaches 5-8 (transposon ORFs, low-complexity peptides), a mitochondrial gene
/// 13-140 on its best stretch (references of the PMAT2 benchmark scanned with the metazoa panel);
/// with 2 genes at 5 hits, 1,514 of 74,216 Gibbula magus reads passed, of 19 mitochondrial ones.
pub const WHOLE_READ_MIN_HITS: u32 = 12;
/// Seed reads kept per organelle in the whole-read scan (in input order).
pub const WHOLE_READ_MAX_SEEDS: usize = 3000;

impl ProteinSeedReport {
    /// Whether an organelle is no deeper than the nucleus in these reads, so recruitment should not
    /// gate on depth either: its seeds came from the whole-read scan and at most half of them pass
    /// the depth gate (the gate's own test, applied to the organelle's own reads).
    pub fn shallow(&self) -> bool {
        self.whole_read_scan.as_ref().is_some_and(WholeReadScan::shallow)
    }

    /// One line for the terminal.
    pub fn summary(&self) -> String {
        let calls: Vec<String> = self
            .clusters
            .iter()
            .map(|c| format!("{:.0}x {} ({} reads, seeds {})", c.depth, c.call, c.reads, c.seeded.join("+")))
            .collect();
        let whole = match &self.whole_read_scan {
            Some(w) => format!(
                "; whole-read scan for {} ({} reads): {} seed reads, {} above the depth gate{}",
                w.organelles.join(" and "),
                w.reads,
                w.seed_reads.values().sum::<usize>(),
                w.seed_reads_above_gate.values().sum::<usize>(),
                if self.shallow() { " (shallow)" } else { "" }
            ),
            None => String::new(),
        };
        format!(
            "[ovasm] protein seeds: {} organelle-like reads in {} clusters [{}]{}; seeds: mitochondrion {}, plastid {} reads ({:.0} s)",
            self.organelle_like_reads,
            self.clusters.len(),
            calls.join(", "),
            whole,
            self.seed_reads.get("mitochondrion").copied().unwrap_or(0),
            self.seed_reads.get("plastid").copied().unwrap_or(0),
            self.elapsed_seconds
        )
    }
}

/// One scanned read: its cluster, its gene stretches, the genes of its reversed control, and its
/// id and sequence when it has stretches (a seed candidate).
struct Scanned {
    cluster: usize,
    stretches: Vec<(u16, usize, usize)>,
    control: Vec<u16>,
    read: Option<(Vec<u8>, Vec<u8>)>,
}

/// The cluster of each organelle-like read: the depth peak nearest its median (log2 units).
fn nearest(peaks: &[f64], median: f64) -> usize {
    let lm = median.log2();
    (0..peaks.len())
        .min_by(|&a, &b| (peaks[a] - lm).abs().total_cmp(&(peaks[b] - lm).abs()))
        .expect("peaks is not empty")
}

/// Which organelle-like reads are scanned: for each cluster, every `step`-th of its reads (in
/// input order) up to `sample`, `step` spreading the sample over the whole cluster.
fn sample_reads(
    candidates: &[(u64, usize)],
    clusters: usize,
    sample: usize,
) -> FxHashMap<u64, usize> {
    let mut size = vec![0usize; clusters];
    for &(_, c) in candidates {
        size[c] += 1;
    }
    let step: Vec<usize> = size.iter().map(|&n| (n / sample.max(1)).max(1)).collect();
    let (mut seen, mut taken) = (vec![0usize; clusters], vec![0usize; clusters]);
    let mut picked = FxHashMap::default();
    for &(ordinal, c) in candidates {
        if seen[c] % step[c] == 0 && taken[c] < sample {
            picked.insert(ordinal, c);
            taken[c] += 1;
        }
        seen[c] += 1;
    }
    picked
}

/// A read's class: the organelle of at least `min` of its gene stretches and more than the other's.
fn read_class(genes: &[u16], panel: &Panel, min: usize) -> Option<&'static str> {
    let mito = genes
        .iter()
        .filter(|&&g| panel.gene_organelle(g) == "mitochondrion")
        .count();
    let plastid = genes.len() - mito;
    if mito >= min && mito > plastid {
        Some("mitochondrion")
    } else if plastid >= min && plastid > mito {
        Some("plastid")
    } else {
        None
    }
}

/// The seed of a read: the span of its stretches of `organelle`'s genes plus `flank` bases each
/// side, None without such stretches.
fn crop(
    len: usize,
    stretches: &[(u16, usize, usize)],
    organelle: &str,
    panel: &Panel,
    flank: usize,
) -> Option<(usize, usize)> {
    let mine = stretches
        .iter()
        .filter(|&&(g, _, _)| panel.gene_organelle(g) == organelle);
    let (lo, hi) = mine.fold((usize::MAX, 0), |(lo, hi), &(_, a, b)| (lo.min(a), hi.max(b)));
    (lo < hi).then(|| (lo.saturating_sub(flank), (hi + flank).min(len)))
}

/// Fewest bases a depth cluster must stand for (its reads' bases over its depth) to seed: a deep
/// nuclear repeat with organelle-derived genes stands for a few hundred bases at most (Orcinus
/// orca, the PMAT2 benchmark's 1 Gb: a 101,070x cluster of 356 reads, 3.6 Mb, about 36 bp, with
/// atp8/nad3/nad4 hits, where the whole read set held 19 mitochondrial reads); the smallest
/// mitochondrial genomes are about 6 kb.
pub const MIN_CLUSTER_SPAN: f64 = 4000.0;

/// Each gene's stretches in excess of the reversed control's: the forward count less the control's
/// where the forward count is at least twice the control's, else none. The reversed reads (not
/// complemented) keep the cluster's composition and low-complexity runs but not real codon order,
/// so their stretches per gene are what chance alone assigns that gene in this cluster; a gene is
/// present only where the reads hold it beyond that level. Chance hits gather on the short,
/// divergent, hydrophobic animal genes the metazoa panel holds most proteins of (atp8 197, nad6
/// 102, nad4l 95, nad2 52, nad3 32 of 551), whose larger k-mer space low-complexity peptides hit
/// more often. Acipenser ruthenus (PMAT2 benchmark, first 1.9 Gb), the 332x cluster with
/// forward/control nad6 52/97, nad4l 17/4, nad2 6/0, atp8 3/3, nad4 3/12: 5 genes forwards, 2 in
/// excess (none of its recruited reads mitochondrial); Coccinella septempunctata 256x: 13 genes,
/// all in excess (12 at 165-402 against at most 2, atp8 15/3); Bugula neritina 91x: 12 genes, 11
/// in excess (nad6 54/39 is not).
fn excess_over_control(forward: &[u64], control: &[u64]) -> Vec<u64> {
    forward
        .iter()
        .zip(control)
        .map(|(&f, &c)| if f >= 2 * c { f - c } else { 0 })
        .collect()
}

/// Whether a class seeds from a depth cluster: its (genes, stretches per read) count
/// (`PanelParams::min_genes`, `min_density`), the cluster stands for at least `MIN_CLUSTER_SPAN`,
/// and the same reads reversed (not complemented) do not count for the class too. Reversal keeps a
/// read's composition and its low-complexity runs but destroys the order of real codons, so a
/// protein-coding organelle gives genes only forwards; a class that counts both ways is chance
/// hits of hydrophobic, low-complexity peptides on repeats, which match the short, divergent
/// animal genes the panel holds most of (atp8, nad6, nad4l). PMAT2 benchmark, 8bb78d1: the three
/// seeded clusters whose control counted were Acipenser ruthenus 140x and 70x (7 and 6 genes, mostly
/// nad6/nad4l/atp8; controls 5 and 9 genes; 646 recruited reads, none mitochondrial, the
/// mitochondrion about 5x) and Canis lupus 49x (its mitochondrion the 332x cluster); the controls
/// of all 76 other seeded clusters (69 samples) counted no class.
fn class_seeds(
    forward: (&ClassReport, f64),
    reversed: (&ClassReport, f64),
    span: f64,
    p: &PanelParams,
) -> bool {
    let counts = |(r, d): (&ClassReport, f64)| r.genes >= p.min_genes && d >= p.min_density;
    counts(forward) && !counts(reversed) && span >= MIN_CLUSTER_SPAN
}

fn panel_gene_sets(panel: &Panel) -> Vec<PanelGeneSet> {
    ORGANELLES
        .iter()
        .map(|&organelle| {
            let mut families: Vec<_> = (0..panel.gene_count())
                .map(|gene| gene as u16)
                .filter(|&gene| panel.gene_organelle(gene) == organelle)
                .map(|gene| panel.gene_name(gene).to_owned())
                .collect();
            families.sort_unstable();
            PanelGeneSet {
                organelle,
                family_count: families.len(),
                families,
            }
        })
        .collect()
}

struct GeneTally {
    stretches: Vec<u64>,
    records: Vec<u64>,
    last_record: Vec<Option<usize>>,
}

impl GeneTally {
    fn new(genes: usize) -> Self {
        Self {
            stretches: vec![0; genes],
            records: vec![0; genes],
            last_record: vec![None; genes],
        }
    }

    fn add(
        &mut self,
        panel: &Panel,
        record: usize,
        class: Option<&str>,
        genes: impl Iterator<Item = u16>,
        mut on_first_in_record: impl FnMut(usize),
    ) {
        for gene in genes {
            if Some(panel.gene_organelle(gene)) != class {
                continue;
            }
            let gene = usize::from(gene);
            self.stretches[gene] += 1;
            if self.last_record[gene] != Some(record) {
                self.records[gene] += 1;
                self.last_record[gene] = Some(record);
                on_first_in_record(gene);
            }
        }
    }
}

fn cluster_gene_counts(
    panel: &Panel,
    forward: &GeneTally,
    reversed: &GeneTally,
    both: &[u64],
) -> Vec<ClusterGeneCount> {
    let mut counts: Vec<_> = (0..panel.gene_count())
        .map(|i| {
            let gene = i as u16;
            ClusterGeneCount {
                organelle: panel.gene_organelle(gene),
                family: panel.gene_name(gene).to_owned(),
                forward_stretches: forward.stretches[i],
                reversed_stretches: reversed.stretches[i],
                forward_records: forward.records[i],
                reversed_records: reversed.records[i],
                both_orientation_records: both[i],
            }
        })
        .collect();
    counts.sort_unstable_by(|a, b| {
        a.organelle.cmp(b.organelle).then_with(|| a.family.cmp(&b.family))
    });
    counts
}

/// Protein seeds for the organelles in `want` (the others are not looked for in the whole-read
/// scan; depth clusters may still seed them).
pub fn run(
    inputs: &[PathBuf],
    panel: &Panel,
    p: &ProteinSeedParams,
    want: &[&'static str],
    out_dir: &Path,
) -> Result<ProteinSeedReport> {
    let t0 = Instant::now();
    std::fs::create_dir_all(out_dir)?;
    // passes 1-2: the depth sketch and the nuclear peak
    let mut gate = DepthGate::long_reads(inputs)?;
    match (p.nuclear_depth, p.genome_size) {
        (Some(d), _) => gate.set_nuclear_depth(d, "given"),
        (None, Some(g)) => {
            let d = ((gate.report.sketch_bases as f64 / g as f64).round() as u32).max(1);
            gate.set_nuclear_depth(d, "genome size")
        }
        (None, None) => {}
    }
    if let Some(f) = p.depth_fold {
        gate.set_depth_fold(f);
    }
    // pass 3: organelle-like reads and their median depths
    let (rx, reader) = read_batches(inputs);
    let mut liked: Vec<(u64, f64, u64)> = Vec::new();
    let mut reads = 0u64;
    for batch in rx {
        let first = reads;
        reads += batch.len() as u64;
        let found: Vec<(u64, f64, u64)> = batch
            .par_iter()
            .enumerate()
            .filter_map(|(i, r)| {
                let d = gate.profile(&r.seq)?;
                d.organelle_like
                    .then_some((first + i as u64, d.median, r.seq.len() as u64))
            })
            .collect();
        liked.extend(found);
    }
    reader.join().expect("reader thread panicked")?;
    let medians: Vec<f64> = liked.iter().map(|x| x.1).collect();
    let mut peaks = depth_peaks(&medians, p.min_cluster_reads);
    peaks.sort_by(|a, b| b.total_cmp(a)); // deepest first (usually the plastid)
    let candidates: Vec<(u64, usize)> = if peaks.is_empty() {
        Vec::new()
    } else {
        liked
            .iter()
            .map(|&(o, m, _)| (o, nearest(&peaks, m)))
            .collect()
    };
    let mut cluster_reads = vec![0u64; peaks.len()];
    let mut cluster_bases = vec![0u64; peaks.len()];
    for (&(_, c), &(_, _, len)) in candidates.iter().zip(&liked) {
        cluster_reads[c] += 1;
        cluster_bases[c] += len;
    }
    let picked = sample_reads(&candidates, peaks.len(), p.sample);
    drop(candidates);

    // pass 4: the sampled reads' gene stretches
    let mut scanned: Vec<Scanned> = Vec::new();
    if !picked.is_empty() {
        let (rx, reader) = read_batches(inputs);
        let mut ordinal = 0u64;
        for batch in rx {
            let first = ordinal;
            ordinal += batch.len() as u64;
            let mine: Vec<(usize, crate::recruit::Rec)> = batch
                .into_iter()
                .enumerate()
                .filter_map(|(i, r)| picked.get(&(first + i as u64)).map(|&c| (c, r)))
                .collect();
            let done: Vec<Scanned> = mine
                .into_par_iter()
                .map_init(
                    || vec![0u32; panel.gene_count()],
                    |hits, (cluster, r)| {
                        let stretches = panel::gene_stretches(&r.seq, panel, &p.panel, hits);
                        let reversed: Vec<u8> = r.seq.iter().rev().copied().collect();
                        let control = panel::gene_stretches(&reversed, panel, &p.panel, hits)
                            .into_iter()
                            .map(|x| x.0)
                            .collect();
                        let read = (!stretches.is_empty()).then(|| (r.id, r.seq));
                        Scanned {
                            cluster,
                            stretches,
                            control,
                            read,
                        }
                    },
                )
                .collect();
            scanned.extend(done);
        }
        reader.join().expect("reader thread panicked")?;
    }

    // calls and seeds
    let mut fasta: FxHashMap<&'static str, std::io::BufWriter<std::fs::File>> =
        FxHashMap::default();
    let mut seeds = FxHashMap::default();
    let (mut seed_reads, mut seed_bases) = (FxHashMap::default(), FxHashMap::default());
    let mut clusters = Vec::new();
    for (c, &peak) in peaks.iter().enumerate() {
        let mine: Vec<&Scanned> = scanned.iter().filter(|s| s.cluster == c).collect();
        let n = mine.len() as u64;
        // only the stretches of each read's own class count, so a homologous gene of the other
        // organelle (or a plastid insertion in a mitochondrial read) does not add to its class
        let (mut segs, mut ctrl) = (
            GeneTally::new(panel.gene_count()),
            GeneTally::new(panel.gene_count()),
        );
        let classes: Vec<Option<&'static str>> = mine
            .iter()
            .map(|s| {
                let genes: Vec<u16> = s.stretches.iter().map(|x| x.0).collect();
                read_class(&genes, panel, p.min_read_stretches)
            })
            .collect();
        let mut both = vec![0u64; panel.gene_count()];
        for (record, (s, class)) in mine.iter().zip(&classes).enumerate() {
            segs.add(panel, record, *class, s.stretches.iter().map(|x| x.0), |_| {});
            let cclass = read_class(&s.control, panel, p.min_read_stretches);
            ctrl.add(panel, record, cclass, s.control.iter().copied(), |gene| {
                both[gene] += u64::from(segs.last_record[gene] == Some(record));
            });
        }
        let gene_evidence = cluster_gene_counts(panel, &segs, &ctrl, &both);
        let (mitochondrion, plastid, md, pd, call) = panel::call_segments(
            panel,
            excess_over_control(&segs.stretches, &ctrl.stretches),
            n,
            &p.panel,
        );
        let (cm, cp, cmd, cpd, ccall) = panel::call_segments(panel, ctrl.stretches, n, &p.panel);
        // a class that counts here (enough genes, enough stretches per read) and not in the
        // reversed control seeds from its reads, also where both count (a cluster holding parts of
        // both genomes at one depth)
        let span = cluster_bases[c] as f64 / 2f64.powf(peak);
        let seeded: Vec<&'static str> = [((&mitochondrion, md), (&cm, cmd)), ((&plastid, pd), (&cp, cpd))]
            .iter()
            .zip(ORGANELLES)
            .filter(|((fwd, rev), _)| class_seeds(*fwd, *rev, span, &p.panel))
            .map(|(_, o)| o)
            .collect();
        let mut given = 0;
        for (s, class) in mine.iter().zip(&classes) {
            let Some(organelle) = *class else { continue };
            if !seeded.contains(&organelle) {
                continue;
            }
            {
                let Some((id, seq)) = &s.read else { continue };
                let Some((a, b)) = crop(seq.len(), &s.stretches, organelle, panel, p.flank)
                else {
                    continue;
                };
                if !fasta.contains_key(organelle) {
                    let path = out_dir.join(format!("protein_seeds.{organelle}.fasta"));
                    fasta.insert(
                        organelle,
                        std::io::BufWriter::new(std::fs::File::create(&path)?),
                    );
                    seeds.insert(organelle, path);
                }
                let w = fasta.get_mut(organelle).expect("just inserted");
                w.write_all(b">")?;
                w.write_all(id.split(|&b| b == b' ' || b == b'\t').next().unwrap_or(id))?;
                writeln!(w, " cluster={} span={}-{}", c + 1, a, b)?;
                w.write_all(&seq[a..b])?;
                w.write_all(b"\n")?;
                given += 1;
                *seed_reads.entry(organelle).or_insert(0) += 1;
                *seed_bases.entry(organelle).or_insert(0u64) += (b - a) as u64;
            }
        }
        clusters.push(ClusterCall {
            depth: 2f64.powf(peak),
            reads: cluster_reads[c],
            span_bp: span.round() as u64,
            scanned: n,
            mitochondrion,
            plastid,
            mitochondrion_density: md,
            plastid_density: pd,
            call,
            seeded,
            reversed_control_call: ccall,
            reversed_control_genes: (cm.genes, cp.genes),
            reversed_control_segments: cm.gene_segments.iter().chain(&cp.gene_segments).cloned().collect(),
            gene_evidence,
            seed_reads: given,
        });
    }
    for (_, mut w) in fasta {
        w.flush()?;
    }
    // organelles no cluster seeded: every read, without the depth gate
    let missing: Vec<&'static str> = want
        .iter()
        .copied()
        .filter(|o| !seeds.contains_key(o))
        .collect();
    let whole_read_scan = if missing.is_empty() {
        None
    } else {
        let w = whole_read_seeds(inputs, panel, p, &missing, out_dir, &gate)?;
        for (&o, &n) in &w.0.seed_reads {
            if n > 0 {
                seeds.insert(o, out_dir.join(format!("protein_seeds.{o}.fasta")));
                seed_reads.insert(o, n);
                seed_bases.insert(o, w.1.get(o).copied().unwrap_or(0));
            }
        }
        Some(w.0)
    };
    let (panel_proteins, panel_genes) = panel.sizes();
    let report = ProteinSeedReport {
        tool: "ovasm find protein-seeds",
        version: env!("CARGO_PKG_VERSION"),
        params: p.clone(),
        panel_proteins,
        panel_genes,
        panel_gene_sets: panel_gene_sets(panel),
        gene_evidence_record_unit: "input_fastx_record",
        gene_evidence_scope: "depth_cluster_records_after_read_classification",
        reverse_control: "reverse_sequence_without_complement",
        panel_max_evalue: panel.max_evalue(),
        gate: gate.report.clone(),
        reads,
        organelle_like_reads: liked.len() as u64,
        clusters,
        seed_reads,
        seed_bases,
        seeds,
        whole_read_scan,
        elapsed_seconds: t0.elapsed().as_secs_f64(),
        peak_rss_mb: peak_rss_mb(),
    };
    std::fs::write(
        out_dir.join("protein_seeds.json"),
        serde_json::to_string_pretty(&report)?,
    )?;
    Ok(report)
}

/// Whole reads as seeds for `organelles`, from a scan of every read: a read is a seed when its
/// stretches make it of the organelle's class (`read_class`) and cover at least
/// `WHOLE_READ_MIN_GENES` distinct genes. For reads where the organelle is not deeper than the
/// nucleus (the PMAT2 benchmark's read amounts: Orcinus orca 19 mitochondrial reads in 1 Gb, about
/// 6x; Sterechinus neumayeri 9 in 2.5 Gb), which the depth gate cannot set apart.
fn whole_read_seeds(
    inputs: &[PathBuf],
    panel: &Panel,
    p: &ProteinSeedParams,
    organelles: &[&'static str],
    out_dir: &Path,
    gate: &DepthGate,
) -> Result<(WholeReadScan, FxHashMap<&'static str, u64>)> {
    let mut above: FxHashMap<&'static str, usize> = FxHashMap::default();
    let t0 = Instant::now();
    let mut writers: FxHashMap<&'static str, std::io::BufWriter<std::fs::File>> =
        FxHashMap::default();
    let (mut n, mut b): (FxHashMap<&'static str, usize>, FxHashMap<&'static str, u64>) =
        Default::default();
    let (mut reads, mut bases) = (0u64, 0u64);
    let (rx, reader) = read_batches(inputs);
    for batch in rx {
        reads += batch.len() as u64;
        bases += batch.iter().map(|r| r.seq.len() as u64).sum::<u64>();
        let strict = PanelParams {
            min_hits: p.panel.min_hits.max(WHOLE_READ_MIN_HITS),
            ..p.panel.clone()
        };
        let hits: Vec<(usize, &'static str, usize, usize)> = batch
            .par_iter()
            .enumerate()
            .map_init(
                || vec![0u32; panel.gene_count()],
                |hits, (i, r)| {
                    let st = panel::gene_stretches(&r.seq, panel, &strict, hits);
                    let genes: Vec<u16> = st.iter().map(|x| x.0).collect();
                    let class = read_class(&genes, panel, p.min_read_stretches)?;
                    if !organelles.contains(&class) {
                        return None;
                    }
                    let mut distinct: Vec<u16> = genes
                        .iter()
                        .copied()
                        .filter(|&g| panel.gene_organelle(g) == class)
                        .collect();
                    distinct.sort_unstable();
                    distinct.dedup();
                    if distinct.len() < WHOLE_READ_MIN_GENES {
                        return None;
                    }
                    // the read's own genes and flanks, as a cluster seed (`crop`): the rest of a
                    // read can be nuclear (Gliophorus irrigatus: 18 of 21 seed reads 26-80%
                    // mitochondrial, the rest a GC-rich nuclear repeat that recruited 125 nuclear
                    // reads next to 23 mitochondrial ones)
                    let (a, b) = crop(r.seq.len(), &st, class, panel, p.flank)?;
                    Some((i, class, a, b))
                },
            )
            .flatten()
            .collect();
        for (i, class, lo, hi) in hits {
            let c = n.entry(class).or_insert(0);
            if *c >= WHOLE_READ_MAX_SEEDS {
                continue;
            }
            if !writers.contains_key(class) {
                let path = out_dir.join(format!("protein_seeds.{class}.fasta"));
                writers.insert(class, std::io::BufWriter::new(std::fs::File::create(&path)?));
            }
            let w = writers.get_mut(class).expect("just inserted");
            let r = &batch[i];
            w.write_all(b">")?;
            w.write_all(r.id.split(|&x| x == b' ' || x == b'\t').next().unwrap_or(&r.id))?;
            writeln!(w, " whole_read_scan span={lo}-{hi}")?;
            w.write_all(&r.seq[lo..hi])?;
            w.write_all(b"\n")?;
            *c += 1;
            *b.entry(class).or_insert(0) += (hi - lo) as u64;
            if seed_above_gate(gate, &r.seq, lo, hi) {
                *above.entry(class).or_insert(0) += 1;
            }
        }
    }
    reader.join().expect("reader thread panicked")?;
    for (_, mut w) in writers {
        w.flush()?;
    }
    for &o in organelles {
        n.entry(o).or_insert(0);
    }
    Ok((
        WholeReadScan {
            organelles: organelles.to_vec(),
            reads,
            bases,
            seed_reads: n,
            seed_reads_above_gate: above,
            gate_without_nucleus: gate.report.peak_window_bases < crate::seedfree::MIN_NUCLEAR_GENOME,
            seconds: t0.elapsed().as_secs_f64(),
        },
        b,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_go_to_the_nearest_depth_peak() {
        let peaks = [10f64.log2(), 100f64.log2()];
        assert_eq!(nearest(&peaks, 12.0), 0);
        assert_eq!(nearest(&peaks, 80.0), 1);
    }

    #[test]
    fn the_sample_is_spread_over_each_cluster() {
        // cluster 0: 10 reads, cluster 1: 3 reads, interleaved
        let cands: Vec<(u64, usize)> = (0..13u64)
            .map(|o| (o, if o % 4 == 3 { 1 } else { 0 }))
            .collect();
        let picked = sample_reads(&cands, 2, 4);
        let of = |c: usize| {
            let mut v: Vec<u64> = picked
                .iter()
                .filter(|(_, &x)| x == c)
                .map(|(&o, _)| o)
                .collect();
            v.sort();
            v
        };
        // cluster 0 (ordinals 0,1,2,4,5,6,8,9,10,12): step 10/4 = 2, every second, 4 at most
        assert_eq!(of(0), vec![0, 2, 5, 8]);
        // cluster 1 (3, 7, 11): fewer than the sample, all of them
        assert_eq!(of(1), vec![3, 7, 11]);
    }

    #[test]
    fn a_read_needs_two_stretches_of_its_class_and_a_majority() {
        let text = ">mito|cox1|Sp\nMAAAAAAAAAAAAAAAAAAAAA\n>plastid|rbcl|Sp\nMCCCCCCCCCCCCCCCCCCCCC\n";
        let panel = Panel::parse(text, 6, &[]).unwrap();
        let (m, pl) = if panel.gene_organelle(0) == "mitochondrion" {
            (0u16, 1u16)
        } else {
            (1, 0)
        };
        assert_eq!(read_class(&[m], &panel, 2), None);
        assert_eq!(read_class(&[m, m], &panel, 2), Some("mitochondrion"));
        assert_eq!(read_class(&[m, m, pl], &panel, 2), Some("mitochondrion"));
        assert_eq!(read_class(&[m, m, pl, pl], &panel, 2), None);
        assert_eq!(read_class(&[pl, pl, pl, m], &panel, 2), Some("plastid"));
    }

    #[test]
    fn a_seed_is_its_class_stretches_plus_flanks() {
        let text = ">mito|cox1|Sp\nMAAAAAAAAAAAAAAAAAAAAA\n>plastid|rbcl|Sp\nMCCCCCCCCCCCCCCCCCCCCC\n";
        let panel = Panel::parse(text, 6, &[]).unwrap();
        let (mito, pltd) = if panel.gene_organelle(0) == "mitochondrion" {
            (0u16, 1u16)
        } else {
            (1, 0)
        };
        let stretches = [(mito, 5000, 6000), (pltd, 100, 400), (mito, 9000, 9500)];
        assert_eq!(
            crop(20_000, &stretches, "mitochondrion", &panel, 2000),
            Some((3000, 11_500))
        );
        // flanks stop at the read's ends
        assert_eq!(
            crop(1000, &stretches, "plastid", &panel, 2000),
            Some((0, 1000))
        );
        assert_eq!(crop(1000, &[(pltd, 1, 2)], "mitochondrion", &panel, 10), None);
    }

    fn class(genes: usize, segments: u64) -> ClassReport {
        ClassReport {
            genes,
            segments,
            gene_segments: Vec::new(),
        }
    }

    #[test]
    fn gene_evidence_reports_effective_per_organelle_panel_sizes() {
        let text = ">mito|cox1|KeepA\nMAAAAAAAAAAAAAAAAAAAAA\n\
                    >mito|cox1|KeepB\nMAAAAAAAAAAAAAAAAAAAAA\n\
                    >mito|nad1|ExcludeA\nMAAAAAAAAAAAAAAAAAAAAA\n\
                    >plastid|cox1|KeepA\nMCCCCCCCCCCCCCCCCCCCCC\n";
        let panel = Panel::parse(text, 6, &[]).unwrap();
        let smaller = Panel::parse(text, 6, &["Exclude".to_owned()]).unwrap();
        let counts = panel_gene_sets(&panel);
        let filtered = panel_gene_sets(&smaller);
        assert_eq!(counts[0].organelle, "mitochondrion");
        assert_eq!(counts[0].family_count, 2);
        assert_eq!(counts[0].families, ["cox1", "nad1"]);
        assert_eq!(counts[1].organelle, "plastid");
        assert_eq!(counts[1].family_count, 1);
        assert_eq!(counts[1].families, ["cox1"]);
        assert_eq!(filtered[0].family_count, 1);
        assert_eq!(filtered[0].families, ["cox1"]);
        assert_eq!(filtered[1].families, ["cox1"]);
    }

    #[test]
    fn gene_evidence_counts_each_orientation_and_record_before_filtering() {
        let text = ">mito|cox1|Sp\nMAAAAAAAAAAAAAAAAAAAAA\n\
                    >mito|nad1|Sp\nMAAAAAAAAAAAAAAAAAAAAA\n\
                    >plastid|rbcL|Sp\nMCCCCCCCCCCCCCCCCCCCCC\n";
        let panel = Panel::parse(text, 6, &[]).unwrap();
        let mut forward = GeneTally::new(panel.gene_count());
        let mut reversed = GeneTally::new(panel.gene_count());
        let mut both = vec![0u64; panel.gene_count()];
        forward.add(&panel, 0, Some("mitochondrion"), [0, 0, 2].into_iter(), |_| {});
        reversed.add(&panel, 0, Some("mitochondrion"), [0, 0, 1, 1].into_iter(), |gene| {
            both[gene] += u64::from(forward.last_record[gene] == Some(0));
        });
        forward.add(&panel, 1, Some("mitochondrion"), [0, 1].into_iter(), |_| {});
        reversed.add(&panel, 1, Some("plastid"), [2].into_iter(), |gene| {
            both[gene] += u64::from(forward.last_record[gene] == Some(1));
        });
        let counts = cluster_gene_counts(&panel, &forward, &reversed, &both);
        assert_eq!(counts[0].family, "cox1");
        assert_eq!(counts[0].forward_stretches, 3);
        assert_eq!(counts[0].forward_records, 2);
        assert_eq!(counts[0].reversed_stretches, 2);
        assert_eq!(counts[0].reversed_records, 1);
        assert_eq!(counts[0].both_orientation_records, 1);
        assert_eq!(counts[1].forward_stretches, 1);
        assert_eq!(counts[1].reversed_stretches, 2);
        assert_eq!(counts[1].reversed_records, 1);
        assert_eq!(counts[1].both_orientation_records, 0);
        assert_eq!(counts[2].family, "rbcL");
        assert_eq!(counts[2].forward_records, 0);
        assert_eq!(counts[2].reversed_records, 1);
        assert_eq!(counts[2].both_orientation_records, 0);
        assert_eq!(
            excess_over_control(&forward.stretches, &reversed.stretches),
            [0, 0, 0]
        );
    }

    #[test]
    fn a_gene_counts_only_beyond_its_reversed_control() {
        // Acipenser ruthenus 332x: nad6, nad4l, nad2, atp8, nad4 forward against reversed
        let fwd = [52, 17, 6, 3, 3];
        let rev = [97, 4, 0, 3, 12];
        assert_eq!(excess_over_control(&fwd, &rev), vec![0, 13, 6, 0, 0]);
        // a mitochondrial gene with a chance background still counts its excess
        assert_eq!(excess_over_control(&[83, 54, 2], &[40, 39, 0]), vec![43, 0, 2]);
        // exactly twice the control: the excess equals the background
        assert_eq!(excess_over_control(&[10], &[5]), vec![5]);
    }

    #[test]
    fn a_class_that_counts_reversed_too_seeds_nothing() {
        let p = PanelParams::default();
        let none = class(0, 0);
        // a mitochondrial cluster: many genes forwards, a few chance ones reversed (Canis lupus 332x)
        assert!(class_seeds((&class(12, 4226), 4.17), (&class(3, 12), 0.01), 18_910.0, &p));
        // Acipenser ruthenus 70x: 7 genes at 0.47 forwards, 9 reversed at the same density
        assert!(!class_seeds((&class(7, 136), 0.47), (&class(9, 140), 0.49), 52_337.0, &p));
        // reversed genes without the density do not count, and do not stop it
        assert!(class_seeds((&class(7, 136), 0.47), (&class(9, 20), 0.07), 52_337.0, &p));
        // the forward class must count, and the cluster stand for enough sequence
        assert!(!class_seeds((&class(4, 400), 2.0), (&none, 0.0), 52_337.0, &p));
        assert!(!class_seeds((&class(12, 400), 0.2), (&none, 0.0), 52_337.0, &p));
        assert!(!class_seeds((&class(5, 762), 2.14), (&none, 0.0), 23.0, &p));
    }

    fn scan(seeds: &[(&'static str, usize, usize)]) -> WholeReadScan {
        WholeReadScan {
            organelles: seeds.iter().map(|x| x.0).collect(),
            reads: 1000,
            bases: 10_000_000,
            seed_reads: seeds.iter().map(|&(o, n, _)| (o, n)).collect(),
            seed_reads_above_gate: seeds.iter().map(|&(o, _, a)| (o, a)).collect(),
            gate_without_nucleus: false,
            seconds: 0.0,
        }
    }

    #[test]
    fn shallow_is_at_most_half_the_seeds_above_the_gate() {
        // Orcinus orca: the mitochondrion no deeper than the nucleus, no seed above the gate
        assert!(scan(&[("mitochondrion", 19, 0)]).shallow());
        // exactly half: shallow
        assert!(scan(&[("mitochondrion", 20, 10)]).shallow());
        // Acanthosoma haemorrhoidale: deep, its seeds above the gate, only its clusters called none
        assert!(!scan(&[("mitochondrion", 20, 11)]).shallow());
        // no seeds: nothing to judge
        assert!(!scan(&[("mitochondrion", 0, 0)]).shallow());
        // any shallow organelle makes the scan shallow
        assert!(scan(&[("mitochondrion", 20, 20), ("plastid", 4, 1)]).shallow());
        // a gate without a nucleus judges nothing: Gliophorus irrigatus, 13 of 21 above it
        let mut s = scan(&[("mitochondrion", 21, 13)]);
        assert!(!s.shallow());
        s.gate_without_nucleus = true;
        assert!(s.shallow());
        // still nothing to judge without seeds
        let mut s = scan(&[("mitochondrion", 0, 0)]);
        s.gate_without_nucleus = true;
        assert!(!s.shallow());
    }

    fn pseudo_random(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize]
            })
            .collect()
    }

    /// Reads of 3 kb every `step` bases along `genome`, as FASTA.
    fn tiled(genome: &[u8], step: usize, tag: &str) -> String {
        (0..genome.len() - 3000)
            .step_by(step)
            .enumerate()
            .map(|(i, st)| {
                format!(">{tag}{i}\n{}\n", std::str::from_utf8(&genome[st..st + 3000]).unwrap())
            })
            .collect()
    }

    #[test]
    fn a_partly_nuclear_seed_read_is_judged_on_its_seed() {
        // nucleus at 30x, organelle at 300x; a read of 3 kb of each, its seed the organelle half
        let nuclear = pseudo_random(200_000, 33);
        let organelle = pseudo_random(20_000, 7);
        let dir = std::env::temp_dir().join(format!("ovasm-seed-gate-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("reads.fa");
        std::fs::write(&fa, tiled(&nuclear, 100, "n") + &tiled(&organelle, 10, "m")).unwrap();
        let gate = DepthGate::long_reads(&[fa]).unwrap();
        let read: Vec<u8> = [&organelle[8000..11_000], &nuclear[50_000..53_000]].concat();
        assert!(
            !gate.profile(&read).is_some_and(|d| d.organelle_like),
            "the whole read is half nuclear"
        );
        assert!(seed_above_gate(&gate, &read, 0, 3000), "its seed is organelle-deep");
        assert!(!seed_above_gate(&gate, &read, 3000, 6000), "a nuclear span is not");
        std::fs::remove_dir_all(&dir).ok();
    }
}
