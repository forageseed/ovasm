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
//! The same reads reversed (not complemented) are scanned and called as a control. An organelle
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
pub struct ClusterCall {
    /// Centre of the depth peak (sketch depth).
    pub depth: f64,
    /// Organelle-like reads nearest this peak.
    pub reads: u64,
    /// The sequence those reads stand for: their bases over the depth (`MIN_CLUSTER_SPAN`).
    pub span_bp: u64,
    pub scanned: u64,
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
    /// Seed reads per organelle (whole reads, not cropped).
    pub seed_reads: FxHashMap<&'static str, usize>,
    pub seconds: f64,
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
    /// Whether some seeds came from the whole-read scan: the organelle is not deeper than the
    /// nucleus in these reads, so recruitment should not gate on depth either.
    pub fn shallow(&self) -> bool {
        self.whole_read_scan
            .as_ref()
            .is_some_and(|w| w.seed_reads.values().any(|&n| n > 0))
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
                "; whole-read scan for {} ({} reads): {} seed reads",
                w.organelles.join(" and "),
                w.reads,
                w.seed_reads.values().sum::<usize>()
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
        let (mut segs, mut ctrl) = (vec![0u64; panel.gene_count()], vec![0u64; panel.gene_count()]);
        let classes: Vec<Option<&'static str>> = mine
            .iter()
            .map(|s| {
                let genes: Vec<u16> = s.stretches.iter().map(|x| x.0).collect();
                read_class(&genes, panel, p.min_read_stretches)
            })
            .collect();
        for (s, class) in mine.iter().zip(&classes) {
            for &(g, _, _) in &s.stretches {
                if Some(panel.gene_organelle(g)) == *class {
                    segs[g as usize] += 1;
                }
            }
            let cclass = read_class(&s.control, panel, p.min_read_stretches);
            for &g in &s.control {
                if Some(panel.gene_organelle(g)) == cclass {
                    ctrl[g as usize] += 1;
                }
            }
        }
        let (mitochondrion, plastid, md, pd, call) =
            panel::call_segments(panel, segs, n, &p.panel);
        let (cm, cp, _, _, ccall) = panel::call_segments(panel, ctrl, n, &p.panel);
        // a class that counts here (enough genes, enough stretches per read) seeds from its reads,
        // also where both count (a cluster holding parts of both genomes at one depth)
        let counts = |r: &ClassReport, d: f64| r.genes >= p.panel.min_genes && d >= p.panel.min_density;
        let span = cluster_bases[c] as f64 / 2f64.powf(peak);
        let seeded: Vec<&'static str> = [(&mitochondrion, md), (&plastid, pd)]
            .iter()
            .zip(ORGANELLES)
            .filter(|((r, d), _)| counts(*r, *d) && span >= MIN_CLUSTER_SPAN)
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
        let w = whole_read_seeds(inputs, panel, p, &missing, out_dir)?;
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
) -> Result<(WholeReadScan, FxHashMap<&'static str, u64>)> {
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
        let hits: Vec<(usize, &'static str)> = batch
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
                    (distinct.len() >= WHOLE_READ_MIN_GENES).then_some((i, class))
                },
            )
            .flatten()
            .collect();
        for (i, class) in hits {
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
            w.write_all(b" whole_read_scan\n")?;
            w.write_all(&r.seq)?;
            w.write_all(b"\n")?;
            *c += 1;
            *b.entry(class).or_insert(0) += r.seq.len() as u64;
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
}
