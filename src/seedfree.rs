//! Seed-free organelle discovery.
//!
//! Organelle genomes are present at many times the nuclear copy number (Nipponbare HiFi:
//! mitochondria ~433x, plastids ~2,665x, nucleus ~80x), so their k-mers stand out in a depth
//! sketch without any reference. A read is organelle-like when its sampled k-mers are deep
//! AND uniformly deep along the read; high-copy nuclear repeats usually sit next to single-copy
//! flanks, which makes their depth profile uneven. Organelle-like reads are then split into
//! depth clusters (plastid and mitochondrion form separate peaks).

use std::path::PathBuf;
use std::time::Instant;

use anyhow::{bail, Result};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::Serialize;

use crate::kmer::CanonicalKmers;
use crate::recruit::{peak_rss_mb, read_batches, unit_hash, write_rec, Rec};

#[derive(Debug, Clone, Serialize)]
pub struct DiscoverParams {
    pub k: usize,
    /// Keep one k-mer in `scale` (FracMinHash) for the depth sketch.
    pub scale: u64,
    /// Bases read to build the sketch; 0 reads the whole input. A prefix stands for the whole
    /// only when reads come in random order, and runs are not always shuffled: in the v4
    /// holdout Vicia run, 7,414 of its 7,449 plastid reads came after the first 3 Gb of 3.29 Gb,
    /// so a 3 Gb sketch saw the plastid at ~5x and judged none of its reads organelle-like.
    pub sketch_bases: u64,
    /// Organelle threshold = nuclear peak depth x this factor.
    pub depth_fold: f64,
    /// Minimum sampled k-mers for a read to be judged.
    pub min_sampled: usize,
    /// Fraction of the read's sampled k-mers that must exceed the organelle threshold.
    pub min_high_frac: f64,
    /// Maximum (Q3 - Q1) / median of sampled k-mer depths: uniformity along the read.
    pub max_spread: f64,
    /// Minimum distinct / total sampled k-mers within a read. Tandem arrays (rDNA,
    /// centromeric satellites) are deep and uniform too, but repeat their k-mers in-read.
    pub min_distinct: f64,
    /// Minimum reads for a depth peak to become a cluster.
    pub min_cluster_reads: usize,
    pub target_depth: Option<f64>,
    pub salt: u64,
    /// Write every read's features to profiles.tsv (threshold tuning).
    pub dump_profiles: bool,
}

#[inline]
fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn sampled(seq: &[u8], k: usize, threshold: u64) -> impl Iterator<Item = u64> + '_ {
    CanonicalKmers::new(seq, k)
        .map(|(_, km)| km)
        .filter(move |&km| mix64(km) < threshold)
}

struct Sketch {
    counts: FxHashMap<u64, u32>,
    bases: u64,
}

/// `find discover`'s long-read defaults, as `DepthGate::long_reads` sketches with them (scale
/// 64, the whole input).
fn long_read_gate_params() -> DiscoverParams {
    DiscoverParams {
        k: 21,
        scale: 64,
        sketch_bases: 0,
        depth_fold: 3.0,
        min_sampled: 20,
        min_high_frac: 0.9,
        max_spread: 1.5,
        min_distinct: 0.8,
        min_cluster_reads: 20,
        target_depth: None,
        salt: 0,
        dump_profiles: false,
    }
}

/// The depth gate's first sketch pass, counted from reads another pass already reads
/// (recruitment's first pass): `sample` per read on any thread, `add` in any order. The counts do
/// not depend on the order, so the gate is the one `DepthGate::long_reads_from` builds.
pub struct GateSketch {
    sketch: Sketch,
    k: usize,
    threshold: u64,
}

impl Default for GateSketch {
    fn default() -> Self {
        let p = long_read_gate_params();
        Self {
            sketch: Sketch {
                counts: FxHashMap::default(),
                bases: 0,
            },
            k: p.k,
            threshold: u64::MAX / p.scale,
        }
    }
}

impl GateSketch {
    /// The sampled k-mers of one read.
    pub fn sample(k: usize, threshold: u64, seq: &[u8]) -> Vec<u64> {
        sampled(seq, k, threshold).collect()
    }

    /// `k` and `threshold` for `sample`.
    pub fn sampling(&self) -> (usize, u64) {
        (self.k, self.threshold)
    }

    /// One batch: the sampled k-mers of each read and the batch's bases.
    pub fn add(&mut self, reads: Vec<Vec<u64>>, bases: u64) {
        self.sketch.bases += bases;
        for v in reads {
            for km in v {
                *self.sketch.counts.entry(km).or_insert(0) += 1;
            }
        }
    }
}

/// Batches of reads and the thread producing them: the input files (`read_batches`), or the
/// recruitment cache that holds their sequences (`Replay::sequences`).
pub type Batches = (
    crossbeam_channel::Receiver<Vec<Rec>>,
    std::thread::JoinHandle<Result<()>>,
);

fn build_sketch(batches: Batches, p: &DiscoverParams, threshold: u64) -> Result<Sketch> {
    let (rx, reader) = batches;
    let mut sketch = Sketch {
        counts: FxHashMap::default(),
        bases: 0,
    };
    for batch in rx.iter() {
        sketch.bases += batch.iter().map(|r| r.seq.len() as u64).sum::<u64>();
        let kmers: Vec<Vec<u64>> = batch
            .par_iter()
            .map(|r| sampled(&r.seq, p.k, threshold).collect())
            .collect();
        for v in kmers {
            for km in v {
                *sketch.counts.entry(km).or_insert(0) += 1;
            }
        }
        if p.sketch_bases > 0 && sketch.bases >= p.sketch_bases {
            break;
        }
    }
    drop(rx); // stops the reader early; its send error is an expected shutdown
    let _ = reader.join();
    Ok(sketch)
}

fn histogram_peak(hist: &FxHashMap<u32, u64>) -> u32 {
    hist.iter()
        .max_by_key(|&(&c, &n)| (n, std::cmp::Reverse(c)))
        .map_or(1, |(&c, _)| c)
}

/// Raw mode retained for audits: recurrent sequencing errors can also occur at count 2.
#[cfg(test)]
fn nuclear_peak(counts: &FxHashMap<u64, u32>) -> u32 {
    let mut hist: FxHashMap<u32, u64> = FxHashMap::default();
    for &c in counts.values() {
        if c >= 2 {
            *hist.entry(c).or_insert(0) += 1;
        }
    }
    histogram_peak(&hist)
}

/// A recurrent k-mer contributes to the nuclear histogram only when its count agrees
/// with the middle half of a carrier read's depth profile. This excludes low-count
/// errors on deeper reads without excluding genuinely shallow genomic sequence.
fn context_supported_keys(keys: Vec<u64>, counts: &FxHashMap<u64, u32>) -> Vec<u64> {
    let mut depths: Vec<u32> = keys.iter().map(|key| counts[key]).collect();
    if depths.is_empty() {
        return Vec::new();
    }
    depths.sort_unstable();
    let q = |f: f64| depths[((depths.len() - 1) as f64 * f).round() as usize];
    let (lo, hi) = (q(0.25), q(0.75));
    keys.into_iter()
        .filter(|key| {
            let count = counts[key];
            count >= 2 && count >= lo && count <= hi
        })
        .collect()
}

fn read_context_histogram(
    batches: Batches,
    sketch: &Sketch,
    p: &DiscoverParams,
    threshold: u64,
) -> Result<FxHashMap<u32, u64>> {
    let (rx, reader) = batches;
    let mut supported = rustc_hash::FxHashSet::default();
    let mut bases = 0;
    for batch in rx.iter() {
        bases += batch.iter().map(|r| r.seq.len() as u64).sum::<u64>();
        let selected: Vec<Vec<u64>> = batch
            .par_iter()
            .map(|r| {
                context_supported_keys(sampled(&r.seq, p.k, threshold).collect(), &sketch.counts)
            })
            .collect();
        // Each distinct k-mer contributes once, even when many reads support it.
        supported.extend(selected.into_iter().flatten());
        if p.sketch_bases > 0 && bases >= p.sketch_bases {
            break;
        }
    }
    drop(rx);
    reader.join().expect("reader thread panicked")?;
    let mut histogram = FxHashMap::default();
    for key in supported {
        *histogram.entry(sketch.counts[&key]).or_insert(0) += 1;
    }
    Ok(histogram)
}

/// Per-read depth features from the sketch.
#[derive(Debug, Clone, Copy)]
struct Features {
    sampled: usize,
    distinct_frac: f64,
    high_frac: f64,
    spread: f64,
    median: f64,
}

fn features(
    seq: &[u8],
    sketch: &Sketch,
    p: &DiscoverParams,
    threshold: u64,
    organelle_min: f64,
) -> Option<Features> {
    features_of(&sampled(seq, p.k, threshold).collect::<Vec<u64>>(), sketch, organelle_min)
}

/// `features` of a read's sampled k-mers.
fn features_of(kmers: &[u64], sketch: &Sketch, organelle_min: f64) -> Option<Features> {
    if kmers.is_empty() {
        return None;
    }
    let distinct: rustc_hash::FxHashSet<u64> = kmers.iter().copied().collect();
    // k-mers seen at most once carry this read's own sequencing errors (even HiFi puts an error
    // in a few percent of 21-mers); they say nothing about genomic depth, so they are left out.
    // Nuclear flanks still count as low depth: their k-mers recur at the nuclear depth (>= 2).
    let mut depths: Vec<u32> = kmers
        .iter()
        .map(|km| sketch.counts.get(km).copied().unwrap_or(0))
        .filter(|&d| d >= 2)
        .collect();
    if depths.is_empty() {
        return None;
    }
    let high = depths
        .iter()
        .filter(|&&d| d as f64 >= organelle_min)
        .count();
    depths.sort_unstable();
    let q = |f: f64| depths[((depths.len() - 1) as f64 * f).round() as usize] as f64;
    let (q1, median, q3) = (q(0.25), q(0.5), q(0.75));
    Some(Features {
        sampled: kmers.len(),
        distinct_frac: distinct.len() as f64 / kmers.len() as f64,
        high_frac: high as f64 / depths.len() as f64,
        spread: if median > 0.0 {
            (q3 - q1) / median
        } else {
            f64::INFINITY
        },
        median,
    })
}

fn passes(f: &Features, p: &DiscoverParams) -> bool {
    f.sampled >= p.min_sampled
        && f.distinct_frac >= p.min_distinct
        && f.high_frac >= p.min_high_frac
        && f.median > 0.0
        && f.spread <= p.max_spread
}

/// Peaks of the log2(median depth) histogram; returns peak centres (log2 units).
fn depth_peaks(medians: &[f64], min_reads: usize) -> Vec<f64> {
    const BIN: f64 = 0.25;
    let mut hist: FxHashMap<i64, usize> = FxHashMap::default();
    for &m in medians {
        *hist.entry((m.log2() / BIN).floor() as i64).or_insert(0) += 1;
    }
    let raw = |b: i64| hist.get(&b).copied().unwrap_or(0);
    // Smoothed count first; the raw count breaks the plateau that a single-bin population
    // leaves in the smoothed profile (its neighbours then share the same smoothed sum).
    let key = |b: i64| ((b - 1..=b + 1).map(raw).sum::<usize>(), raw(b));
    let mut bins: Vec<i64> = hist.keys().copied().collect();
    bins.sort_unstable();
    let mut peaks = Vec::new();
    for &b in &bins {
        let here = key(b);
        if here.0 >= min_reads && here >= key(b - 1) && here > key(b + 1) {
            peaks.push((b as f64 + 0.5) * BIN);
        }
    }
    // Merge peaks closer than 1 log2 unit (2-fold): one genome, not two.
    peaks.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut merged: Vec<f64> = Vec::new();
    for pk in peaks {
        match merged.last_mut() {
            Some(last) if pk - *last < 1.0 => *last = (*last + pk) / 2.0,
            _ => merged.push(pk),
        }
    }
    merged
}

/// The organelle-like read test of `find discover`, as a gate for reads found another way.
pub struct DepthGate {
    sketch: Sketch,
    params: DiscoverParams,
    threshold: u64,
    pub nuclear_peak_depth: u32,
    pub organelle_min_depth: f64,
    pub report: GateReport,
}

/// How the gate's nuclear peak came about, for `recruit.json`.
#[derive(Debug, Clone, Serialize)]
pub struct GateReport {
    /// "inputs" or "recruitment cache": what the two sketch passes read.
    pub source: &'static str,
    pub sketch_bases: u64,
    pub sampled_kmers: usize,
    pub nuclear_peak_depth: u32,
    pub organelle_min_depth: f64,
    /// The organelle threshold over the nuclear peak (3, or `find recruit --gate-fold`).
    pub depth_fold: f64,
    /// "reads" (`nuclear_read_depth`), or "kmers" where no read-depth peak stands for a nuclear
    /// genome (the mode of the read-context k-mer histogram).
    pub nuclear_source: &'static str,
    /// Reads whose median depth could be judged (see `nuclear_read_depth`).
    pub judged_reads: u64,
    /// The genome size the peak stands for: the bases of the reads within half a log2 unit of it
    /// over its depth (for "kmers", the peak's distinct k-mers times the sampling scale).
    pub peak_window_bases: u64,
    /// The mode of the read-context k-mer histogram, the nuclear peak up to 0.2.1.
    pub kmer_peak_depth: u32,
    /// Judged reads per quarter-log2 bin of their median depth: (bin centre, reads).
    pub read_depth_histogram: Vec<(f64, u64)>,
    pub seconds: f64,
}

/// Width of the read-depth bins, in log2 units: 19% of the depth, so a bin holds the same share of
/// a Poisson-spread depth peak whatever its depth.
const DEPTH_BIN: f64 = 0.25;
/// Reads with a lower median depth are not judged: their recurrent k-mers are mostly at count 2,
/// recurrent errors or sequence present once or twice.
const MIN_JUDGED_DEPTH: u32 = 3;
/// A depth peak whose reads stand for less sequence than this is no nuclear genome's (the smallest
/// land-plant nuclear genomes are ~60 Mb, mitochondria at most ~12 Mb).
pub const MIN_NUCLEAR_GENOME: u64 = 20_000_000;

/// The nuclear depth and the genome size its reads stand for: of the `DEPTH_BIN` bins of the
/// judged reads' median depths (each summed with its neighbours, to smooth), the fullest whose
/// reads stand for `min_genome` (`MIN_NUCLEAR_GENOME`), as the median depth of the reads within
/// two bins of it; none in a library too shallow for the nucleus to be judged (below 3x), or one of
/// organelle reads alone. Also the binned histogram.
///
/// Reads, not k-mers. The k-mer estimate (the mode of `read_context_histogram`) counts distinct
/// k-mers, and a few percent of shallow reads (contaminants, sequence at 1-3x) bring more distinct
/// k-mers at count 2 than any single count of the nuclear peak holds: on the 2026-10-21 benchmark
/// libraries 176,546 against at most 140,853 (Carex echinata, at 57x), 211,677 against 187,296
/// (Medicago polymorpha) and 4.2 M against a peak spread over 15-30x (Lophophora williamsii). The
/// peak came out at 2x, the organelle threshold at 6x, and about 62% of all Carex reads passed the
/// gate. Most reads of a library are nuclear, so their depth mode is the nuclear depth however few
/// distinct k-mers it has; integer bins would favour the shallow end, where a Poisson peak is
/// narrow (Lophophora: 186,439 reads at median 2, 85-90 thousand at each of 19-24x). The read
/// estimate gives 61x for Carex, 55x for Medicago, 22x for Lophophora and 87x for Nipponbare,
/// where the k-mer estimate was right (84x).
fn nuclear_read_depth(
    medians: &FxHashMap<u32, (u64, u64)>,
    min_genome: u64,
) -> (Option<(u32, u64)>, Vec<(f64, u64)>) {
    let bin = |d: u32| ((d as f64).log2() / DEPTH_BIN).floor() as i64;
    let mut bins: FxHashMap<i64, u64> = FxHashMap::default();
    for (&d, &(reads, _)) in medians.iter().filter(|(&d, _)| d >= MIN_JUDGED_DEPTH) {
        *bins.entry(bin(d)).or_insert(0) += reads;
    }
    let mut histogram: Vec<(f64, u64)> = bins
        .iter()
        .map(|(&b, &n)| (2f64.powf((b as f64 + 0.5) * DEPTH_BIN), n))
        .collect();
    histogram.sort_by(|a, b| a.0.total_cmp(&b.0));
    let raw = |b: i64| bins.get(&b).copied().unwrap_or(0);
    let mut fullest: Vec<i64> = bins.keys().copied().collect();
    fullest.sort_by_key(|&b| {
        std::cmp::Reverse((raw(b - 1) + raw(b) + raw(b + 1), raw(b), std::cmp::Reverse(b)))
    });
    for peak in fullest {
        let mut window: Vec<(u32, u64, u64)> = medians
            .iter()
            .filter(|(&d, _)| d >= MIN_JUDGED_DEPTH && (bin(d) - peak).abs() <= 2)
            .map(|(&d, &(reads, bases))| (d, reads, bases))
            .collect();
        window.sort_unstable();
        let total: u64 = window.iter().map(|w| w.1).sum();
        let (mut seen, mut depth) = (0, window[0].0);
        for &(d, reads, _) in &window {
            seen += reads;
            depth = d;
            if 2 * seen >= total {
                break;
            }
        }
        let size = window.iter().map(|w| w.2).sum::<u64>() / u64::from(depth);
        if size >= min_genome {
            return (Some((depth, size)), histogram);
        }
    }
    (None, histogram)
}

/// One pass for the nuclear depth: the read-context k-mer histogram (the estimate up to 0.2.1,
/// kept for the report), and the reads and bases at each integer median depth among the reads
/// that can be judged (`min_sampled` sampled k-mers).
type NuclearPass = (FxHashMap<u32, u64>, FxHashMap<u32, (u64, u64)>);

fn nuclear_pass(batches: Batches, sketch: &Sketch, p: &DiscoverParams, threshold: u64) -> Result<NuclearPass> {
    let (rx, reader) = batches;
    let mut supported = rustc_hash::FxHashSet::default();
    let mut medians: FxHashMap<u32, (u64, u64)> = FxHashMap::default();
    for batch in rx.iter() {
        let per_read: Vec<(Vec<u64>, Option<u32>, u64)> = batch
            .par_iter()
            .map(|r| {
                let keys: Vec<u64> = sampled(&r.seq, p.k, threshold).collect();
                let median = features_of(&keys, sketch, 0.0)
                    .filter(|f| f.sampled >= p.min_sampled)
                    .map(|f| f.median as u32);
                (context_supported_keys(keys, &sketch.counts), median, r.seq.len() as u64)
            })
            .collect();
        for (keys, median, len) in per_read {
            supported.extend(keys);
            if let Some(m) = median {
                let e = medians.entry(m).or_insert((0, 0));
                e.0 += 1;
                e.1 += len;
            }
        }
    }
    reader.join().expect("reader thread panicked")?;
    let mut histogram = FxHashMap::default();
    for key in supported {
        *histogram.entry(sketch.counts[&key]).or_insert(0) += 1;
    }
    Ok((histogram, medians))
}

/// The gate's view of one read.
#[derive(Debug, Clone, Copy)]
pub struct ReadDepth {
    /// Deep and uniformly deep along the read (`organelle_like`).
    pub organelle_like: bool,
    /// Median sampled k-mer depth (sketch depth; about the read depth for HiFi reads).
    pub median: f64,
}

impl DepthGate {
    /// Sketch the inputs with `find discover`'s long-read defaults (scale 64, the whole input).
    pub fn long_reads(inputs: &[PathBuf]) -> Result<Self> {
        Self::long_reads_from(|| read_batches(inputs), "inputs")
    }

    /// `long_reads` over the batches `source` gives; it is called once per sketch pass.
    pub fn long_reads_from(source: impl Fn() -> Batches, label: &'static str) -> Result<Self> {
        Self::long_reads_with(source, label, MIN_NUCLEAR_GENOME)
    }

    /// `long_reads_from` with the sketch already counted (`GateSketch`, during recruitment's first
    /// pass over the inputs): only the nuclear pass reads `source`.
    pub fn long_reads_from_sketch(
        sketch: GateSketch,
        source: impl Fn() -> Batches,
        label: &'static str,
    ) -> Result<Self> {
        Self::long_reads_with_sketch(Some(sketch.sketch), source, label, MIN_NUCLEAR_GENOME)
    }

    fn long_reads_with(
        source: impl Fn() -> Batches,
        label: &'static str,
        min_genome: u64,
    ) -> Result<Self> {
        Self::long_reads_with_sketch(None, source, label, min_genome)
    }

    fn long_reads_with_sketch(
        sketch: Option<Sketch>,
        source: impl Fn() -> Batches,
        label: &'static str,
        min_genome: u64,
    ) -> Result<Self> {
        let start = Instant::now();
        let params = long_read_gate_params();
        let threshold = u64::MAX / params.scale;
        let sketch = match sketch {
            Some(sketch) => sketch,
            None => build_sketch(source(), &params, threshold)?,
        };
        let (kmers, medians) = nuclear_pass(source(), &sketch, &params, threshold)?;
        let kmer_peak = histogram_peak(&kmers);
        let (by_reads, read_depth_histogram) = nuclear_read_depth(&medians, min_genome);
        // No judged depth peak stands for a nuclear genome: the nucleus is below 3x (the k-mer
        // estimate is then right), or the reads are organelle reads only.
        let (nuclear, window_bases, nuclear_source) = match by_reads {
            Some((depth, size)) => (depth, size, "reads"),
            None => {
                let window: u64 = kmers
                    .iter()
                    .filter(|&(&c, _)| 2 * c >= kmer_peak && c <= 2 * kmer_peak)
                    .map(|(_, &n)| n)
                    .sum();
                (kmer_peak, window * params.scale, "kmers")
            }
        };
        let organelle_min = nuclear as f64 * params.depth_fold;
        let report = GateReport {
            source: label,
            sketch_bases: sketch.bases,
            sampled_kmers: sketch.counts.len(),
            nuclear_peak_depth: nuclear,
            organelle_min_depth: organelle_min,
            depth_fold: params.depth_fold,
            nuclear_source,
            judged_reads: read_depth_histogram.iter().map(|h| h.1).sum(),
            peak_window_bases: window_bases,
            kmer_peak_depth: kmer_peak,
            read_depth_histogram,
            seconds: start.elapsed().as_secs_f64(),
        };
        Ok(Self {
            sketch,
            threshold,
            nuclear_peak_depth: nuclear,
            organelle_min_depth: organelle_min,
            params,
            report,
        })
    }

    /// Deep and uniformly deep along the whole read: organelle, not a nuclear insertion.
    pub fn organelle_like(&self, seq: &[u8]) -> bool {
        self.profile(seq).is_some_and(|d| d.organelle_like)
    }

    /// Median sampled k-mer depth along the read (sketch depth), if any k-mer was sampled.
    pub fn median_depth(&self, seq: &[u8]) -> Option<f64> {
        self.profile(seq).map(|d| d.median)
    }

    /// The organelle threshold at `fold` times the nuclear peak instead of `discover`'s 3.
    pub fn set_depth_fold(&mut self, fold: f64) {
        self.params.depth_fold = fold;
        self.organelle_min_depth = self.nuclear_peak_depth as f64 * fold;
        self.report.organelle_min_depth = self.organelle_min_depth;
        self.report.depth_fold = fold;
    }

    /// Both of the above from one look at the read's sampled k-mers.
    pub fn profile(&self, seq: &[u8]) -> Option<ReadDepth> {
        features(
            seq,
            &self.sketch,
            &self.params,
            self.threshold,
            self.organelle_min_depth,
        )
        .map(|f| ReadDepth {
            organelle_like: passes(&f, &self.params),
            median: f.median,
        })
    }

    /// `profile` of the read's sampled k-mers that `keep` selects, its deep ones counted from
    /// `deep` up (`organelle_min_depth` in `profile`).
    pub fn profile_of(
        &self,
        seq: &[u8],
        deep: f64,
        keep: impl Fn(u64) -> bool,
    ) -> Option<ReadDepth> {
        let kmers: Vec<u64> = sampled(seq, self.params.k, self.threshold)
            .filter(|&km| keep(km))
            .collect();
        features_of(&kmers, &self.sketch, deep).map(|f| ReadDepth {
            organelle_like: passes(&f, &self.params),
            median: f.median,
        })
    }

    /// Word length of the sketch.
    pub fn k(&self) -> usize {
        self.params.k
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct ClusterOutput {
    pub name: String,
    /// Median sampled k-mer depth in the sketch prefix.
    pub sketch_depth: f64,
    /// Scaled to the full input (sketch depth x total bases / sketch bases).
    pub estimated_depth: f64,
    pub reads: u64,
    pub bases: u64,
    pub sampling_fraction: f64,
    pub kept_reads: u64,
    pub output: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct DiscoverReport {
    pub tool: &'static str,
    pub version: &'static str,
    pub params: DiscoverParams,
    pub inputs: Vec<String>,
    pub sketch_bases: u64,
    pub sketch_kmers: usize,
    pub nuclear_peak_depth: u32,
    pub organelle_min_depth: f64,
    pub reads_scanned: u64,
    pub bases_scanned: u64,
    pub organelle_like_reads: u64,
    pub clusters: Vec<ClusterOutput>,
    pub elapsed_seconds: f64,
    pub peak_rss_mb: Option<f64>,
}

pub fn run(inputs: Vec<PathBuf>, p: DiscoverParams, out_dir: PathBuf) -> Result<DiscoverReport> {
    let start = Instant::now();
    std::fs::create_dir_all(&out_dir)?;
    if p.scale == 0 {
        bail!("--scale must be positive");
    }
    let threshold = u64::MAX / p.scale;
    let sketch = build_sketch(read_batches(&inputs), &p, threshold)?;
    let nuclear = histogram_peak(&read_context_histogram(
        read_batches(&inputs),
        &sketch,
        &p,
        threshold,
    )?);
    let organelle_min = nuclear as f64 * p.depth_fold;
    eprintln!(
        "[ovasm] sketch: {:.2} Gbp, {} sampled k-mers, nuclear peak {}x, organelle threshold {:.0}x (sketch depth)",
        sketch.bases as f64 / 1e9, sketch.counts.len(), nuclear, organelle_min
    );

    let (rx, reader) = read_batches(&inputs);
    let (mut reads, mut bases) = (0u64, 0u64);
    // Organelle-like reads go to disk as they are found, their median depths stay in memory:
    // a plastid at thousands of x is gigabases of reads, more than a laptop holds.
    let candidates_path = out_dir.join(".candidates.fastq");
    let mut candidates_w =
        std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(&candidates_path)?);
    let mut medians: Vec<f64> = Vec::new();
    let mut dump = if p.dump_profiles {
        let mut w = std::io::BufWriter::new(std::fs::File::create(out_dir.join("profiles.tsv"))?);
        use std::io::Write;
        writeln!(
            w,
            "read_id	sampled	distinct_frac	high_frac	spread	median	pass"
        )?;
        Some(w)
    } else {
        None
    };
    for batch in rx {
        reads += batch.len() as u64;
        bases += batch.iter().map(|r| r.seq.len() as u64).sum::<u64>();
        let scored: Vec<(Rec, Option<Features>)> = batch
            .into_par_iter()
            .map(|r| {
                let f = features(&r.seq, &sketch, &p, threshold, organelle_min);
                (r, f)
            })
            .collect();
        for (r, f) in scored {
            let Some(f) = f else { continue };
            let ok = passes(&f, &p);
            if let Some(w) = dump.as_mut() {
                use std::io::Write;
                writeln!(
                    w,
                    "{}	{}	{:.3}	{:.3}	{:.3}	{}	{}",
                    String::from_utf8_lossy(&r.id),
                    f.sampled,
                    f.distinct_frac,
                    f.high_frac,
                    f.spread,
                    f.median,
                    ok as u8
                )?;
            }
            if ok {
                write_rec(&mut candidates_w, &r)?;
                medians.push(f.median);
            }
        }
    }
    if let Some(mut w) = dump {
        use std::io::Write;
        w.flush()?;
    }
    reader.join().expect("reader thread panicked")?;
    {
        use std::io::Write;
        candidates_w.flush()?;
    }
    drop(candidates_w);

    let mut peaks = depth_peaks(&medians, p.min_cluster_reads);
    peaks.sort_by(|a, b| b.partial_cmp(a).unwrap()); // highest depth first (usually plastid)
    let scale_to_full = if sketch.bases > 0 {
        bases as f64 / sketch.bases as f64
    } else {
        1.0
    };
    let fractions: Vec<f64> = peaks
        .iter()
        .map(|peak| match p.target_depth {
            Some(d) if 2f64.powf(*peak) * scale_to_full > d => {
                d / (2f64.powf(*peak) * scale_to_full)
            }
            _ => 1.0,
        })
        .collect();
    let paths: Vec<PathBuf> = (0..peaks.len())
        .map(|i| out_dir.join(format!("cluster{}.fastq", i + 1)))
        .collect();
    let mut writers = paths
        .iter()
        .map(|path| {
            Ok(std::io::BufWriter::with_capacity(
                1 << 20,
                std::fs::File::create(path)?,
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    // (reads, bases, kept reads) per cluster
    let mut tally = vec![(0u64, 0u64, 0u64); peaks.len()];
    if !peaks.is_empty() {
        let mut parser = needletail::parse_fastx_file(&candidates_path)?;
        let mut n = 0usize;
        while let Some(record) = parser.next() {
            let record = record?;
            let lm = medians[n].log2();
            n += 1;
            let i = (0..peaks.len())
                .min_by(|&a, &b| (peaks[a] - lm).abs().total_cmp(&(peaks[b] - lm).abs()))
                .expect("peaks is not empty");
            tally[i].0 += 1;
            tally[i].1 += record.num_bases() as u64;
            if fractions[i] >= 1.0 || unit_hash(record.id(), p.salt) < fractions[i] {
                tally[i].2 += 1;
                write_rec(
                    &mut writers[i],
                    &Rec {
                        id: record.id().to_vec(),
                        seq: record.seq().into_owned(),
                        qual: record.qual().map(<[u8]>::to_vec),
                    },
                )?;
            }
        }
    }
    for w in &mut writers {
        use std::io::Write;
        w.flush()?;
    }
    std::fs::remove_file(&candidates_path)?;
    let mut clusters = Vec::new();
    for (i, peak) in peaks.iter().enumerate() {
        let sketch_depth = 2f64.powf(*peak);
        clusters.push(ClusterOutput {
            name: format!("cluster{}", i + 1),
            sketch_depth,
            estimated_depth: sketch_depth * scale_to_full,
            reads: tally[i].0,
            bases: tally[i].1,
            sampling_fraction: fractions[i],
            kept_reads: tally[i].2,
            output: paths[i].display().to_string(),
        });
    }
    let report = DiscoverReport {
        tool: "ovasm find discover",
        version: env!("CARGO_PKG_VERSION"),
        inputs: inputs.iter().map(|p| p.display().to_string()).collect(),
        sketch_bases: sketch.bases,
        sketch_kmers: sketch.counts.len(),
        nuclear_peak_depth: nuclear,
        organelle_min_depth: organelle_min,
        reads_scanned: reads,
        bases_scanned: bases,
        organelle_like_reads: medians.len() as u64,
        clusters,
        elapsed_seconds: start.elapsed().as_secs_f64(),
        peak_rss_mb: peak_rss_mb(),
        params: p,
    };
    std::fs::write(
        out_dir.join("discover.json"),
        serde_json::to_string_pretty(&report)?,
    )?;
    Ok(report)
}

#[cfg(test)]
#[path = "seedfree_diagnostics.rs"]
mod diagnostic_tests;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nuclear_peak_skips_the_error_spike() {
        let mut counts = FxHashMap::default();
        for i in 0..1000u64 {
            counts.insert(i, 1); // sequencing errors
        }
        for i in 1000..1300u64 {
            counts.insert(i, 7); // nuclear
        }
        for i in 1300..1350u64 {
            counts.insert(i, 40); // organelle
        }
        assert_eq!(nuclear_peak(&counts), 7);
    }

    #[test]
    fn read_context_rejects_recurrent_errors_without_a_fixed_depth_floor() {
        let mut counts = FxHashMap::default();
        let mut supported = rustc_hash::FxHashSet::default();
        // Three recurrent error keys per read dominate the raw distinct-key mode.
        // Their carrier reads otherwise agree at depth 17.
        for i in 0..100u64 {
            counts.insert(i, 17);
        }
        for i in 100..400u64 {
            counts.insert(i, 2);
        }
        for i in 400..450u64 {
            counts.insert(i, 100); // organelle: deep, but fewer distinct keys
        }
        assert_eq!(nuclear_peak(&counts), 2);
        for i in 0..100u64 {
            let mut keys: Vec<_> = (0..20).map(|j| (i + j) % 100).collect();
            keys.extend(100 + i * 3..103 + i * 3);
            supported.extend(context_supported_keys(keys, &counts));
        }
        supported.extend(context_supported_keys((400..450).collect(), &counts));
        let mut hist = FxHashMap::default();
        for key in supported {
            *hist.entry(counts[&key]).or_insert(0) += 1;
        }
        assert_eq!(histogram_peak(&hist), 17);
        assert!(!hist.contains_key(&2));

        // Truly shallow genomic reads remain usable; count 2 is not an error label.
        let shallow = context_supported_keys((100..400).collect(), &counts);
        assert_eq!(shallow.len(), 300);
        assert!(context_supported_keys(Vec::new(), &counts).is_empty());
    }

    #[test]
    fn the_nuclear_depth_is_the_mode_of_read_depths_in_log_bins() {
        // As in Lophophora williamsii: shallow reads at 2-3x, more of them than at any single
        // depth of a broad nuclear peak at 15-30x, and the organelles far deeper.
        let mut medians: FxHashMap<u32, (u64, u64)> = FxHashMap::default();
        medians.insert(2, (186_000, 186_000 * 12_000));
        medians.insert(3, (91_000, 91_000 * 12_000));
        for d in 15..=30u32 {
            let n = 85_000 - 3_000 * u64::from(d.abs_diff(22));
            medians.insert(d, (n, n * 12_000));
        }
        medians.insert(114, (6_700, 6_700 * 12_000));
        medians.insert(442, (10_000, 10_000 * 12_000));
        let (found, histogram) = nuclear_read_depth(&medians, MIN_NUCLEAR_GENOME);
        let (depth, size) = found.unwrap();
        assert_eq!(depth, 22);
        // the window's reads over that depth: a genome of hundreds of megabases
        assert!(size > 500_000_000, "{size}");
        assert!(histogram.iter().all(|&(centre, _)| centre >= 3.0));
        // nothing judged, or organelle reads alone (a 20 kb plastid at 300x): no nuclear peak
        let none = |m: &[(u32, (u64, u64))]| {
            nuclear_read_depth(&FxHashMap::from_iter(m.iter().copied()), MIN_NUCLEAR_GENOME).0
        };
        assert_eq!(none(&[(2, (10, 10))]), None);
        assert_eq!(none(&[(300, (2000, 2000 * 3000)), (100, (400, 400 * 3000))]), None);
    }

    #[test]
    fn shallow_reads_no_longer_make_the_nuclear_peak() {
        // A 300 kb nuclear genome at 20x beside 3 Mb of a contaminant at 2x: the contaminant holds
        // ten times the distinct k-mers, all at count 2, and so the read-context k-mer peak; the
        // reads are mostly nuclear.
        let genome = |n: usize, seed: u64| -> Vec<u8> {
            let mut x = seed;
            (0..n)
                .map(|_| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    b"ACGT"[(x % 4) as usize]
                })
                .collect()
        };
        let (nuclear, other) = (genome(300_000, 5), genome(3_000_000, 6));
        let mut text = String::new();
        for (tag, g, step) in [("n", &nuclear, 150), ("c", &other, 1500)] {
            for (i, st) in (0..g.len() - 3000).step_by(step).enumerate() {
                text.push_str(&format!(">{tag}{i}\n{}\n", std::str::from_utf8(&g[st..st + 3000]).unwrap()));
            }
        }
        let dir = std::env::temp_dir().join(format!("ovasm-gate-depth-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("reads.fa");
        std::fs::write(&fa, text).unwrap();
        // a 300 kb genome stands in for a nuclear one here
        let gate = DepthGate::long_reads_with(|| read_batches(&[fa.clone()]), "inputs", 200_000)
            .unwrap();
        assert_eq!(gate.report.kmer_peak_depth, 2, "the k-mer estimate");
        assert_eq!(gate.report.nuclear_source, "reads");
        assert!(
            (17..=21).contains(&gate.nuclear_peak_depth),
            "{}",
            gate.nuclear_peak_depth
        );
        assert!(gate.organelle_min_depth >= 51.0);
        // at the real size floor the 300 kb genome is no nucleus, and the k-mer estimate stands
        let small = DepthGate::long_reads(&[fa]).unwrap();
        assert_eq!(
            (small.nuclear_peak_depth, small.report.nuclear_source),
            (2, "kmers")
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn histogram_peak_breaks_ties_toward_the_lower_depth() {
        let hist = FxHashMap::from_iter([(3, 100), (7, 100), (30, 20)]);
        assert_eq!(histogram_peak(&hist), 3);
        assert_eq!(histogram_peak(&FxHashMap::default()), 1);
    }

    #[test]
    fn two_depth_populations_give_two_peaks() {
        let mut medians = vec![40.0; 300];
        medians.extend(vec![240.0; 800]);
        medians.extend([35.0, 45.0, 250.0, 230.0]);
        let mut peaks = depth_peaks(&medians, 20);
        peaks.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(peaks.len(), 2, "{peaks:?}");
        assert!((2f64.powf(peaks[0]) - 40.0).abs() < 10.0);
        assert!((2f64.powf(peaks[1]) - 240.0).abs() < 60.0);
    }
}
