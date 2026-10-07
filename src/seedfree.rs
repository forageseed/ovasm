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

fn build_sketch(inputs: &[PathBuf], p: &DiscoverParams, threshold: u64) -> Result<Sketch> {
    let (rx, reader) = read_batches(inputs);
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

/// Mode of the k-mer count histogram above the error spike (count 1).
fn nuclear_peak(counts: &FxHashMap<u64, u32>) -> u32 {
    let mut hist: FxHashMap<u32, u64> = FxHashMap::default();
    for &c in counts.values() {
        if c >= 2 {
            *hist.entry(c).or_insert(0) += 1;
        }
    }
    hist.into_iter()
        .max_by_key(|&(c, n)| (n, std::cmp::Reverse(c)))
        .map(|(c, _)| c)
        .unwrap_or(1)
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
    let kmers: Vec<u64> = sampled(seq, p.k, threshold).collect();
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

/// The organelle-like read test of `discover`, as a gate for reads found another way.
pub struct DepthGate {
    sketch: Sketch,
    params: DiscoverParams,
    threshold: u64,
    pub nuclear_peak_depth: u32,
    pub organelle_min_depth: f64,
}

impl DepthGate {
    /// Sketch the inputs with `discover`'s long-read defaults (scale 64, the whole input).
    pub fn long_reads(inputs: &[PathBuf]) -> Result<Self> {
        let params = DiscoverParams {
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
        };
        let threshold = u64::MAX / params.scale;
        let sketch = build_sketch(inputs, &params, threshold)?;
        let nuclear = nuclear_peak(&sketch.counts);
        Ok(Self {
            sketch,
            threshold,
            nuclear_peak_depth: nuclear,
            organelle_min_depth: nuclear as f64 * params.depth_fold,
            params,
        })
    }

    /// Deep and uniformly deep along the whole read: organelle, not a nuclear insertion.
    pub fn organelle_like(&self, seq: &[u8]) -> bool {
        features(
            seq,
            &self.sketch,
            &self.params,
            self.threshold,
            self.organelle_min_depth,
        )
        .is_some_and(|f| passes(&f, &self.params))
    }

    /// Median sampled k-mer depth along the read (sketch depth), if any k-mer was sampled.
    pub fn median_depth(&self, seq: &[u8]) -> Option<f64> {
        features(
            seq,
            &self.sketch,
            &self.params,
            self.threshold,
            self.organelle_min_depth,
        )
        .map(|f| f.median)
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
    let sketch = build_sketch(&inputs, &p, threshold)?;
    let nuclear = nuclear_peak(&sketch.counts);
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
        tool: "ovasm discover",
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
