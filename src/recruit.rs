//! Single-pass organelle read recruitment with iterative baiting, NUMT flagging and
//! deterministic depth-targeted downsampling.

use std::fs::File;
use std::hash::Hasher;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result};
use rayon::prelude::*;
use rustc_hash::{FxHashSet, FxHasher};
use serde::Serialize;

use crate::bait::BaitSet;
use crate::kmer::CanonicalKmers;
use crate::seedfree::DepthGate;

#[derive(Debug, Clone, Serialize)]
pub struct ClassifyParams {
    /// Minimum fraction of the read spanned by organelle-like blocks to call it organellar.
    /// Blocks bridge gaps up to `block_gap`, so sequencing errors (ONT) and seed/sample SNVs
    /// (short reads) do not break coverage, while NUMT flanks (kilobases of nuclear DNA) do.
    pub min_cover: f64,
    /// Minimum fraction of the read's k-mers that hit the bait set.
    pub min_hit_frac: f64,
    /// A read below `min_cover` but holding a contiguous organelle-like block at least this
    /// long is flagged as a NUMT/NUPT candidate (organelle segment with non-organelle flanks).
    pub numt_min_block: usize,
    /// Hits closer than this many bases still extend the same block (tolerates SNPs/indels).
    pub block_gap: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Call {
    /// Bitmask of targets: several bits when the read fits them equally well (e.g. a short
    /// read inside an MTPT identical in both genomes); it then belongs to each assembly.
    Organelle(u8),
    NumtCandidate(usize),
    Other,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct ReadScore {
    pub cover_frac: f64,
    pub hit_frac: f64,
    pub longest_block: usize,
}

/// Bait hits of one read, per target.
struct TargetHits {
    /// Bait k-mers hit.
    hits: [usize; 8],
    /// Summed length of gap-bridged blocks: the coverage used for the organelle call.
    block_cover: [usize; 8],
    longest: [usize; 8],
    /// k-mers in the read.
    total: usize,
}

fn target_hits(seq: &[u8], bait: &BaitSet, block_gap: usize) -> TargetHits {
    let n = bait.targets.len();
    let k = bait.k;
    let mut h = TargetHits {
        hits: [0; 8],
        block_cover: [0; 8],
        longest: [0; 8],
        total: 0,
    };
    let mut last_end = [0usize; 8];
    let mut block_start = [0usize; 8];
    for (pos, kmer) in CanonicalKmers::new(seq, k) {
        h.total += 1;
        let mask = bait.get(kmer);
        if mask == 0 {
            continue;
        }
        let end = pos + k;
        for t in 0..n {
            if mask & (1 << t) == 0 {
                continue;
            }
            if h.hits[t] == 0 {
                block_start[t] = pos;
            } else if pos > last_end[t] + block_gap {
                let len = last_end[t] - block_start[t];
                h.block_cover[t] += len;
                h.longest[t] = h.longest[t].max(len);
                block_start[t] = pos;
            }
            last_end[t] = last_end[t].max(end);
            h.hits[t] += 1;
        }
    }
    for t in 0..n {
        if h.hits[t] > 0 {
            let len = last_end[t] - block_start[t];
            h.block_cover[t] += len;
            h.longest[t] = h.longest[t].max(len);
        }
    }
    h
}

/// Score a read against every target and decide what it is.
pub fn classify(seq: &[u8], bait: &BaitSet, p: &ClassifyParams) -> (Call, ReadScore) {
    let n = bait.targets.len();
    let TargetHits {
        hits,
        block_cover,
        longest,
        total,
    } = target_hits(seq, bait, p.block_gap);

    // Gap-bridged coverage can tie (short reads inside MTPTs cover both genomes fully once
    // SNVs are bridged); exact k-mer hits then separate the true source from the donor copy.
    let key = |t: usize| (block_cover[t], hits[t]);
    let best = (0..n).max_by_key(|&t| key(t)).unwrap_or(0);
    let tied: u8 = (0..n)
        .filter(|&t| key(t) == key(best))
        .fold(0, |m, t| m | (1 << t));
    if n == 0 || total == 0 || hits[best] == 0 {
        return (Call::Other, ReadScore::default());
    }
    let score = ReadScore {
        cover_frac: block_cover[best] as f64 / seq.len() as f64,
        hit_frac: hits[best] as f64 / total as f64,
        longest_block: longest[best],
    };
    let call = if score.cover_frac >= p.min_cover && score.hit_frac >= p.min_hit_frac {
        Call::Organelle(tied)
    } else if score.longest_block >= p.numt_min_block {
        Call::NumtCandidate(best)
    } else {
        Call::Other
    };
    (call, score)
}

#[derive(Clone)]
pub struct Rec {
    pub id: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct TargetCounts {
    pub reads: u64,
    pub bases: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct IterationStats {
    pub iteration: usize,
    pub seconds: f64,
    pub reads_scanned: u64,
    pub bases_scanned: u64,
    pub bait_kmers: usize,
    pub recruited: Vec<TargetCounts>,
    pub numt_candidates: TargetCounts,
    /// Reads written to more than one target (tied fit, e.g. identical MTPT/NUPT).
    pub multi_target_reads: u64,
    /// k-mers added to each target's bait set after this pass (iterative baiting).
    pub bait_kmers_added: Vec<usize>,
    /// Reads each target's bait was extended from after this pass: only reads not used in an
    /// earlier pass, subsampled to `extend_depth`.
    pub extended_from_reads: Vec<u64>,
    /// Targets no longer extended: their recruited reads grew by less than `saturation`.
    pub saturated: Vec<bool>,
    /// Partially matching, organelle-deep reads each target's bait was bootstrapped from after
    /// this pass (`--bootstrap`, see `bootstrap_candidates`).
    pub bootstrapped_from_reads: Vec<u64>,
    /// Reads left out of each target's extension by the depth gate (`extend_gate`).
    pub extension_gated_out: Vec<u64>,
    /// Organelle-like reads left out of each target's extension for being shallower than
    /// `EXTEND_CORE_FLOOR` of its core depth (`extend_gate`).
    pub extension_shallow_out: Vec<u64>,
    /// And for being deeper than `HIGH_COPY_FOLD` times its core depth: a high-copy repeat.
    pub extension_deep_out: Vec<u64>,
    /// Sketch depth of each target's reads of this pass (`target_depth` of a sample of
    /// `DEPTH_SAMPLE`), when the depth sketch exists.
    pub read_depth: Vec<Option<f64>>,
    /// The same for the reads new to the target in this pass.
    pub fresh_read_depth: Vec<Option<f64>>,
    /// Recruited bases over the read depth: the genome size the target's reads stand for.
    pub implied_size: Vec<Option<u64>>,
}

/// Counts of one pass; its recruited reads went to the pass files (see `scan`).
struct ScanResult {
    recruited: Vec<TargetCounts>,
    /// Recruited reads not yet used to extend the target (`used` in `run`).
    fresh: Vec<TargetCounts>,
    /// Depth samples of each target's reads: all of them, and the fresh ones.
    all_sample: Vec<Sample>,
    fresh_sample: Vec<Sample>,
    numt: TargetCounts,
    reads: u64,
    bases: u64,
    multi_target: u64,
    reclassified: u64,
}

/// Reads per depth estimate. A median of 2,000 read depths is within a few percent, and the
/// sample costs one pass's memory for about 2,000 reads per target.
const DEPTH_SAMPLE: usize = 2000;

/// The `DEPTH_SAMPLE` reads with the smallest name hashes among those offered: the same reads
/// whatever the input order, in memory for 2 x `DEPTH_SAMPLE` reads at most.
#[derive(Default)]
struct Sample {
    kept: Vec<(u64, Vec<u8>)>,
    /// Hashes from here on cannot be among the smallest any more.
    bound: Option<u64>,
}

impl Sample {
    fn offer(&mut self, key: u64, seq: &[u8]) {
        if self.bound.is_some_and(|b| key >= b) {
            return;
        }
        self.kept.push((key, seq.to_vec()));
        if self.kept.len() == 2 * DEPTH_SAMPLE {
            self.kept.select_nth_unstable_by_key(DEPTH_SAMPLE, |&(k, _)| k);
            self.bound = Some(self.kept[DEPTH_SAMPLE].0);
            self.kept.truncate(DEPTH_SAMPLE);
        }
    }

    /// The sampled reads, smallest hash first.
    fn reads(&self) -> Vec<&[u8]> {
        let mut kept: Vec<&(u64, Vec<u8>)> = self.kept.iter().collect();
        kept.sort_unstable_by_key(|&&(k, _)| k);
        kept.truncate(DEPTH_SAMPLE);
        kept.into_iter().map(|(_, s)| s.as_slice()).collect()
    }
}

/// Sketch depth of the sampled reads (`target_depth`) and how many had one.
fn sample_depth(gate: &DepthGate, sample: &Sample) -> (Option<f64>, usize) {
    let mut reads: Vec<(f64, usize)> = sample
        .reads()
        .into_par_iter()
        .filter_map(|s| gate.median_depth(s).map(|d| (d, s.len())))
        .collect();
    if reads.is_empty() {
        return (None, 0);
    }
    let n = reads.len();
    (Some(target_depth(&mut reads, gate.organelle_min_depth)), n)
}

/// Reads this many times deeper than a target's own are a high-copy repeat, not the target: its
/// depth leaves them out (`target_depth`), and so does its extension (`run`).
const HIGH_COPY_FOLD: f64 = 10.0;

/// With the gate, a target that gained at least this share of the reads it gained the pass before
/// is still extending along its edges, not saturating, however small that is beside what it has:
/// a pass extends the bait by about a tenth of a read length past each edge, so a large genome
/// grows by about as many reads every pass. The Lophophora williamsii mitochondrion (2.2 Mb by
/// whole-library HiMT) gained 250-330 reads a pass and the 5% rule stopped it at pass 23, at
/// 418 kb. A target converging gains less every pass: on the four benchmark libraries every other
/// target stopped at the same pass with this rule as without it.
const STEADY_GROWTH: f64 = 0.5;

/// A target's depth from its reads' (sketch depth, length): their median, unless that is more
/// than `HIGH_COPY_FOLD` times the depth most of the target's sequence has, which is then taken.
/// A read stands for its length over its depth, so the reads of a high-copy repeat stand for
/// almost none: the first pass of the Lophophora williamsii plastid recruited 2,899 reads, most of
/// them a repeat at ~700,000x beside the plastid's own at ~440x; their median gave a 24 bp
/// "genome", extension took one read in 3,000 (none) and the output kept one read of 3,199 bp.
/// The sequence's depth is the span-weighted median of the organelle-deep reads (from
/// `organelle_min`), not of all: shallow reads stand for much sequence each. A mitochondrion with
/// plastid-deep reads among its own keeps its median (Medicago polymorpha: 3,212x, its sequence at
/// 1,045x), as does a plastid with its inverted repeat.
fn target_depth(reads: &mut [(f64, usize)], organelle_min: f64) -> f64 {
    reads.sort_by(|a, b| a.0.total_cmp(&b.0));
    let median = reads[reads.len() / 2].0;
    let deep: Vec<(f64, f64)> = reads
        .iter()
        .filter(|r| r.0 > 0.0 && r.0 >= organelle_min)
        .map(|&(d, len)| (d, len as f64 / d))
        .collect();
    let total: f64 = deep.iter().map(|r| r.1).sum();
    let mut acc = 0.0;
    for &(d, span) in &deep {
        acc += span;
        if 2.0 * acc >= total {
            return if median > HIGH_COPY_FOLD * d { d } else { median };
        }
    }
    median
}

/// Which of a target's recruited reads a pass writes to its pass file.
#[derive(Clone, Copy, PartialEq, Debug)]
enum Keep {
    /// None: the target no longer extends.
    Nothing,
    /// The reads new to it, the only ones extension reads.
    Fresh,
    /// All of them (a bootstrapped target recounts all its reads every pass).
    All,
}

/// A target's recruited reads of the current pass, in input order.
fn pass_path(out_dir: &Path, target: &str) -> PathBuf {
    out_dir.join(format!(".pass.{target}.fastq"))
}

/// Every record of a FASTQ written by `write_rec`, in order. A pass that recruited nothing
/// leaves an empty file, which the parser rejects, so it is read as no records.
pub(crate) fn for_each_rec(path: &Path, mut f: impl FnMut(Rec)) -> Result<()> {
    if std::fs::metadata(path)
        .with_context(|| format!("cannot read {}", path.display()))?
        .len()
        == 0
    {
        return Ok(());
    }
    let mut parser = needletail::parse_fastx_file(path)
        .with_context(|| format!("cannot read {}", path.display()))?;
    while let Some(record) = parser.next() {
        let record = record.with_context(|| format!("malformed record in {}", path.display()))?;
        f(Rec {
            id: record.id().to_vec(),
            seq: record.seq().into_owned(),
            qual: record.qual().map(<[u8]>::to_vec),
        });
    }
    Ok(())
}

/// Bases per batch. Batching by bases, not reads, keeps the work per batch the same for 150 bp
/// short reads and 15 kb long reads: 2,048 short reads were too little work to spread over 32
/// threads (8 threads classified faster than 32).
pub const BATCH_BASES: usize = 8 << 20;

/// Stream reads from all inputs in batches.
///
/// A gzip stream only decompresses serially, which bounds a pass over short reads (one reader
/// ran at `zcat` speed), so each input file gets its own reader thread: paired files decompress
/// in parallel. Their batches are merged round-robin (one batch of each open file in turn), so
/// the order depends only on the inputs and every result stays reproducible.
pub fn read_batches(
    inputs: &[PathBuf],
) -> (
    crossbeam_channel::Receiver<Vec<Rec>>,
    std::thread::JoinHandle<Result<()>>,
) {
    let (tx, rx) = crossbeam_channel::bounded::<Vec<Rec>>(8);
    let paths = inputs.to_vec();
    let merger = std::thread::spawn(move || -> Result<()> {
        let files: Vec<_> = paths
            .into_iter()
            .map(|path| {
                let (ftx, frx) = crossbeam_channel::bounded::<Vec<Rec>>(4);
                (frx, std::thread::spawn(move || read_file(&path, &ftx)))
            })
            .collect();
        let mut open = vec![true; files.len()];
        while open.iter().any(|&o| o) {
            for (i, (frx, _)) in files.iter().enumerate() {
                if !open[i] {
                    continue;
                }
                match frx.recv() {
                    Ok(batch) => {
                        if tx.send(batch).is_err() {
                            // the consumer stopped; dropping the channels stops the readers
                            return Ok(());
                        }
                    }
                    Err(_) => open[i] = false,
                }
            }
        }
        for (_, h) in files {
            h.join().expect("reader thread panicked")?;
        }
        Ok(())
    });
    (rx, merger)
}

/// Read one FASTA/FASTQ(.gz) file into batches of about `BATCH_BASES` bases.
fn read_file(path: &Path, tx: &crossbeam_channel::Sender<Vec<Rec>>) -> Result<()> {
    let mut parser = needletail::parse_fastx_file(path)
        .with_context(|| format!("cannot read {}", path.display()))?;
    let (mut batch, mut bases) = (Vec::new(), 0usize);
    while let Some(record) = parser.next() {
        let record = record.with_context(|| format!("malformed record in {}", path.display()))?;
        bases += record.num_bases();
        batch.push(Rec {
            id: record.id().to_vec(),
            seq: record.seq().into_owned(),
            qual: record.qual().map(<[u8]>::to_vec),
        });
        if bases >= BATCH_BASES {
            if tx.send(std::mem::take(&mut batch)).is_err() {
                return Ok(());
            }
            bases = 0;
        }
    }
    if !batch.is_empty() {
        let _ = tx.send(batch);
    }
    Ok(())
}

/// What a pass writes besides its counts and depth samples.
struct Writes<'a> {
    out_dir: &'a Path,
    numt_path: &'a Path,
    /// Reads already used to extend each target: the others are fresh.
    used: &'a [FxHashSet<u64>],
    /// `None`: every recruited read and every NUMT candidate, as a pass whose files may be the
    /// last. A replayed pass is never the last (one original-input pass follows for the
    /// qualities), so it writes per target only what extension reads, and no NUMT candidates:
    /// the cache gives those when a bootstrap needs them. A deep plastid saturates early and
    /// its gigabases are then no longer rewritten every pass, nor are the NUMT candidates
    /// (Nipponbare: 3 Gb a pass).
    keep: Option<&'a [Keep]>,
}

/// One pass over the inputs. Recruited reads are written, in input order, to each target's
/// pass file and NUMT candidates to `numt_path`, not kept: a plastid at 20,000x recruits
/// gigabases (Col-0 HiFi: 3.2 Gb, 14 GB of memory when held), more than a laptop has.
fn scan(
    inputs: &[PathBuf],
    bait: &BaitSet,
    p: &ClassifyParams,
    w: &Writes,
    mut replay: Option<&mut crate::replay::Replay>,
    delta: Option<&BaitSet>,
    final_outputs: bool,
) -> Result<ScanResult> {
    let (out_dir, numt_path) = (w.out_dir, w.numt_path);
    let replaying = replay.as_ref().is_some_and(|r| r.ready);
    let (rx, reader) = if replaying && !final_outputs {
        replay.as_ref().unwrap().batches()
    } else {
        crate::replay::source_batches(inputs)
    };
    let mut writers = bait
        .targets
        .iter()
        .map(|t| {
            Ok(BufWriter::with_capacity(
                1 << 20,
                File::create(pass_path(out_dir, &t.name))?,
            ))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut numt_w = BufWriter::with_capacity(1 << 20, File::create(numt_path)?);
    let expected_input = if final_outputs {
        replay.as_ref().map(|cache| cache.batches())
    } else {
        None
    };
    let n_targets = bait.targets.len();
    let mut out = ScanResult {
        recruited: vec![TargetCounts::default(); n_targets],
        fresh: vec![TargetCounts::default(); n_targets],
        all_sample: (0..n_targets).map(|_| Sample::default()).collect(),
        fresh_sample: (0..n_targets).map(|_| Sample::default()).collect(),
        numt: TargetCounts::default(),
        reads: 0,
        bases: 0,
        multi_target: 0,
        reclassified: 0,
    };
    for batch in rx {
        if let Some((expected, _)) = &expected_input {
            let cached = expected
                .recv()
                .context("original input has more batches than its recruitment cache")?;
            anyhow::ensure!(
                cached.len() == batch.len()
                    && cached
                        .iter()
                        .zip(&batch)
                        .all(|(a, b)| a.rec.id == b.rec.id && a.rec.seq == b.rec.seq),
                "original identifiers or sequences changed during recruitment replay"
            );
        }
        if !replaying {
            if let Some(cache) = replay.as_mut() {
                cache.store(&batch)?;
            }
        }
        let offset = out.reads as usize;
        let previous = if replaying {
            Some(
                replay
                    .as_ref()
                    .unwrap()
                    .calls
                    .get(offset..offset + batch.len())
                    .context("input has more records than the recruitment cache")?,
            )
        } else {
            None
        };
        out.reads += batch.len() as u64;
        out.bases += batch.iter().map(|r| r.rec.seq.len() as u64).sum::<u64>();
        let calls: Vec<(Call, crate::replay::Entry, bool)> = batch
            .into_par_iter()
            .enumerate()
            .map(|(i, rec)| {
                let changed = previous.is_none()
                    || (!final_outputs
                        && CanonicalKmers::new(&rec.rec.seq, bait.k)
                            .any(|(_, key)| delta.expect("replay has change set").get(key) != 0));
                let call = if changed {
                    classify(&rec.rec.seq, bait, p).0
                } else {
                    previous.unwrap()[i]
                };
                (call, rec, changed)
            })
            .collect();
        for (i, (call, entry, changed)) in calls.into_iter().enumerate() {
            out.reclassified += u64::from(changed);
            if let Some(cache) = replay.as_mut() {
                if replaying {
                    cache.calls[offset + i] = call;
                } else {
                    cache.calls.push(call);
                }
            }
            if call == Call::Other {
                continue;
            }
            let rec = entry.rec;
            let len = rec.seq.len() as u64;
            match call {
                Call::Organelle(mask) => {
                    if mask.count_ones() > 1 {
                        out.multi_target += 1;
                    }
                    let key = name_key(&rec.id);
                    for t in (0..8).filter(|t| mask & (1 << t) != 0) {
                        let fresh = !w.used[t].contains(&key);
                        let write = match w.keep.map(|k| k[t]) {
                            None | Some(Keep::All) => true,
                            Some(Keep::Fresh) => fresh,
                            Some(Keep::Nothing) => false,
                        };
                        if write {
                            write_rec(&mut writers[t], &rec)?;
                        }
                        out.recruited[t].reads += 1;
                        out.recruited[t].bases += len;
                        out.all_sample[t].offer(key, &rec.seq);
                        if fresh {
                            out.fresh[t].reads += 1;
                            out.fresh[t].bases += len;
                            out.fresh_sample[t].offer(key, &rec.seq);
                        }
                    }
                }
                Call::NumtCandidate(_) => {
                    if w.keep.is_none() {
                        write_rec(&mut numt_w, &rec)?;
                    }
                    out.numt.reads += 1;
                    out.numt.bases += len;
                }
                Call::Other => {}
            }
        }
    }
    reader.join().expect("reader thread panicked")?;
    if let Some((expected, reader)) = expected_input {
        anyhow::ensure!(
            expected.recv().is_err(),
            "original input has fewer batches than its recruitment cache"
        );
        reader.join().expect("cache reader panicked")?;
    }
    if let Some(cache) = replay {
        if !replaying {
            cache.finish()?;
        }
        if !final_outputs {
            cache.reclassified_reads.push(out.reclassified);
        }
    }
    for w in &mut writers {
        w.flush()?;
    }
    numt_w.flush()?;
    Ok(out)
}

/// 64-bit key of a read name (splitmix64 over FxHash).
fn name_key(id: &[u8]) -> u64 {
    let mut h = FxHasher::default();
    h.write(id);
    let mut z = h.finish();
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Uniform in [0, 1) from the read name and a salt; identical for identical inputs.
/// Salt of the subsample used for bait extension (independent of the output downsampling).
const EXTEND_SALT: u64 = 0xB417_E87E;

pub fn unit_hash(id: &[u8], salt: u64) -> f64 {
    let mut h = FxHasher::default();
    h.write(id);
    let mut z = h.finish() ^ salt.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    // splitmix64 finaliser: FxHash alone is not uniform enough for sampling.
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    (z >> 11) as f64 / (1u64 << 53) as f64
}

pub(crate) fn write_rec(w: &mut impl Write, r: &Rec) -> Result<()> {
    w.write_all(b"@")?;
    w.write_all(&r.id)?;
    w.write_all(b"\n")?;
    w.write_all(&r.seq)?;
    w.write_all(b"\n+\n")?;
    match &r.qual {
        Some(q) => w.write_all(q)?,
        // FASTA input: emit a neutral quality so every output is valid FASTQ.
        None => w.write_all(&vec![b'I'; r.seq.len()])?,
    }
    w.write_all(b"\n")?;
    Ok(())
}

#[derive(Debug, Clone, Serialize)]
pub struct TargetOutput {
    pub name: String,
    pub recruited_reads: u64,
    pub recruited_bases: u64,
    pub seed_length: u64,
    /// recruited bases / seed length; approximate when seeds come from a related species.
    pub estimated_depth: f64,
    /// Median sketch depth of the recruited reads (with the depth sketch): the organelle's
    /// read depth, whatever the seeds.
    pub read_depth: Option<f64>,
    /// Recruited bases / read depth: the genome size the reads stand for.
    pub implied_size: Option<u64>,
    /// Why the target stopped extending: seeded_from_reads, saturated (grew by less than
    /// --saturation), bait_complete (extension added nothing), size_cap (stood for more than
    /// `size_cap`), iteration_limit (still growing when the passes ran out) or single_pass.
    pub stop_reason: &'static str,
    /// False for iteration_limit and size_cap: the reads may be incomplete or hold others.
    pub converged: bool,
    pub sampling_fraction: f64,
    /// The rule that set `sampling_fraction` below 1: target_depth, max_read_depth or max_bases.
    pub sampled_by: Option<&'static str>,
    pub kept_reads: u64,
    pub kept_bases: u64,
    pub output: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct RecruitReport {
    pub tool: &'static str,
    pub version: &'static str,
    pub k: usize,
    pub params: ClassifyParams,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input_replay: Option<crate::replay::ReplayStats>,
    pub iterations_requested: usize,
    /// Passes actually run: iteration stops early once no target's bait grows.
    pub iterations_run: usize,
    pub extend_min_count: u32,
    pub extend_depth: f64,
    pub saturation: f64,
    pub target_depth: Option<f64>,
    pub max_read_depth: Option<f64>,
    pub max_bases: Option<u64>,
    pub salt: u64,
    pub inputs: Vec<String>,
    pub targets: Vec<crate::bait::Target>,
    pub iterations: Vec<IterationStats>,
    pub outputs: Vec<TargetOutput>,
    pub numt_candidates: u64,
    pub numt_output: String,
    pub source_candidates: serde_json::Value,
    pub bootstrap: bool,
    pub extend_gate: bool,
    /// Depth sketch of the bootstrap or extension gate, when one was built.
    pub bootstrap_nuclear_peak_depth: Option<u32>,
    pub bootstrap_organelle_min_depth: Option<f64>,
    pub depth_gate: Option<crate::seedfree::GateReport>,
    /// Each target's core depth (see `EXTEND_CORE_FLOOR`), once set.
    pub core_depth: Vec<Option<f64>>,
    /// What a caller should know about the result: a doubtful depth gate, a target stopped by
    /// `size_cap`, or one still growing when the passes ran out.
    pub warnings: Vec<String>,
    pub elapsed_seconds: f64,
    pub peak_rss_mb: Option<f64>,
}

pub struct RecruitConfig {
    /// Optional new directory for lossless, process-local input replay.
    pub replay_dir: Option<PathBuf>,
    pub inputs: Vec<PathBuf>,
    pub params: ClassifyParams,
    pub iterations: usize,
    pub extend_min_count: u32,
    /// Depth the new reads of a target are subsampled to before counting bait k-mers.
    pub extend_depth: f64,
    /// A target whose recruited reads grow by less than this fraction stops being extended.
    pub saturation: f64,
    pub target_depth: Option<f64>,
    pub salt: u64,
    pub out_dir: PathBuf,
    /// Bootstrap a target whose seed is too distant to recruit its reads (see
    /// `bootstrap_candidates`). Long reads only: the depth gate needs whole-read depth profiles.
    pub bootstrap: bool,
    /// Extend a target's bait only from reads that are organelle-deep along their whole length
    /// (`discover`'s test). Iterative extension follows any read the bait covers by
    /// `min_cover`; one bad seed sequence then walks into the nucleus pass after pass. In the
    /// Carex holdout library, extension to the fixed point added 75,858 reads to the plastid,
    /// 3% of them organelle-like; the 1,141 reads it added to the Buddleja mitochondrion were
    /// 100% organelle-like. Long reads only.
    ///
    /// The gate also keeps out reads much shallower than the target's own (`EXTEND_CORE_FLOOR`)
    /// or far deeper (`HIGH_COPY_FOLD`), and stops a target whose reads stand for more sequence
    /// than any such genome has (`size_cap`).
    pub extend_gate: bool,
    /// Downsample each target's output to about this read depth: the median sketch depth of
    /// its reads, so the recruited bases over it are the genome size they stand for. Long reads
    /// only (the depth sketch judges whole long reads).
    pub max_read_depth: Option<f64>,
    /// Downsample each target's output to about this many bases at most.
    pub max_bases: Option<u64>,
}

/// A bootstrap read must not be this much shallower than the target's own reads (when the
/// target has some): a mitochondrial read carrying a plastid insertion is ~10-30x shallower
/// than plastid reads, so it cannot bootstrap the plastid, and the reverse needs the other
/// target to be complete (see `bootstrap_candidates`).
const BOOTSTRAP_DEPTH_FLOOR: f64 = 0.5;
/// Fewest own reads for a depth floor.
const BOOTSTRAP_DEPTH_MIN_READS: usize = 10;
/// A bootstrap read must not be covered this much by another target's bait.
const BOOTSTRAP_MAX_OTHER_COVER: f64 = 0.5;

/// With the extension gate, a read extends a target only from this share of the target's core
/// depth up: the median depth of its reads in the first pass that recruited `CORE_MIN_READS`.
/// The gate's own threshold is relative to the nuclear depth, and that estimate can fail: on
/// three of the four 2026-10-21 benchmark libraries (Carex echinata, Lophophora williamsii,
/// Medicago polymorpha) it came out at 2x, the threshold at 6x, and every nuclear read
/// passed. A plastid at thousands of x walking into mitochondrial reads (tens to hundreds of
/// x, through a shared MTPT) is kept out the same way.
pub const EXTEND_CORE_FLOOR: f64 = 0.25;
/// Fewest reads with a depth that set a target's core depth.
const CORE_MIN_READS: usize = 20;

/// Most sequence a target named for an organelle can stand for (recruited bases over read
/// depth) before its extension stops: four times the largest land-plant plastome (Pelargonium,
/// 242 kb), and nearly twice the largest mitochondrial genome (Silene conica, 11.3 Mb). The
/// unguarded Carex echinata plastid stood for tens of megabases of nuclear reads.
pub fn size_cap(name: &str) -> Option<u64> {
    match name {
        "plastid" => Some(1_000_000),
        "mitochondrion" => Some(20_000_000),
        _ => None,
    }
}

/// The size a target's reads stand for and its cap, when the size is over it.
fn over_size_cap(name: &str, implied_size: Option<u64>) -> Option<(u64, u64)> {
    let (size, cap) = (implied_size?, size_cap(name)?);
    (size > cap).then_some((size, cap))
}

/// NUMT candidates of the last pass that bootstrap target `t`.
///
/// A seed from distant species matches a read only at conserved genes: with the land-plant
/// SeedDB, 874 of 1,916 Ajuga mitochondrial HiFi reads held a >= 1 kb seed block but fewer than
/// 2% were recruited, so iterative baiting never started, and with the mitochondrion seeded
/// alone, plastid reads with MTPT-like blocks took its place. Those partial reads are the NUMT
/// candidates, where nuclear insertions also land; depth separates the two: a mitochondrial read
/// is deep along its whole length, a nuclear one only at its organellar block. On the three
/// frozen holdout samples, `discover`'s organelle-like test passed 979 of 989 mitochondrial
/// candidates and 16 of 351 others. A candidate must also fit `t` better than any other target
/// and not be mostly covered by another target, and be at least `floor` deep; and it is only
/// used once every other target's bait has stopped growing, so their own reads are recruited,
/// not left as candidates.
fn bootstrap_candidates(
    numt_path: &Path,
    floor: f64,
    t: usize,
    bait: &BaitSet,
    p: &ClassifyParams,
    gate: &DepthGate,
    used: &FxHashSet<u64>,
) -> Result<(Vec<Vec<u8>>, Vec<u64>)> {
    let n = bait.targets.len();
    let (mut picked, mut keys) = (Vec::new(), Vec::new());
    for_each_rec(numt_path, |r| {
        let key = name_key(&r.id);
        if used.contains(&key) {
            return;
        }
        let h = target_hits(&r.seq, bait, p.block_gap);
        let fit = |u: usize| (h.block_cover[u], h.hits[u]);
        let len = r.seq.len() as f64;
        if h.longest[t] < p.numt_min_block
            || (0..n).any(|u| u != t && fit(u) >= fit(t))
            || (0..n).any(|u| u != t && h.block_cover[u] as f64 >= BOOTSTRAP_MAX_OTHER_COVER * len)
            || !gate.organelle_like(&r.seq)
            || gate.median_depth(&r.seq).unwrap_or(0.0) < floor
        {
            return;
        }
        keys.push(key);
        picked.push(r.seq);
    })?;
    Ok((picked, keys))
}

/// The depth sketch: from the recruitment cache when there is one, which saves its two passes
/// over the original input (each a full decompression).
fn build_gate(
    inputs: &[PathBuf],
    replay: Option<&crate::replay::Replay>,
    warnings: &mut Vec<String>,
) -> Result<DepthGate> {
    let gate = match replay.filter(|r| r.ready) {
        Some(cache) => DepthGate::long_reads_from(|| cache.sequences(), "recruitment cache")?,
        None => DepthGate::long_reads(inputs)?,
    };
    eprintln!(
        "[ovasm] depth gate: nuclear peak {}x, organelle threshold {:.0}x (sketch depth), from the {} in {:.0}s",
        gate.nuclear_peak_depth,
        gate.organelle_min_depth,
        gate.report.source,
        gate.report.seconds
    );
    // the nuclear peak stands for no nuclear genome: the organelle threshold is doubtful
    if gate.report.peak_window_bases < crate::seedfree::MIN_NUCLEAR_GENOME {
        warnings.push(format!(
            "the depth sketch's nuclear peak ({}x) stands for {:.1} Mb of sequence, too little \
             for a nuclear genome: the library may be too shallow for the nucleus to show, and \
             the organelle threshold ({:.0}x) is doubtful",
            gate.nuclear_peak_depth,
            gate.report.peak_window_bases as f64 / 1e6,
            gate.organelle_min_depth
        ));
    }
    Ok(gate)
}

/// The NUMT candidates of the last pass, from the cache: a replayed pass does not write them.
fn numt_from_cache(cache: &crate::replay::Replay, path: &Path) -> Result<()> {
    let mut w = BufWriter::with_capacity(1 << 20, File::create(path)?);
    let (rx, reader) = cache.batches();
    let mut i = 0;
    for batch in rx {
        for entry in batch {
            if matches!(cache.calls[i], Call::NumtCandidate(_)) {
                write_rec(&mut w, &entry.rec)?;
            }
            i += 1;
        }
    }
    reader.join().expect("cache reader panicked")?;
    w.flush()?;
    Ok(())
}

/// Depths of one pass's targets (`None` without the depth sketch).
struct PassDepths {
    read: Vec<Option<f64>>,
    fresh: Vec<Option<f64>>,
    /// Reads of the target that had a depth.
    measured: Vec<usize>,
    implied_size: Vec<Option<u64>>,
}

fn pass_depths(gate: Option<&DepthGate>, res: &ScanResult) -> PassDepths {
    let n = res.recruited.len();
    let mut d = PassDepths {
        read: vec![None; n],
        fresh: vec![None; n],
        measured: vec![0; n],
        implied_size: vec![None; n],
    };
    if let Some(g) = gate {
        for t in 0..n {
            (d.read[t], d.measured[t]) = sample_depth(g, &res.all_sample[t]);
            d.fresh[t] = sample_depth(g, &res.fresh_sample[t]).0;
            d.implied_size[t] = d.read[t]
                .filter(|&x| x > 0.0)
                .map(|x| (res.recruited[t].bases as f64 / x) as u64);
        }
    }
    d
}

pub fn run(mut bait: BaitSet, cfg: RecruitConfig) -> Result<RecruitReport> {
    let start = Instant::now();
    std::fs::create_dir_all(&cfg.out_dir)?;
    let mut iterations = Vec::new();
    let mut replay = cfg
        .replay_dir
        .as_deref()
        .map(crate::replay::Replay::new)
        .transpose()?;
    if replay.is_some() {
        bait.track_changes();
    }
    let mut result = None;
    let input_metadata = if replay.is_some() {
        cfg.inputs
            .iter()
            .map(|path| {
                let meta = std::fs::metadata(path)?;
                Ok((meta.len(), meta.modified()?))
            })
            .collect::<Result<Vec<_>>>()?
    } else {
        Vec::new()
    };
    let n_targets = bait.targets.len();
    let numt_path = cfg.out_dir.join("numt_candidates.fastq");
    // Reads already used to extend each target's bait: later passes count only new reads, so a
    // deep, already complete target (a plastid at 20,000x) is not re-counted every pass.
    // (64-bit hashes of the read names: millions of names, 8 bytes each instead of a string.)
    let mut used: Vec<FxHashSet<u64>> = vec![FxHashSet::default(); n_targets];
    // A target seeded from reads of the sample holds its whole genome already: extending it
    // only adds passes over the input (each a full decompression) for nothing.
    let mut saturated: Vec<bool> = bait.targets.iter().map(|t| t.from_reads).collect();
    // Why each target no longer extends, once it does not.
    let mut stopped: Vec<Option<&'static str>> = bait
        .targets
        .iter()
        .map(|t| t.from_reads.then_some("seeded_from_reads"))
        .collect();
    let mut previous = vec![0u64; n_targets];
    // Reads each target gained in the last pass (`STEADY_GROWTH`).
    let mut last_growth = vec![0u64; n_targets];
    // Depth sketch, built after the first pass when the gate or the downsampling needs it, or
    // the first time a target needs bootstrapping.
    let mut gate: Option<DepthGate> = None;
    let wants_gate = cfg.extend_gate || cfg.max_read_depth.is_some();
    let mut core_depth: Vec<Option<f64>> = vec![None; n_targets];
    let mut warnings: Vec<String> = Vec::new();
    let mut boot_used: Vec<FxHashSet<u64>> = vec![FxHashSet::default(); n_targets];
    for it in 0..cfg.iterations.max(1) {
        let t0 = Instant::now();
        let bait_kmers = bait.len();
        bait.prepare_lookup();
        let replaying = replay.as_ref().is_some_and(|r| r.ready);
        let delta = replay
            .as_ref()
            .filter(|r| r.ready)
            .map(|_| bait.take_changes());
        let keep: Vec<Keep> = (0..n_targets)
            .map(|t| {
                if saturated[t] {
                    Keep::Nothing
                } else if !boot_used[t].is_empty() {
                    Keep::All
                } else {
                    Keep::Fresh
                }
            })
            .collect();
        let res = scan(
            &cfg.inputs,
            &bait,
            &cfg.params,
            &Writes {
                out_dir: &cfg.out_dir,
                numt_path: &numt_path,
                used: &used,
                keep: replaying.then_some(keep.as_slice()),
            },
            replay.as_mut(),
            delta.as_ref(),
            false,
        )?;
        let last = it + 1 == cfg.iterations.max(1);
        if gate.is_none() && wants_gate && (!last || cfg.max_read_depth.is_some()) {
            gate = Some(build_gate(&cfg.inputs, replay.as_ref(), &mut warnings)?);
        }
        let mut depths = pass_depths(gate.as_ref(), &res);
        for t in 0..n_targets {
            if core_depth[t].is_none() && depths.measured[t] >= CORE_MIN_READS {
                core_depth[t] = depths.read[t];
            }
        }
        if cfg.extend_gate && !last {
            for t in 0..n_targets {
                let Some((size, cap)) =
                    over_size_cap(&bait.targets[t].name, depths.implied_size[t])
                else {
                    continue;
                };
                if saturated[t] {
                    continue;
                }
                saturated[t] = true;
                stopped[t] = Some("size_cap");
                warnings.push(format!(
                    "{}: after pass {} its {} reads ({:.1} Mb at {:.0}x) stand for {:.1} Mb, more \
                     than any {} ({:.1} Mb); extension stopped there",
                    bait.targets[t].name,
                    it + 1,
                    res.recruited[t].reads,
                    res.recruited[t].bases as f64 / 1e6,
                    depths.read[t].unwrap_or(0.0),
                    size as f64 / 1e6,
                    bait.targets[t].name,
                    cap as f64 / 1e6
                ));
            }
        }
        let mut added = vec![0; n_targets];
        let mut extended_from = vec![0u64; n_targets];
        let mut gated_out = vec![0u64; n_targets];
        let mut shallow_out = vec![0u64; n_targets];
        let mut deep_out = vec![0u64; n_targets];
        if !last {
            for t in 0..n_targets {
                if saturated[t] {
                    continue;
                }
                // A bootstrapped target counts all its reads every pass: it grows from scattered
                // genes at mitochondrial depth, where a pass adds a handful of reads, too few for
                // a k-mer to reach `extend_min_count` among them alone (Alisma, ~14x: the bait
                // stalled at 258 of 316 reads, 55 of the rest partial candidates at its edges).
                // For the same reason it grows by a few percent a pass for many passes, so it
                // ends when its bait stops growing, not by `saturation`.
                let recount = !boot_used[t].is_empty();
                let now = res.recruited[t].reads;
                let grown = now.saturating_sub(previous[t]);
                let steady = cfg.extend_gate && grown as f64 >= STEADY_GROWTH * last_growth[t] as f64;
                if !recount
                    && previous[t] > 0
                    && (grown as f64) < cfg.saturation * previous[t] as f64
                    && !steady
                {
                    saturated[t] = true;
                    stopped[t] = Some("saturated");
                }
                last_growth[t] = grown;
                previous[t] = now;
                if saturated[t] {
                    continue;
                }
                // New reads of this pass, from its file: their depth sets the subsample
                // fraction (at most `extend_depth`). The genome size is the one the target's reads
                // stand for when the depth sketch gives it, else the seed length: a multi-species
                // library's summed length (6.4 Mb of plastids) is no genome size, and the
                // Medicago polymorpha plastid (~10,000x) extended from 128,386 reads.
                let path = pass_path(&cfg.out_dir, &bait.targets[t].name);
                let considered = if recount {
                    &res.recruited[t]
                } else {
                    &res.fresh[t]
                };
                let genome = depths.implied_size[t]
                    .filter(|_| cfg.extend_gate)
                    .unwrap_or(bait.targets[t].seed_length)
                    .max(1) as f64;
                let depth = considered.bases as f64 / genome;
                let fraction = if depth > cfg.extend_depth {
                    cfg.extend_depth / depth
                } else {
                    1.0
                };
                let extend_gate = gate.as_ref().filter(|_| cfg.extend_gate);
                let floor = core_depth[t].map(|c| EXTEND_CORE_FLOOR * c);
                let ceiling = core_depth[t].map(|c| HIGH_COPY_FOLD * c);
                let (mut fresh, mut picked) = (Vec::new(), Vec::<Vec<u8>>::new());
                for_each_rec(&path, |r| {
                    let key = name_key(&r.id);
                    if !recount && used[t].contains(&key) {
                        return;
                    }
                    fresh.push(key);
                    if fraction < 1.0 && unit_hash(&r.id, EXTEND_SALT) >= fraction {
                        return;
                    }
                    match extend_gate.map(|g| g.profile(&r.seq)) {
                        None => picked.push(r.seq),
                        Some(Some(d)) if d.organelle_like => {
                            if floor.is_some_and(|f| d.median < f) {
                                shallow_out[t] += 1;
                            } else if ceiling.is_some_and(|c| d.median > c) {
                                deep_out[t] += 1;
                            } else {
                                picked.push(r.seq);
                            }
                        }
                        Some(_) => gated_out[t] += 1,
                    }
                })?;
                extended_from[t] = picked.len() as u64;
                added[t] = bait.extend_from_reads(
                    t,
                    picked.iter().map(Vec::as_slice),
                    cfg.extend_min_count,
                );
                used[t].extend(fresh);
            }
        }
        let mut bootstrapped = vec![0u64; n_targets];
        let mut numt_ready = !replaying;
        if cfg.bootstrap && !last {
            for t in 0..n_targets {
                // complete: saturated, or its bait gained nothing from this pass
                let others_complete = (0..n_targets)
                    .filter(|&u| u != t)
                    .all(|u| saturated[u] || added[u] == 0);
                if bait.targets[t].from_reads
                    || added[t] > 0
                    || !others_complete
                    || !boot_used[t].is_empty()
                    || stopped[t] == Some("size_cap")
                {
                    continue;
                }
                if gate.is_none() {
                    gate = Some(build_gate(&cfg.inputs, replay.as_ref(), &mut warnings)?);
                    depths = pass_depths(gate.as_ref(), &res);
                }
                if !numt_ready {
                    numt_from_cache(replay.as_ref().expect("a replayed pass"), &numt_path)?;
                    numt_ready = true;
                }
                let g = gate.as_ref().expect("built above");
                let floor = if depths.measured[t] >= BOOTSTRAP_DEPTH_MIN_READS {
                    BOOTSTRAP_DEPTH_FLOOR * depths.read[t].unwrap_or(0.0)
                } else {
                    0.0
                };
                let (picked, keys) = bootstrap_candidates(
                    &numt_path,
                    floor,
                    t,
                    &bait,
                    &cfg.params,
                    g,
                    &boot_used[t],
                )?;
                // Only a seed that failed is bootstrapped: one that recruits more reads than are
                // left as deep partial candidates has its genome already. Its leftovers are
                // other deep sequence (Antitrichia: 7 nuclear repeat reads beside the complete
                // 592-read mitochondrion, from which extension ran away into the nucleus).
                if picked.is_empty() || picked.len() as u64 <= res.recruited[t].reads {
                    continue;
                }
                for seq in &picked {
                    added[t] += bait.add_sequence_counted(t, seq);
                }
                bootstrapped[t] = picked.len() as u64;
                boot_used[t].extend(keys);
                // its bait changed: extend it again from what the next pass recruits, and its
                // core depth is that of the reads it now recruits
                saturated[t] = false;
                stopped[t] = None;
                core_depth[t] = None;
            }
        }
        iterations.push(IterationStats {
            iteration: it + 1,
            seconds: t0.elapsed().as_secs_f64(),
            reads_scanned: res.reads,
            bases_scanned: res.bases,
            bait_kmers,
            recruited: res.recruited.clone(),
            numt_candidates: res.numt.clone(),
            multi_target_reads: res.multi_target,
            bait_kmers_added: added.clone(),
            extended_from_reads: extended_from,
            saturated: saturated.clone(),
            bootstrapped_from_reads: bootstrapped,
            extension_gated_out: gated_out,
            extension_shallow_out: shallow_out,
            extension_deep_out: deep_out,
            read_depth: depths.read,
            fresh_read_depth: depths.fresh,
            implied_size: depths.implied_size,
        });
        result = Some(res);
        // Nothing new to bait with: this pass already used the final bait set.
        if !last && added.iter().all(|&a| a == 0) {
            break;
        }
    }
    // Targets still extending when the passes ran out had not converged.
    let ran_out = iterations.len() == cfg.iterations.max(1) && cfg.iterations > 1;
    for t in 0..n_targets {
        if stopped[t].is_some() {
            continue;
        }
        let grew = iterations.len() >= 2 && iterations[iterations.len() - 2].bait_kmers_added[t] > 0;
        stopped[t] = Some(if cfg.iterations <= 1 {
            "single_pass"
        } else if ran_out && grew {
            "iteration_limit"
        } else {
            "bait_complete"
        });
        if stopped[t] == Some("iteration_limit") {
            let n = iterations.len();
            let (a, b) = (
                iterations[n - 2].recruited[t].reads,
                iterations[n - 1].recruited[t].reads,
            );
            warnings.push(format!(
                "{}: still growing after the last of {} passes ({} reads, +{:.1}% in the last \
                 pass): it may be incomplete; more passes would extend it further \
                 (--iterations; `run --recruit-iterations`)",
                bait.targets[t].name,
                n,
                b,
                100.0 * b.saturating_sub(a) as f64 / a.max(1) as f64
            ));
        }
    }
    let mut res = result.expect("at least one iteration");
    // All intermediate bait operations use only identifiers and sequence. Restore
    // actual qualities in one original-input pass after the final calls are known.
    if replay.is_some() && iterations.len() > 1 {
        for (path, &(size, modified)) in cfg.inputs.iter().zip(&input_metadata) {
            let meta = std::fs::metadata(path)?;
            anyhow::ensure!(
                meta.len() == size && meta.modified()? == modified,
                "original input metadata changed during recruitment replay: {}",
                path.display()
            );
        }
        let final_res = scan(
            &cfg.inputs,
            &bait,
            &cfg.params,
            &Writes {
                out_dir: &cfg.out_dir,
                numt_path: &numt_path,
                used: &used,
                keep: None,
            },
            replay.as_mut(),
            None,
            true,
        )?;
        anyhow::ensure!(
            final_res.reads == res.reads && final_res.bases == res.bases,
            "original input changed during recruitment replay"
        );
        res = final_res;
    }

    let source_candidates =
        crate::retained::recruitment(&cfg.out_dir, &numt_path, &bait, &cfg.params)?;
    let final_depths = pass_depths(gate.as_ref(), &res);
    // Final outputs from the last pass's files, downsampled to `target_depth`, `max_read_depth`
    // and `max_bases`, whichever keeps fewest.
    let mut outputs = Vec::new();
    for (t, counts) in res.recruited.iter().enumerate() {
        let target = &bait.targets[t];
        let depth = if target.seed_length > 0 {
            counts.bases as f64 / target.seed_length as f64
        } else {
            0.0
        };
        let mut fraction = 1.0;
        let mut sampled_by = None;
        let mut cap = |f: f64, rule: &'static str| {
            if f < fraction {
                fraction = f;
                sampled_by = Some(rule);
            }
        };
        if let Some(d) = cfg.target_depth.filter(|&d| depth > d) {
            cap(d / depth, "target_depth");
        }
        if let (Some(max), Some(d)) = (cfg.max_read_depth, final_depths.read[t]) {
            if d > max {
                cap(max / d, "max_read_depth");
            }
        }
        if let Some(max) = cfg.max_bases.filter(|&m| counts.bases > m) {
            cap(max as f64 / counts.bases as f64, "max_bases");
        }
        let path = cfg.out_dir.join(format!("{}.fastq", target.name));
        let pass = pass_path(&cfg.out_dir, &target.name);
        let mut w = BufWriter::with_capacity(1 << 20, File::create(&path)?);
        let (mut kept_reads, mut kept_bases) = (0u64, 0u64);
        let mut failed = Ok(());
        for_each_rec(&pass, |r| {
            if failed.is_ok() && (fraction >= 1.0 || unit_hash(&r.id, cfg.salt) < fraction) {
                kept_reads += 1;
                kept_bases += r.seq.len() as u64;
                failed = write_rec(&mut w, &r);
            }
        })?;
        failed?;
        w.flush()?;
        std::fs::remove_file(&pass)?;
        let stop_reason = stopped[t].unwrap_or("bait_complete");
        outputs.push(TargetOutput {
            name: target.name.clone(),
            recruited_reads: counts.reads,
            recruited_bases: counts.bases,
            seed_length: target.seed_length,
            estimated_depth: depth,
            read_depth: final_depths.read[t],
            implied_size: final_depths.implied_size[t],
            stop_reason,
            converged: !matches!(stop_reason, "iteration_limit" | "size_cap"),
            sampling_fraction: fraction,
            sampled_by,
            kept_reads,
            kept_bases,
            output: path.display().to_string(),
        });
    }
    for w in &warnings {
        eprintln!("[ovasm] recruit warning: {w}");
    }

    let report = RecruitReport {
        tool: "ovasm recruit",
        version: env!("CARGO_PKG_VERSION"),
        k: bait.k,
        params: cfg.params.clone(),
        input_replay: replay
            .as_ref()
            .map(|r| r.stats(iterations.len()))
            .transpose()?,
        iterations_requested: cfg.iterations,
        iterations_run: iterations.len(),
        extend_min_count: cfg.extend_min_count,
        extend_depth: cfg.extend_depth,
        saturation: cfg.saturation,
        target_depth: cfg.target_depth,
        max_read_depth: cfg.max_read_depth,
        max_bases: cfg.max_bases,
        salt: cfg.salt,
        inputs: cfg.inputs.iter().map(|p| p.display().to_string()).collect(),
        targets: bait.targets.clone(),
        iterations,
        outputs,
        numt_candidates: res.numt.reads,
        numt_output: numt_path.display().to_string(),
        source_candidates,
        bootstrap: cfg.bootstrap,
        extend_gate: cfg.extend_gate,
        bootstrap_nuclear_peak_depth: gate.as_ref().map(|g| g.nuclear_peak_depth),
        bootstrap_organelle_min_depth: gate.as_ref().map(|g| g.organelle_min_depth),
        depth_gate: gate.as_ref().map(|g| g.report.clone()),
        core_depth,
        warnings,
        elapsed_seconds: start.elapsed().as_secs_f64(),
        peak_rss_mb: peak_rss_mb(),
    };
    std::fs::write(
        cfg.out_dir.join("recruit.json"),
        serde_json::to_string_pretty(&report)?,
    )?;
    if let Some(cache) = replay {
        cache.remove()?;
    }
    Ok(report)
}

/// Peak resident memory from /proc (Linux); `None` elsewhere.
pub fn peak_rss_mb() -> Option<f64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    let line = status.lines().find(|l| l.starts_with("VmHWM:"))?;
    let kb: f64 = line.split_whitespace().nth(1)?.parse().ok()?;
    Some(kb / 1024.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unchanged_keys_preserve_calls_across_competing_target_growth() {
        let a = pseudo_random(2000, 601);
        let b = pseudo_random(2000, 602);
        let mut bait = bait_from(&a[..400], 21);
        bait.target("plastid").unwrap();
        bait.track_changes();
        let mut flanked = pseudo_random(1500, 603);
        flanked.extend_from_slice(&a);
        flanked.extend(pseudo_random(1500, 604));
        let reads = [&a[..], &b[..], &a[..300], &flanked[..], &b[400..800]];
        for (target, sequence) in [(0, &a[..]), (1, &a[..]), (1, &b[..])] {
            let old: Vec<_> = reads
                .iter()
                .map(|r| classify(r, &bait, &params()).0)
                .collect();
            bait.add_sequence(target, sequence);
            let delta = bait.take_changes();
            for (read, previous) in reads.iter().zip(old) {
                let changed = CanonicalKmers::new(read, 21).any(|(_, key)| delta.get(key) != 0);
                let replayed = if changed {
                    classify(read, &bait, &params()).0
                } else {
                    previous
                };
                assert_eq!(replayed, classify(read, &bait, &params()).0);
            }
        }
    }

    fn bait_from(seq: &[u8], k: usize) -> BaitSet {
        let mut b = BaitSet::new(k);
        let t = b.target("mito").unwrap();
        b.add_sequence(t, seq);
        b.targets[t].seed_length = seq.len() as u64;
        b
    }

    fn params() -> ClassifyParams {
        ClassifyParams {
            min_cover: 0.8,
            min_hit_frac: 0.2,
            numt_min_block: 40,
            block_gap: 10,
        }
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

    #[test]
    fn organelle_read_is_recruited() {
        let genome = pseudo_random(2000, 7);
        let bait = bait_from(&genome, 15);
        let (call, score) = classify(&genome[100..600], &bait, &params());
        assert_eq!(call, Call::Organelle(1));
        assert!(score.cover_frac > 0.99);
    }

    #[test]
    fn presence_index_preserves_calls_and_scores_across_read_contexts() {
        for (k, gap, hit, numt) in [
            (21, 50, 0.2, 1000),
            (15, 200, 0.05, 2000),
            (21, 40, 0.2, usize::MAX),
            (31, 0, 0.2, 40),
        ] {
            let shared = pseudo_random(2500, 39);
            let mut exact = BaitSet::new(k);
            let mut reads = vec![Vec::new(), vec![b'N'; 300], pseudo_random(5000, 99)];
            for i in 0..8 {
                let t = exact.target(&format!("target{i}")).unwrap();
                let genome = pseudo_random(5000, 7 + i);
                exact.add_sequence(t, &genome);
                exact.add_sequence(t, &shared);
                reads.push(genome.clone());
                let mut flanked = pseudo_random(2500, 101);
                flanked.extend_from_slice(&genome[1000..3500]);
                flanked.extend(pseudo_random(2500, 103));
                reads.push(flanked);
                let mut with_gap = genome;
                with_gap[2000..2500].fill(b'N');
                reads.push(with_gap);
            }
            reads.push(shared);
            let reverse_reads: Vec<Vec<u8>> = reads
                .iter()
                .map(|r| {
                    r.iter()
                        .rev()
                        .map(|b| match b {
                            b'A' => b'T',
                            b'C' => b'G',
                            b'G' => b'C',
                            b'T' => b'A',
                            _ => b'N',
                        })
                        .collect()
                })
                .collect();
            reads.extend(reverse_reads);
            let mut indexed = exact.clone();
            indexed.prepare_lookup();
            let p = ClassifyParams {
                min_cover: 0.8,
                min_hit_frac: hit,
                numt_min_block: numt,
                block_gap: gap,
            };
            for read in &reads {
                let (a, x) = classify(read, &exact, &p);
                let (b, y) = classify(read, &indexed, &p);
                assert_eq!(a, b);
                assert_eq!(x.cover_frac.to_bits(), y.cover_frac.to_bits());
                assert_eq!(x.hit_frac.to_bits(), y.hit_frac.to_bits());
                assert_eq!(x.longest_block, y.longest_block);
            }
        }
    }

    #[test]
    fn nuclear_read_is_ignored() {
        let genome = pseudo_random(2000, 7);
        let bait = bait_from(&genome, 15);
        let nuclear = pseudo_random(500, 99);
        assert_eq!(classify(&nuclear, &bait, &params()).0, Call::Other);
    }

    #[test]
    fn organelle_segment_with_foreign_flanks_is_a_numt_candidate() {
        let genome = pseudo_random(2000, 7);
        let bait = bait_from(&genome, 15);
        let mut read = pseudo_random(200, 11);
        read.extend_from_slice(&genome[500..600]);
        read.extend(pseudo_random(200, 13));
        let (call, score) = classify(&read, &bait, &params());
        assert_eq!(call, Call::NumtCandidate(0));
        assert!(score.longest_block >= 90);
    }

    #[test]
    fn error_clusters_are_bridged_but_not_without_gap_tolerance() {
        // Isolated errors cost one base of block coverage; what breaks coverage on raw ONT is
        // error-dense stretches where no error-free k-mer fits. Every 100 bp here carries six
        // errors 7 bp apart (positions 40..75), leaving a ~36 bp hit-free hole.
        let genome = pseudo_random(3000, 7);
        let bait = bait_from(&genome, 15);
        let mut noisy = genome[500..1500].to_vec();
        for seg in (0..noisy.len()).step_by(100) {
            for off in [40, 47, 54, 61, 68, 75] {
                if let Some(b) = noisy.get_mut(seg + off) {
                    *b = if *b == b'A' { b'C' } else { b'A' };
                }
            }
        }
        let bridged = ClassifyParams {
            min_cover: 0.8,
            min_hit_frac: 0.05,
            numt_min_block: 5000,
            block_gap: 50,
        };
        assert_eq!(classify(&noisy, &bait, &bridged).0, Call::Organelle(1));
        let strict = ClassifyParams {
            block_gap: 0,
            ..bridged
        };
        assert_eq!(classify(&noisy, &bait, &strict).0, Call::Other);
    }

    #[test]
    fn exact_hits_break_coverage_ties_and_identical_copies_go_to_both() {
        let mito = pseudo_random(2000, 7);
        let mut plastid = pseudo_random(2000, 21);
        plastid[800..1000].copy_from_slice(&mito[800..1000]); // MTPT-like shared segment
        let snv = 900; // the plastid copy differs by one SNV
        plastid[snv] = if mito[snv] == b'A' { b'C' } else { b'A' };
        let mut bait = BaitSet::new(15);
        let (m, c) = (
            bait.target("mito").unwrap(),
            bait.target("plastid").unwrap(),
        );
        bait.add_sequence(m, &mito);
        bait.add_sequence(c, &plastid);
        let p = ClassifyParams {
            min_cover: 0.8,
            min_hit_frac: 0.2,
            numt_min_block: 5000,
            block_gap: 40,
        };
        // A mito read spanning the SNV: both targets are fully covered once the gap is bridged,
        // but the mito copy has more exact hits.
        assert_eq!(
            classify(&mito[850..950], &bait, &p).0,
            Call::Organelle(1 << m)
        );
        // A read from the identical part fits both equally and is assigned to both.
        assert_eq!(
            classify(&mito[810..880], &bait, &p).0,
            Call::Organelle((1 << m) | (1 << c))
        );
    }

    #[test]
    fn an_empty_pass_file_has_no_records() {
        let dir = std::env::temp_dir().join(format!("ovasm-empty-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(".pass.mitochondrion.fastq");
        std::fs::write(&path, b"").unwrap();
        let mut n = 0;
        for_each_rec(&path, |_| n += 1).unwrap();
        assert_eq!(n, 0);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_distant_seed_is_bootstrapped_from_deep_partial_reads_only() {
        // The mitochondrion's seed is one 1.5 kb gene, so no 3 kb read reaches 80% cover; the
        // plastid seed is complete. A nuclear insertion of the same gene (NUMT) at nuclear depth
        // must not bootstrap the mitochondrion, nor may plastid reads.
        let mito = pseudo_random(12_000, 7);
        let plastid = pseudo_random(20_000, 21);
        let mut nuclear = pseudo_random(300_000, 33);
        let gene = mito[5000..6500].to_vec();
        nuclear[150_000..151_500].copy_from_slice(&gene);
        let circular_reads = |g: &[u8], step: usize, tag: &str| -> Vec<(String, Vec<u8>)> {
            let circ = [g, &g[..3000]].concat();
            (0..g.len())
                .step_by(step)
                .enumerate()
                .map(|(i, st)| (format!("{tag}{i}"), circ[st..st + 3000].to_vec()))
                .collect()
        };
        let mut reads = circular_reads(&mito, 30, "m"); // 100x
        reads.extend(circular_reads(&plastid, 10, "p")); // 300x
        reads.extend(
            (0..nuclear.len() - 3000)
                .step_by(1500)
                .enumerate()
                .map(|(i, st)| (format!("n{i}"), nuclear[st..st + 3000].to_vec())),
        ); // 2x, three reads span the NUMT
        let dir = std::env::temp_dir().join(format!("ovasm-boot-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("reads.fa");
        let text: String = reads
            .iter()
            .map(|(id, s)| format!(">{id}\n{}\n", std::str::from_utf8(s).unwrap()))
            .collect();
        std::fs::write(&fa, text).unwrap();
        let mut bait = BaitSet::new(21);
        let (m, c) = (
            bait.target("mitochondrion").unwrap(),
            bait.target("plastid").unwrap(),
        );
        bait.add_sequence(m, &gene);
        bait.add_sequence(c, &plastid);
        bait.targets[m].seed_length = mito.len() as u64;
        bait.targets[c].seed_length = plastid.len() as u64;
        let cfg = |bootstrap| RecruitConfig {
            replay_dir: None,
            inputs: vec![fa.clone()],
            params: ClassifyParams {
                min_cover: 0.8,
                min_hit_frac: 0.2,
                numt_min_block: 1000,
                block_gap: 50,
            },
            iterations: 40,
            extend_min_count: 3,
            extend_depth: 200.0,
            saturation: 0.05,
            target_depth: None,
            salt: 0,
            out_dir: dir.join(format!("boot{bootstrap}")),
            bootstrap,
            extend_gate: false,
            max_read_depth: None,
            max_bases: None,
        };
        let ids = |path: &str| {
            let mut v = Vec::new();
            for_each_rec(Path::new(path), |r| {
                v.push(String::from_utf8(r.id).unwrap())
            })
            .unwrap();
            v
        };
        let plain = run(bait.clone(), cfg(false)).unwrap();
        assert_eq!(
            plain.outputs[m].recruited_reads, 0,
            "the gene alone recruits nothing"
        );
        // replayed passes write no NUMT candidates: the bootstrap takes them from the cache
        let mut replayed = cfg(true);
        replayed.out_dir = dir.join("boot-replayed");
        replayed.replay_dir = Some(dir.join("boot-cache"));
        let cached = run(bait.clone(), replayed).unwrap();
        let boot = run(bait, cfg(true)).unwrap();
        assert_same_recruitment(&boot, &cached);
        let got = ids(&boot.outputs[m].output);
        assert_eq!(got.len(), 400, "every mitochondrial read");
        assert!(got.iter().all(|id| id.starts_with('m')), "{got:?}");
        assert!(boot
            .iterations
            .iter()
            .any(|it| it.bootstrapped_from_reads[m] > 0));
        assert!(boot
            .iterations
            .iter()
            .all(|it| it.bootstrapped_from_reads[c] == 0));
        let pt = ids(&boot.outputs[c].output);
        assert_eq!(pt.len(), 2000);
        assert!(pt.iter().all(|id| id.starts_with('p')));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn the_extension_gate_keeps_a_stray_seed_out_of_the_nucleus() {
        // The seed holds the 12 kb mitochondrion (120x) plus 20 kb of a nuclear genome at 30x,
        // deep enough for its k-mers to pass extension's depth-relative floor. Extending to the
        // fixed point walks along that genome from the stray stretch; with the gate, its reads
        // (the dominant depth in the sketch, so not organelle-deep) never extend the bait.
        let mito = pseudo_random(12_000, 7);
        let nuclear = pseudo_random(200_000, 33);
        let circ = [mito.clone(), mito[..3000].to_vec()].concat();
        let mut text = String::new();
        for (i, st) in (0..mito.len()).step_by(25).enumerate() {
            text.push_str(&format!(
                ">m{i}\n{}\n",
                std::str::from_utf8(&circ[st..st + 3000]).unwrap()
            ));
        }
        for (i, st) in (0..nuclear.len() - 3000).step_by(100).enumerate() {
            text.push_str(&format!(
                ">n{i}\n{}\n",
                std::str::from_utf8(&nuclear[st..st + 3000]).unwrap()
            ));
        }
        let dir = std::env::temp_dir().join(format!("ovasm-gate-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("reads.fa");
        std::fs::write(&fa, text).unwrap();
        let mut bait = bait_from(&mito, 21);
        bait.add_sequence(0, &nuclear[100_000..120_000]);
        let cfg = |extend_gate| RecruitConfig {
            replay_dir: None,
            inputs: vec![fa.clone()],
            params: ClassifyParams {
                min_cover: 0.8,
                min_hit_frac: 0.2,
                numt_min_block: 1000,
                block_gap: 50,
            },
            iterations: 30,
            extend_min_count: 3,
            extend_depth: 200.0,
            saturation: 0.0,
            target_depth: None,
            salt: 0,
            out_dir: dir.join(format!("gate{extend_gate}")),
            bootstrap: false,
            extend_gate,
            max_read_depth: None,
            max_bases: None,
        };
        let nuclear_reads = |r: &RecruitReport| {
            let mut n = 0;
            for_each_rec(Path::new(&r.outputs[0].output), |rec| {
                n += usize::from(rec.id.starts_with(b"n"));
            })
            .unwrap();
            n
        };
        let plain = run(bait.clone(), cfg(false)).unwrap();
        // the depth sketch from the cache is the one from the input
        let mut replayed = cfg(true);
        replayed.out_dir = dir.join("gate-replayed");
        replayed.replay_dir = Some(dir.join("gate-cache"));
        let cached = run(bait.clone(), replayed).unwrap();
        let gated = run(bait, cfg(true)).unwrap();
        assert_same_recruitment(&gated, &cached);
        assert_eq!(gated.depth_gate.as_ref().unwrap().source, "inputs");
        assert_eq!(
            cached.depth_gate.as_ref().unwrap().source,
            "recruitment cache"
        );
        assert!(
            nuclear_reads(&plain) > nuclear_reads(&gated) + 100,
            "extension walks into the nucleus: {}",
            nuclear_reads(&plain)
        );
        assert!(
            nuclear_reads(&gated) <= 200,
            "only reads the stray seed covers: {}",
            nuclear_reads(&gated)
        );
        assert_eq!(
            gated.outputs[0].recruited_reads as usize - nuclear_reads(&gated),
            480
        );
        assert!(gated
            .iterations
            .iter()
            .any(|it| it.extension_gated_out[0] > 0));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn unit_hash_is_deterministic_uniform_and_salted() {
        let ids: Vec<Vec<u8>> = (0..20000)
            .map(|i| format!("read{i}").into_bytes())
            .collect();
        let a: Vec<f64> = ids.iter().map(|i| unit_hash(i, 0)).collect();
        let b: Vec<f64> = ids.iter().map(|i| unit_hash(i, 0)).collect();
        assert_eq!(a, b);
        let kept = a.iter().filter(|&&u| u < 0.25).count() as f64 / a.len() as f64;
        assert!((kept - 0.25).abs() < 0.02, "fraction {kept}");
        let c: Vec<f64> = ids.iter().map(|i| unit_hash(i, 1)).collect();
        assert_ne!(a, c);
    }

    #[test]
    fn files_are_read_in_parallel_and_merged_round_robin_by_bases() {
        // 1 MiB reads: eight fill a batch. A: 20 reads -> 8, 8, 4; B: 10 reads -> 8, 2.
        let dir = std::env::temp_dir().join(format!("ovasm-batches-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let body = "ACGT".repeat(1 << 18);
        let write = |name: &str, n: usize| {
            let p = dir.join(format!("{name}.fa"));
            let text: String = (0..n).map(|i| format!(">{name}{i}\n{body}\n")).collect();
            std::fs::write(&p, text).unwrap();
            p
        };
        let (a, b) = (write("a", 20), write("b", 10));
        let (rx, reader) = read_batches(&[a, b]);
        let batches: Vec<Vec<String>> = rx
            .iter()
            .map(|batch| {
                batch
                    .iter()
                    .map(|r| String::from_utf8(r.id.clone()).unwrap())
                    .collect()
            })
            .collect();
        reader.join().unwrap().unwrap();
        let names = |p: &str, r: std::ops::Range<usize>| -> Vec<String> {
            r.map(|i| format!("{p}{i}")).collect()
        };
        assert_eq!(
            batches,
            vec![
                names("a", 0..8),
                names("b", 0..8),
                names("a", 8..16),
                names("b", 8..10),
                names("a", 16..20),
            ]
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn iterative_baiting_grows_from_a_partial_seed_and_stops_when_saturated() {
        // 12 kb circular genome, seeded with its first 5 kb only; 2 kb reads every 25 bp,
        // plus unrelated nuclear reads that must never be recruited.
        let genome = pseudo_random(12_000, 7);
        let circ = [genome.clone(), genome[..2000].to_vec()].concat();
        let mut fasta = String::new();
        for (i, st) in (0..genome.len()).step_by(25).enumerate() {
            let r = std::str::from_utf8(&circ[st..st + 2000])
                .unwrap()
                .to_string();
            fasta.push_str(&format!(
                ">g{i}
{r}
"
            ));
        }
        for i in 0..200 {
            let r = pseudo_random(2000, 1000 + i);
            fasta.push_str(&format!(
                ">n{i}
{}
",
                std::str::from_utf8(&r).unwrap()
            ));
        }
        let dir = std::env::temp_dir().join(format!("ovasm-iter-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let reads = dir.join("reads.fa");
        let mut fastq = String::new();
        let lines: Vec<_> = fasta.lines().collect();
        for (i, pair) in lines.chunks(2).enumerate() {
            let qual: String = (0..pair[1].len())
                .map(|j| (33 + ((i + j) % 94) as u8) as char)
                .collect();
            fastq.push_str(&format!("@{}\n{}\n+\n{}\n", &pair[0][1..], pair[1], qual));
        }
        std::fs::write(&reads, fastq).unwrap();
        let mut bait = bait_from(&genome[..5000], 21);
        bait.targets[0].seed_length = genome.len() as u64;
        let cfg = |iterations| RecruitConfig {
            replay_dir: None,
            inputs: vec![reads.clone()],
            params: ClassifyParams {
                min_cover: 0.8,
                min_hit_frac: 0.2,
                numt_min_block: 1000,
                block_gap: 50,
            },
            iterations,
            extend_min_count: 3,
            extend_depth: 200.0,
            saturation: 0.05,
            target_depth: None,
            salt: 0,
            out_dir: dir.join(format!("out{iterations}")),
            bootstrap: false,
            extend_gate: false,
            max_read_depth: None,
            max_bases: None,
        };
        let one = run(bait.clone(), cfg(1)).unwrap();
        let mut single_cfg = cfg(1);
        single_cfg.out_dir = dir.join("replayed-one");
        single_cfg.replay_dir = Some(dir.join("cache-one"));
        let single = run(bait.clone(), single_cfg).unwrap();
        assert_eq!(
            single
                .input_replay
                .as_ref()
                .unwrap()
                .final_output_input_passes,
            0
        );
        assert_eq!(
            std::fs::read(&one.outputs[0].output).unwrap(),
            std::fs::read(&single.outputs[0].output).unwrap()
        );
        let many = run(bait.clone(), cfg(40)).unwrap();
        let mut replay_cfg = cfg(40);
        replay_cfg.out_dir = dir.join("replayed");
        replay_cfg.replay_dir = Some(dir.join("cache"));
        let cached = run(bait, replay_cfg).unwrap();
        assert_eq!(many.iterations_run, cached.iterations_run);
        for (a, b) in many.iterations.iter().zip(&cached.iterations) {
            let mut a = serde_json::to_value(a).unwrap();
            let mut b = serde_json::to_value(b).unwrap();
            a.as_object_mut().unwrap().remove("seconds");
            b.as_object_mut().unwrap().remove("seconds");
            assert_eq!(a, b);
        }
        for (a, b) in many.outputs.iter().zip(&cached.outputs) {
            assert_eq!(
                std::fs::read(&a.output).unwrap(),
                std::fs::read(&b.output).unwrap()
            );
        }
        assert_eq!(
            std::fs::read(&many.numt_output).unwrap(),
            std::fs::read(&cached.numt_output).unwrap()
        );
        assert_eq!(
            cached
                .input_replay
                .as_ref()
                .unwrap()
                .source_classification_passes,
            1
        );
        assert!(!dir.join("cache").exists());
        let genome_reads = (genome.len() / 25) as u64;
        assert!(
            one.outputs[0].recruited_reads < genome_reads / 2,
            "a partial seed alone"
        );
        assert_eq!(
            many.outputs[0].recruited_reads, genome_reads,
            "iteration reaches every read"
        );
        assert!(
            many.iterations_run < 40,
            "stops once the bait no longer grows"
        );
        let last = many.iterations.last().unwrap();
        assert!(last.bait_kmers_added.iter().all(|&a| a == 0));
        // later passes extend only from reads not used before
        let used: u64 = many
            .iterations
            .iter()
            .map(|it| it.extended_from_reads[0])
            .sum();
        assert!(
            used <= genome_reads,
            "each read extends the bait at most once: {used}"
        );
    }

    /// The same biological result, pass by pass and file by file (timings aside).
    fn assert_same_recruitment(a: &RecruitReport, b: &RecruitReport) {
        assert_eq!(a.iterations_run, b.iterations_run);
        for (x, y) in a.iterations.iter().zip(&b.iterations) {
            let mut x = serde_json::to_value(x).unwrap();
            let mut y = serde_json::to_value(y).unwrap();
            x.as_object_mut().unwrap().remove("seconds");
            y.as_object_mut().unwrap().remove("seconds");
            assert_eq!(x, y);
        }
        for (x, y) in a.outputs.iter().zip(&b.outputs) {
            assert_eq!(
                std::fs::read(&x.output).unwrap(),
                std::fs::read(&y.output).unwrap()
            );
            assert_eq!(
                (x.stop_reason, x.kept_reads),
                (y.stop_reason, y.kept_reads)
            );
        }
        assert_eq!(
            std::fs::read(&a.numt_output).unwrap(),
            std::fs::read(&b.numt_output).unwrap()
        );
        assert_eq!(a.core_depth, b.core_depth);
    }

    /// FASTA of 3 kb reads every `step` bp around the circular `genome`, named `tag<i>`.
    fn circle_reads(genome: &[u8], step: usize, tag: &str) -> String {
        let circ = [genome, &genome[..3000]].concat();
        (0..genome.len())
            .step_by(step)
            .enumerate()
            .map(|(i, st)| {
                format!(
                    ">{tag}{i}\n{}\n",
                    std::str::from_utf8(&circ[st..st + 3000]).unwrap()
                )
            })
            .collect()
    }

    /// 3 kb reads every `step` bp along a linear genome.
    fn linear_reads(genome: &[u8], step: usize, tag: &str) -> String {
        (0..genome.len() - 3000)
            .step_by(step)
            .enumerate()
            .map(|(i, st)| {
                format!(
                    ">{tag}{i}\n{}\n",
                    std::str::from_utf8(&genome[st..st + 3000]).unwrap()
                )
            })
            .collect()
    }

    fn count_reads(path: &str, tag: u8) -> usize {
        let mut n = 0;
        for_each_rec(Path::new(path), |r| n += usize::from(r.id[0] == tag)).unwrap();
        n
    }

    fn test_config(fa: &Path, out: PathBuf, iterations: usize, extend_gate: bool) -> RecruitConfig {
        RecruitConfig {
            replay_dir: None,
            inputs: vec![fa.to_path_buf()],
            params: ClassifyParams {
                min_cover: 0.8,
                min_hit_frac: 0.2,
                numt_min_block: 1000,
                block_gap: 50,
            },
            iterations,
            extend_min_count: 3,
            // no subsampling: every new read extends
            extend_depth: 1e9,
            saturation: 0.0,
            target_depth: None,
            salt: 0,
            out_dir: out,
            bootstrap: false,
            extend_gate,
            max_read_depth: None,
            max_bases: None,
        }
    }

    #[test]
    fn the_gate_keeps_extension_off_organelle_like_reads_far_shallower_than_the_target() {
        // As in `the_extension_gate_keeps_a_stray_seed_out_of_the_nucleus` (a 12 kb
        // mitochondrion at 120x seeded with 20 kb of another genome at 25x), but beside a nuclear
        // genome at 5x: the other genome is then five times the nuclear depth and passes the
        // nuclear-relative test, as everything did where the nuclear peak came out at 2x. A
        // quarter of the mitochondrion's own depth (30x) it does not reach.
        let mito = pseudo_random(12_000, 7);
        let other = pseudo_random(200_000, 41);
        let nuclear = pseudo_random(300_000, 33);
        let dir = std::env::temp_dir().join(format!("ovasm-floor-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("reads.fa");
        let text = circle_reads(&mito, 25, "m")
            + &linear_reads(&other, 120, "o")
            + &linear_reads(&nuclear, 600, "n");
        std::fs::write(&fa, text).unwrap();
        let mut bait = bait_from(&mito, 21);
        bait.add_sequence(0, &other[100_000..120_000]);
        let plain = run(bait.clone(), test_config(&fa, dir.join("plain"), 30, false)).unwrap();
        let gated = run(bait, test_config(&fa, dir.join("gated"), 30, true)).unwrap();
        // the reads of the other genome that the stray stretch covers by 80%
        let stray = (0..other.len() - 3000)
            .step_by(120)
            .filter(|&st| (99_400..=117_600).contains(&st))
            .count();
        let other_reads = |r: &RecruitReport| count_reads(&r.outputs[0].output, b'o');
        assert!(
            other_reads(&plain) > stray + 100,
            "ungated, extension walks along the other genome: {} of {stray}",
            other_reads(&plain)
        );
        assert_eq!(other_reads(&gated), stray);
        assert_eq!(count_reads(&gated.outputs[0].output, b'm'), 480);
        assert_eq!(count_reads(&gated.outputs[0].output, b'n'), 0);
        let core = gated.core_depth[0].unwrap();
        assert!((core - 119.0).abs() < 3.0, "core depth {core}");
        let first = &gated.iterations[0];
        assert_eq!(first.extension_gated_out[0], 0, "organelle-like by depth");
        assert_eq!(first.extension_shallow_out[0], stray as u64);
        assert!(gated.outputs[0].converged);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn the_gate_keeps_extension_off_a_high_copy_repeat_far_deeper_than_the_target() {
        // A 12 kb mitochondrion at 120x seeded with 3 kb of a 6 kb repeat at 1,500x, beside a
        // nuclear genome at 5x. The first pass recruits more repeat reads than mitochondrial
        // ones, yet the target's depth is the mitochondrion's (`target_depth`), and the repeat,
        // over ten times deeper, does not extend it.
        let mito = pseudo_random(12_000, 7);
        let repeat = pseudo_random(6_000, 43);
        let nuclear = pseudo_random(300_000, 33);
        let dir = std::env::temp_dir().join(format!("ovasm-ceiling-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("reads.fa");
        let text = circle_reads(&mito, 25, "m")
            + &linear_reads(&repeat, 2, "r")
            + &linear_reads(&nuclear, 600, "n");
        std::fs::write(&fa, text).unwrap();
        let mut bait = bait_from(&mito, 21);
        bait.add_sequence(0, &repeat[1500..4500]);
        let plain = run(bait.clone(), test_config(&fa, dir.join("plain"), 30, false)).unwrap();
        let gated = run(bait, test_config(&fa, dir.join("gated"), 30, true)).unwrap();
        let repeat_reads = |r: &RecruitReport| count_reads(&r.outputs[0].output, b'r');
        // the repeat's reads that the seeded stretch covers by 80%: starts 900 to 2100
        let seeded = (0..3000).step_by(2).filter(|&st| (900..=2100).contains(&st)).count();
        assert_eq!(repeat_reads(&plain), 1500, "ungated, extension walks the whole repeat");
        assert_eq!(repeat_reads(&gated), seeded);
        assert_eq!(count_reads(&gated.outputs[0].output, b'm'), 480);
        let core = gated.core_depth[0].unwrap();
        assert!((core - 119.0).abs() < 5.0, "core depth {core}");
        assert!(gated.iterations[0].extension_deep_out[0] > 0);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_gated_target_growing_steadily_is_not_saturated() {
        // A 24 kb mitochondrion at 120x seeded with 4 kb of itself grows by the same 32 reads
        // every pass, a smaller share each time: under a 20% saturation the ungated target stops
        // after a few passes; the gated one keeps going until it closes the circle. (The first
        // pass's reads count as its growth, so a close seed still stops at its second pass.)
        let mito = pseudo_random(24_000, 9);
        let nuclear = pseudo_random(300_000, 33);
        let dir = std::env::temp_dir().join(format!("ovasm-steady-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("reads.fa");
        std::fs::write(&fa, circle_reads(&mito, 25, "m") + &linear_reads(&nuclear, 600, "n")).unwrap();
        let mut bait = BaitSet::new(21);
        let t = bait.target("mitochondrion").unwrap();
        bait.add_sequence(t, &mito[..4000]);
        let config = |name: &str, gate: bool| {
            let mut c = test_config(&fa, dir.join(name), 40, gate);
            c.saturation = 0.2;
            c
        };
        let plain = run(bait.clone(), config("plain", false)).unwrap();
        let gated = run(bait, config("gated", true)).unwrap();
        let mito_reads = |r: &RecruitReport| count_reads(&r.outputs[0].output, b'm');
        assert_eq!(plain.outputs[0].stop_reason, "saturated");
        assert!(mito_reads(&plain) < 480, "{}", mito_reads(&plain));
        assert_eq!(mito_reads(&gated), 960);
        assert!(gated.outputs[0].converged);
        assert_eq!(count_reads(&gated.outputs[0].output, b'n'), 0);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn a_high_copy_repeat_does_not_set_a_targets_depth() {
        // the first pass of the Lophophora williamsii plastid: most reads a repeat at ~700,000x,
        // the plastid's own at ~440x, some nuclear reads at 22x below the organelle threshold
        let mut reads: Vec<(f64, usize)> = (0..1800)
            .map(|i| (650_000.0 + 100.0 * i as f64, 10_000))
            .collect();
        reads.extend((0..900).map(|i| (400.0 + (i % 80) as f64, 10_000)));
        reads.extend((0..100).map(|_| (22.0, 10_000)));
        let d = target_depth(&mut reads, 66.0);
        assert!((400.0..480.0).contains(&d), "{d}");
        // a mitochondrion with plastid-deep reads among its own keeps its median, as Medicago
        // polymorpha's did (3,212x, its sequence at 1,045x)
        let mut mito: Vec<(f64, usize)> = (0..2000)
            .map(|i| (900.0 + (i % 200) as f64, 10_000))
            .collect();
        mito.extend((0..2500).map(|i| (2_000.0 + 11.2 * i as f64, 10_000)));
        let mut sorted = mito.clone();
        sorted.sort_by(|a, b| a.0.total_cmp(&b.0));
        let median = sorted[sorted.len() / 2].0;
        assert!(median > 4.0 * 1100.0, "{median}");
        assert_eq!(target_depth(&mut mito, 165.0), median);
        // no organelle-deep read: the median
        assert_eq!(
            target_depth(&mut [(20.0, 5000), (22.0, 5000), (25.0, 5000)], 66.0),
            22.0
        );
    }

    #[test]
    fn outputs_are_downsampled_to_the_read_depth_and_to_the_base_cap() {
        // a 20 kb plastid at 600x (4,000 reads of 3 kb, 12 Mb) beside a nuclear genome at 10x
        let plastid = pseudo_random(20_000, 21);
        let nuclear = pseudo_random(200_000, 33);
        let dir = std::env::temp_dir().join(format!("ovasm-cap-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("reads.fa");
        std::fs::write(
            &fa,
            circle_reads(&plastid, 5, "p") + &linear_reads(&nuclear, 300, "n"),
        )
        .unwrap();
        let mut bait = BaitSet::new(21);
        let t = bait.target("plastid").unwrap();
        bait.add_sequence(t, &plastid);
        bait.targets[t].seed_length = plastid.len() as u64;
        let mut by_depth = test_config(&fa, dir.join("depth"), 1, false);
        by_depth.max_read_depth = Some(150.0);
        let d = run(bait.clone(), by_depth).unwrap();
        let o = &d.outputs[t];
        assert_eq!(o.recruited_reads, 4000);
        let depth = o.read_depth.unwrap();
        assert!((depth - 596.0).abs() < 5.0, "read depth {depth}");
        assert!((o.implied_size.unwrap() as f64 - 20_000.0).abs() < 500.0);
        assert_eq!(o.sampled_by, Some("max_read_depth"));
        assert!(
            (o.kept_reads as f64 - 4000.0 * 150.0 / depth).abs() < 100.0,
            "{}",
            o.kept_reads
        );
        assert_eq!(o.stop_reason, "single_pass");

        let mut by_bases = test_config(&fa, dir.join("bases"), 1, false);
        by_bases.max_read_depth = Some(1000.0);
        by_bases.max_bases = Some(3_000_000);
        let b = run(bait, by_bases).unwrap();
        let o = &b.outputs[t];
        assert_eq!(o.sampled_by, Some("max_bases"));
        assert!((o.kept_bases as f64 - 3e6).abs() < 2e5, "{}", o.kept_bases);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn only_plastids_and_mitochondria_have_a_size_cap() {
        assert_eq!(over_size_cap("plastid", Some(1_000_000)), None);
        assert_eq!(
            over_size_cap("plastid", Some(1_000_001)),
            Some((1_000_001, 1_000_000))
        );
        assert_eq!(over_size_cap("mitochondrion", Some(19_000_000)), None);
        assert!(over_size_cap("mitochondrion", Some(21_000_000)).is_some());
        assert_eq!(over_size_cap("mito", Some(u64::MAX)), None);
        assert_eq!(over_size_cap("plastid", None), None);
    }

    #[test]
    fn a_depth_sample_is_the_same_whatever_the_order() {
        let reads: Vec<(u64, Vec<u8>)> = (0..10_000u64)
            .map(|i| (name_key(format!("r{i}").as_bytes()), pseudo_random(50, i)))
            .collect();
        let mut forward = Sample::default();
        let mut backward = Sample::default();
        for (k, s) in &reads {
            forward.offer(*k, s);
        }
        for (k, s) in reads.iter().rev() {
            backward.offer(*k, s);
        }
        assert_eq!(forward.reads(), backward.reads());
        let mut smallest: Vec<&(u64, Vec<u8>)> = reads.iter().collect();
        smallest.sort_unstable_by_key(|r| r.0);
        let expected: Vec<&[u8]> = smallest[..DEPTH_SAMPLE]
            .iter()
            .map(|r| r.1.as_slice())
            .collect();
        assert_eq!(forward.reads(), expected);
    }
}
