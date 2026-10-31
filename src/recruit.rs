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
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
/// Profiling build: nanoseconds summed over threads.
pub static PROF_DELTA_NS: AtomicU64 = AtomicU64::new(0);
pub static PROF_CLASSIFY_NS: AtomicU64 = AtomicU64::new(0);
pub static PROF_READER_WORK_NS: AtomicU64 = AtomicU64::new(0);
pub static PROF_READER_SEND_NS: AtomicU64 = AtomicU64::new(0);
pub static PROF_PARSE_NS: AtomicU64 = AtomicU64::new(0);
fn prof_take(a: &AtomicU64) -> f64 { a.swap(0, Relaxed) as f64 / 1e9 }
/// `OVASM_PROFILE` set: each recruitment pass reports where its time went.
pub fn profiling() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("OVASM_PROFILE").is_some())
}
use crate::kmer::CanonicalKmers;
use crate::seedfree::{DepthGate, ReadDepth};

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
    /// A read the bait covers whose own k-mers are mostly repeats of one another (a smaller share
    /// distinct) is a tandem array, not organelle sequence, and is not recruited
    /// (`HIFI_MIN_DISTINCT`); 0 recruits it.
    pub min_distinct: f64,
}

/// `ClassifyParams::min_distinct` for HiFi reads. Once a satellite's unit is in a target's bait,
/// every read of its arrays is covered: Medicago truncatula A17's mitochondrion recruited 115,110
/// reads in its second pass, nine tenths of them a nuclear tandem repeat at about its own depth
/// (21-mers of a read 6-16% distinct; organelle reads near 100%), and assembling them took 365 GB.
/// A read across two tandem copies of an organelle segment is still about two thirds distinct.
pub const HIFI_MIN_DISTINCT: f64 = 0.5;

/// Share of a read's k-mers that are distinct, over an eighth of them chosen by hash (a repeated
/// k-mer is chosen every time); 1 for a read too short to judge.
pub fn distinct_share(seq: &[u8], k: usize) -> f64 {
    let mut seen = FxHashSet::default();
    let mut n = 0usize;
    for (_, kmer) in CanonicalKmers::new(seq, k) {
        if kmer.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 61 != 0 {
            continue;
        }
        n += 1;
        seen.insert(kmer);
    }
    if n < 20 {
        1.0
    } else {
        seen.len() as f64 / n as f64
    }
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
    /// Covered, but a tandem array (`ClassifyParams::min_distinct`): called `Other`.
    pub tandem: bool,
    /// The best target's gap-bridged coverage in bases, the most of any target (`still_other`).
    pub cover: usize,
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
    bait.scan_masks(seq, |pos, mask| {
        h.total += 1;
        if mask == 0 {
            return;
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
    });
    for t in 0..n {
        if h.hits[t] > 0 {
            let len = last_end[t] - block_start[t];
            h.block_cover[t] += len;
            h.longest[t] = h.longest[t].max(len);
        }
    }
    h
}

/// What a bait change can add to a read's coverage, for `still_other`: `None` when no k-mer of the
/// read changed, else an upper bound on how much any target's gap-bridged coverage
/// (`TargetHits::block_cover`) can grow. A target's new hits are its old ones plus some of the
/// changed k-mers. Coverage only grows as hits are added, and the changed k-mers' own blocks
/// (bridged as `target_hits` bridges) add at most their length plus one bridged gap of at most
/// `block_gap` on each side to a target's blocks: hence their length plus `2 * block_gap` a block.
fn added_cover(seq: &[u8], delta: &BaitSet, block_gap: usize) -> Option<usize> {
    let k = delta.k;
    let (mut any, mut start, mut last_end, mut cover, mut blocks) = (false, 0, 0, 0, 0);
    delta.scan_masks(seq, |pos, mask| {
        if mask == 0 {
            return;
        }
        if !any {
            any = true;
            start = pos;
        } else if pos > last_end + block_gap {
            cover += last_end - start;
            blocks += 1;
            start = pos;
        }
        last_end = last_end.max(pos + k);
    });
    any.then(|| cover + (last_end - start) + 2 * block_gap * (blocks + 1))
}

/// A read last called `Other` with at most `bound` coverage by any target stays `Other` when even
/// `bound + added` reaches neither the organelle call's `min_cover` nor a NUMT candidate's block
/// (a block is never longer than the coverage): then the new bound, and the read is not
/// reclassified. `u16::MAX` stands for an unknown bound (coverage past 65,535 bp).
///
/// Bait extension reclassified reads holding any changed k-mer, and the k-mers of a repeat
/// reach most reads of a library: on Eleocharis quinqueflora HiFi two of sixteen replayed passes
/// reclassified 2.25 M of 2.25 M reads, at 160-175 s each, for a few organelle reads.
fn still_other(call: Call, bound: u16, added: usize, len: usize, p: &ClassifyParams) -> Option<u16> {
    if call != Call::Other || bound == u16::MAX {
        return None;
    }
    let new = bound as usize + added;
    let below = new < p.numt_min_block && (new as f64) < p.min_cover * len as f64;
    (below && new < u16::MAX as usize).then_some(new as u16)
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
    let mut score = ReadScore {
        cover_frac: block_cover[best] as f64 / seq.len() as f64,
        hit_frac: hits[best] as f64 / total as f64,
        longest_block: longest[best],
        tandem: false,
        cover: block_cover[best],
    };
    let call = if score.cover_frac >= p.min_cover && score.hit_frac >= p.min_hit_frac {
        if p.min_distinct > 0.0 && distinct_share(seq, bait.k) < p.min_distinct {
            score.tandem = true;
            Call::Other
        } else {
            Call::Organelle(tied)
        }
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
    /// Partial reads at each target's edges that its grown bait came to cover by `min_cover`
    /// within the pass, and extended it again (`FRONTIER_STEPS`); the next pass recruits them.
    pub frontier_reads: Vec<u64>,
    /// Steps of that within-pass extension.
    pub frontier_steps: Vec<u64>,
    /// Reads of this pass's extension (and its frontier) the gate turned down as a whole but took
    /// by the sequence they add, the k-mers the bait lacked (`extends`).
    pub extension_new_part: Vec<u64>,
    /// Reads the bait covered but that are tandem arrays (`ClassifyParams::min_distinct`), not
    /// recruited; among the reads this pass classified (a replayed pass reclassifies only those
    /// the bait's changes touch).
    pub tandem_reads: u64,
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
    /// Each extending target's NUMT candidates that this pass classified anew: the partial reads
    /// at the edges its bait last grew to (`FRONTIER_STEPS`).
    frontier: Vec<Sample<FRONTIER_READS>>,
    numt: TargetCounts,
    reads: u64,
    bases: u64,
    multi_target: u64,
    reclassified: u64,
    /// Reads classified in this pass that the bait covered but are tandem arrays.
    tandem: u64,
}

/// Reads per depth estimate. A median of 2,000 read depths is within a few percent, and the
/// sample costs one pass's memory for about 2,000 reads per target.
const DEPTH_SAMPLE: usize = 2000;

/// Most frontier reads a pass keeps per target (`ScanResult::frontier`): about 150 Mb of HiFi
/// reads, a few hundred at each of dozens of edges; a deep target's are a uniform subsample.
const FRONTIER_READS: usize = 10_000;

/// The `N` reads with the smallest name hashes among those offered: the same reads whatever the
/// input order, in memory for 2 x `N` reads at most.
#[derive(Default)]
struct Sample<const N: usize = DEPTH_SAMPLE> {
    kept: Vec<(u64, Vec<u8>)>,
    /// Hashes from here on cannot be among the smallest any more.
    bound: Option<u64>,
}

impl<const N: usize> Sample<N> {
    fn offer(&mut self, key: u64, seq: &[u8]) {
        if self.bound.is_some_and(|b| key >= b) {
            return;
        }
        self.kept.push((key, seq.to_vec()));
        if self.kept.len() == 2 * N {
            self.kept.select_nth_unstable_by_key(N, |&(k, _)| k);
            self.bound = Some(self.kept[N].0);
            self.kept.truncate(N);
        }
    }

    /// The sampled reads, smallest hash first.
    fn reads(&self) -> Vec<&[u8]> {
        let mut kept: Vec<&(u64, Vec<u8>)> = self.kept.iter().collect();
        kept.sort_unstable_by_key(|&&(k, _)| k);
        kept.truncate(N);
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

/// Whether one measure of a target has settled over a pass: it grew by less than `saturation`
/// of what it had and, with the gate, by less than `STEADY_GROWTH` of its growth the pass before.
fn settled(now: u64, previous: u64, last_growth: u64, saturation: f64, gated: bool) -> bool {
    let grown = now.saturating_sub(previous);
    let steady = gated && grown > 0 && grown as f64 >= STEADY_GROWTH * last_growth as f64;
    previous > 0 && (grown as f64) < saturation * previous as f64 && !steady
}

/// Whether the genome a gated target's reads stand for has settled over a pass: `settled`, and
/// grown by less than a read length (`read_len`, its reads' mean). Growth falls by half or more
/// from one pass to the next when edges meet or reach the end of a linear piece, and a genome
/// still growing by more than a read is still extending at its other edges: the Lophophora
/// williamsii mitochondrion gained 1,917 reads in a pass and 833 (3.6%) in the next, under half,
/// while its genome grew by 65 kb (eight of its 8.1 kb reads), and recruitment stopped there with
/// 358 kb of the whole-library graph never recruited. Every other target of the benchmark pilots
/// settled having grown by at most 7.4 kb.
fn genome_settled(
    now: u64,
    previous: u64,
    last_growth: u64,
    saturation: f64,
    read_len: u64,
) -> bool {
    let grown = now.saturating_sub(previous);
    settled(now, previous, last_growth, saturation, true) && grown < read_len
}

/// Fewest reads a target's extension is subsampled to. A distant seed recruits only a few of the
/// reads at each place in the first passes, so its reads' depth over the genome they stand for
/// overstates how deep they cover any place: the Schoenoplectus lacustris plastid (5,355x, seed
/// containment 0.061) recruited 289 reads in its first pass, extension took 9 of them and 4-34 a
/// pass after, and the plastid stopped at 9 kb (0.2.1, extending from every read: 169 kb).
const MIN_EXTEND_READS: f64 = 5000.0;

/// Share of a target's reads its extension takes: down to `extend_depth` over the genome they
/// stand for, but not below `MIN_EXTEND_READS` reads.
fn extension_fraction(reads: u64, bases: u64, genome: f64, extend_depth: f64) -> f64 {
    let depth = bases as f64 / genome;
    if depth <= extend_depth {
        return 1.0;
    }
    (extend_depth / depth)
        .max(MIN_EXTEND_READS / reads.max(1) as f64)
        .min(1.0)
}

/// Most steps of within-pass extension from a target's frontier reads (`ScanResult::frontier`).
/// A read joins a target only when its bait covers 80% of it, so one pass moves an edge by a tenth
/// of a read length or so and a large mitochondrion takes a pass, a full scan of the input, for
/// every step: Lophophora williamsii's (2.2 Mb) stood at 1.1 Mb after 60 passes. The pass's
/// partial reads at the edges reach most of a read length past them; reclassified against the
/// grown bait, those it covers by 80% are what the next pass would recruit, and the bait grows
/// from them again under the same rules, until no more qualify.
const FRONTIER_STEPS: usize = 10;

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
    let mut tp = Instant::now();
    while let Some(record) = parser.next() {
        let record = record.with_context(|| format!("malformed record in {}", path.display()))?;
        bases += record.num_bases();
        batch.push(Rec {
            id: record.id().to_vec(),
            seq: record.seq().into_owned(),
            qual: record.qual().map(<[u8]>::to_vec),
        });
        if bases >= BATCH_BASES {
            PROF_PARSE_NS.fetch_add(tp.elapsed().as_nanos() as u64, Relaxed);
            let ts = Instant::now();
            if tx.send(std::mem::take(&mut batch)).is_err() {
                return Ok(());
            }
            PROF_READER_SEND_NS.fetch_add(ts.elapsed().as_nanos() as u64, Relaxed);
            tp = Instant::now();
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
    gate_sketch: Option<&mut crate::seedfree::GateSketch>,
) -> Result<ScanResult> {
    let (out_dir, numt_path) = (w.out_dir, w.numt_path);
    let replaying = replay.as_ref().is_some_and(|r| r.ready);
    // the depth gate's sketch, counted on its own thread from the k-mers the workers sample
    let sampling = gate_sketch.as_ref().map(|s| s.sampling());
    let sketcher = gate_sketch.map(|s| {
        let (tx, rx) = crossbeam_channel::bounded::<(Vec<Vec<u64>>, u64)>(4);
        let mut owned = std::mem::take(s);
        let thread = std::thread::spawn(move || {
            for (reads, bases) in rx {
                owned.add(reads, bases);
            }
            owned
        });
        (s, tx, thread)
    });
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
        frontier: (0..n_targets).map(|_| Sample::default()).collect(),
        numt: TargetCounts::default(),
        reads: 0,
        bases: 0,
        multi_target: 0,
        reclassified: 0,
        tandem: 0,
    };
    let (mut t_wait, mut t_par, mut t_serial, mut t_check, mut n_batches) = (0f64, 0f64, 0f64, 0f64, 0u64);
    let scan_start = Instant::now();
    let _ = (prof_take(&PROF_DELTA_NS), prof_take(&PROF_CLASSIFY_NS), prof_take(&PROF_READER_WORK_NS), prof_take(&PROF_READER_SEND_NS), prof_take(&PROF_PARSE_NS));
    let mut rx_iter = rx.into_iter();
    loop {
        let tw = Instant::now();
        let Some(batch) = rx_iter.next() else { break };
        t_wait += tw.elapsed().as_secs_f64();
        n_batches += 1;
        let tc = Instant::now();
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
        t_check += tc.elapsed().as_secs_f64();
        let offset = out.reads as usize;
        let previous = if replaying {
            let cache = replay.as_ref().unwrap();
            let range = offset..offset + batch.len();
            Some((
                cache
                    .calls
                    .get(range.clone())
                    .context("input has more records than the recruitment cache")?,
                &cache.cover_bounds[range],
            ))
        } else {
            None
        };
        out.reads += batch.len() as u64;
        let batch_bases = batch.iter().map(|r| r.rec.seq.len() as u64).sum::<u64>();
        out.bases += batch_bases;
        let tpar = Instant::now();
        if let (Some((k, threshold)), Some((_, tx, _))) = (sampling, sketcher.as_ref()) {
            let reads = batch
                .par_iter()
                .map(|r| crate::seedfree::GateSketch::sample(k, threshold, &r.rec.seq))
                .collect();
            tx.send((reads, batch_bases)).context("depth sketch thread stopped")?;
        }
        // the first pass encodes the batch for the cache while it classifies it
        let storing = !replaying && replay.is_some();
        let (encoded, calls) = rayon::join(
            || storing.then(|| crate::replay::Replay::encode(&batch)),
            || -> Vec<(Call, bool, bool, u16)> {
                batch
                    .par_iter()
                    .enumerate()
                    .map(|(i, rec)| {
                        let seq = &rec.rec.seq;
                        let classified = || {
                            let tcl = Instant::now();
                            let (call, score) = classify(seq, bait, p);
                            PROF_CLASSIFY_NS.fetch_add(tcl.elapsed().as_nanos() as u64, Relaxed);
                            let bound = u16::try_from(score.cover).unwrap_or(u16::MAX);
                            (call, true, score.tandem, bound)
                        };
                        let Some((calls, bounds)) = previous else {
                            return classified();
                        };
                        let (call, bound) = (calls[i], bounds[i]);
                        if final_outputs {
                            return (call, false, false, bound);
                        }
                        let delta = delta.expect("replay has change set");
                        let td = Instant::now();
                        let changed = CanonicalKmers::new(seq, bait.k).any(|(_, key)| delta.get(key) != 0);
                        // unchanged: as it was; changed: kept with a new bound only when it
                        // cannot leave `Other` (the bound's full look, for such reads alone)
                        let kept = if !changed {
                            Some(bound)
                        } else if call == Call::Other && bound != u16::MAX {
                            added_cover(seq, delta, p.block_gap)
                                .and_then(|added| still_other(call, bound, added, seq.len(), p))
                        } else {
                            None
                        };
                        PROF_DELTA_NS.fetch_add(td.elapsed().as_nanos() as u64, Relaxed);
                        match kept {
                            Some(bound) => (call, false, false, bound),
                            None => classified(),
                        }
                    })
                    .collect()
            },
        );
        if let Some(encoded) = encoded {
            replay.as_mut().expect("storing into the cache").write(encoded)?;
        }
        t_par += tpar.elapsed().as_secs_f64();
        let tser = Instant::now();
        for (i, ((call, changed, tandem, bound), entry)) in calls.into_iter().zip(batch).enumerate() {
            out.reclassified += u64::from(changed);
            out.tandem += u64::from(tandem);
            if let Some(cache) = replay.as_mut() {
                if replaying {
                    cache.calls[offset + i] = call;
                    cache.cover_bounds[offset + i] = bound;
                } else {
                    cache.calls.push(call);
                    cache.cover_bounds.push(bound);
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
                Call::NumtCandidate(t) => {
                    if w.keep.is_none() {
                        write_rec(&mut numt_w, &rec)?;
                    }
                    out.numt.reads += 1;
                    out.numt.bases += len;
                    if changed && w.keep.map_or(true, |k| k[t] != Keep::Nothing) {
                        out.frontier[t].offer(name_key(&rec.id), &rec.seq);
                    }
                }
                Call::Other => {}
            }
        }
        t_serial += tser.elapsed().as_secs_f64();
    }
    let t_loop = scan_start.elapsed().as_secs_f64();
    reader.join().expect("reader thread panicked")?;
    if let Some((expected, reader)) = expected_input {
        anyhow::ensure!(
            expected.recv().is_err(),
            "original input has fewer batches than its recruitment cache"
        );
        reader.join().expect("cache reader panicked")?;
    }
    if let Some((s, tx, thread)) = sketcher {
        drop(tx);
        *s = thread.join().expect("depth sketch thread panicked");
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
    if profiling() { eprintln!(
        "[ovasm] profile scan kind={} batches={} reads={} reclassified={} wall={:.1}s wait={:.1}s check={:.1}s par={:.1}s serial={:.1}s | thread-summed delta={:.1}s classify={:.1}s | reader work={:.1}s send_blocked={:.1}s gz_parse={:.1}s",
        if final_outputs { "final" } else if replaying { "replay" } else { "source" },
        n_batches, out.reads, out.reclassified, scan_start.elapsed().as_secs_f64(), t_wait, t_check, t_par, t_serial,
        prof_take(&PROF_DELTA_NS), prof_take(&PROF_CLASSIFY_NS),
        prof_take(&PROF_READER_WORK_NS), prof_take(&PROF_READER_SEND_NS), prof_take(&PROF_PARSE_NS)
    ); }
    let _ = t_loop;
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

/// What a read is downsampled by: its whole name, or with `mates_together` the name's first word.
fn sample_name(id: &[u8], mates_together: bool) -> &[u8] {
    if mates_together {
        id.split(|b| b.is_ascii_whitespace()).next().unwrap_or(id)
    } else {
        id
    }
}

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
    /// (`find discover`'s test). Iterative extension follows any read the bait covers by
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
    /// Downsample by the first word of each read's name, which a pair's mates share, so mates
    /// are kept or left out together (`sample_name`). Reconstruction closes circles with read
    /// pairs: sampled by whole names, ZM4's plastid at 11.6% kept both mates of 1.3% of its pairs,
    /// and of 21 pairs across its join one was left.
    pub mates_together: bool,
    /// Short reads: count the library's k-mers once (count-min, after the first pass) and let a
    /// read extend a target only when the k-mers it would add are at the target's depth
    /// (`kmer_extends`). Iterative extension follows any read the bait covers by `min_cover`,
    /// and from a NUMT's edge into the nucleus: WL354's mitochondrion grew from 0.72 M to 54 M
    /// reads in 20 passes, CRR1151719's to 104 M, and assembly ran out of memory or time.
    pub kmer_gate: bool,
}

/// Counters of the k-mer gate's table: 2^31 of 16 bits, 4 GiB.
const KMER_GATE_SLOTS_LOG2: u32 = 31;
/// A read extends a target when the median count of the k-mers it would add is from this share
/// of the target's own k-mer depth up to `KMER_GATE_CEILING` times it: the nucleus is far
/// shallower than an organelle, a plastid far deeper than a mitochondrion.
const KMER_GATE_FLOOR: f64 = 0.25;
const KMER_GATE_CEILING: f64 = 10.0;
/// Fewer new k-mers than this: the read adds nothing to judge, and extends.
const KMER_GATE_MIN_NEW: usize = 5;

/// Whether a read may extend target `t` (`RecruitConfig::kmer_gate`): the k-mers it would add,
/// those the target's bait lacks, at the target's depth `depth` by the table's counts.
fn kmer_extends(table: &crate::frontier::KmerCounts, bait: &BaitSet, t: usize, seq: &[u8], depth: f64) -> bool {
    let bit = 1u8 << t;
    let mut new: Vec<u32> = CanonicalKmers::new(seq, bait.k)
        .filter(|&(_, km)| bait.get(km) & bit == 0)
        .map(|(_, km)| table.get(km))
        .collect();
    if new.len() < KMER_GATE_MIN_NEW {
        return true;
    }
    new.sort_unstable();
    let m = f64::from(new[new.len() / 2]);
    m >= KMER_GATE_FLOOR * depth && m <= KMER_GATE_CEILING * depth
}

/// A target's k-mer depth: the median, over its sampled reads, of the median count of each
/// read's k-mers that are in the target's bait.
fn kmer_depth(table: &crate::frontier::KmerCounts, bait: &BaitSet, t: usize, reads: &[&[u8]]) -> Option<f64> {
    let bit = 1u8 << t;
    let mut per_read: Vec<u32> = reads
        .par_iter()
        .filter_map(|s| {
            let mut v: Vec<u32> = CanonicalKmers::new(s, bait.k)
                .filter(|&(_, km)| bait.get(km) & bit != 0)
                .map(|(_, km)| table.get(km))
                .collect();
            if v.is_empty() {
                return None;
            }
            v.sort_unstable();
            Some(v[v.len() / 2])
        })
        .collect();
    if per_read.is_empty() {
        return None;
    }
    per_read.sort_unstable();
    Some(f64::from(per_read[per_read.len() / 2]))
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

/// How the extension gate takes a read, by its depth `profile`.
#[derive(Clone, Copy, PartialEq)]
enum Gauge {
    Extends,
    /// Not organelle-like (`DepthGate::profile`).
    Gated,
    /// Shallower than `EXTEND_CORE_FLOOR` of the target's core depth.
    Shallow,
    /// Deeper than `HIGH_COPY_FOLD` times it.
    Deep,
}

fn gauge(d: Option<ReadDepth>, floor: Option<f64>, ceiling: Option<f64>) -> Gauge {
    match d {
        Some(d) if d.organelle_like => {
            if floor.is_some_and(|f| d.median < f) {
                Gauge::Shallow
            } else if ceiling.is_some_and(|c| d.median > c) {
                Gauge::Deep
            } else {
                Gauge::Extends
            }
        }
        _ => Gauge::Gated,
    }
}

/// Depth from which a sampled k-mer of a read counts as deep for extending a target with core
/// depth `core`: the gate's organelle threshold, three times the nuclear depth, or halfway
/// between the nuclear depth and the core depth (log scale) when that is lower, for a target
/// less than nine times as deep as the nucleus. A nuclear k-mer stays far below either, while a
/// target's own sequence a little shallower than the rest of it keeps clear of the threshold:
/// the Lophophora williamsii mitochondrion (core 104x, nucleus 22x, threshold 66x) has 236 kb
/// at its depth between two stretches at 66-90x; reads there kept a tenth of their k-mers below
/// 66x, the gate turned them down, and recruitment stopped on both sides of the 236 kb.
fn extension_deep(g: &DepthGate, core: Option<f64>) -> f64 {
    let organelle = g.organelle_min_depth;
    core.map_or(organelle, |c| {
        organelle.min((f64::from(g.nuclear_peak_depth) * c).sqrt())
    })
}

/// The extension gate on a read (deep k-mers from `deep`, `extension_deep`): on the whole of it,
/// and when that fails, on the sequence it would add to target `t`, its k-mers the bait lacks.
/// The bait holds the rest of the read already, and it may be the organelle's own repeat: the
/// Schoenoplectus lacustris plastid's inverted repeat holds tandem arrays, 30-40 times as deep
/// as the rest of it, between its outer part and the part that leads to the small single-copy
/// region; reads across them repeat k-mers and vary in depth along their length, the gate turned
/// down most of them (and 30% of the reads its second pass drew for extension), and the plastid
/// stopped there, a 165 kb molecule without its 10.8 kb SSC. Returns the verdict and whether it
/// was the added sequence that passed.
fn extends(
    g: &DepthGate,
    bait: &BaitSet,
    t: usize,
    seq: &[u8],
    deep: f64,
    (floor, ceiling): (Option<f64>, Option<f64>),
) -> (Gauge, bool) {
    let whole = gauge(g.profile_of(seq, deep, |_| true), floor, ceiling);
    if whole == Gauge::Extends || g.k() != bait.k {
        return (whole, false);
    }
    let bit = 1u8 << t;
    let new = gauge(g.profile_of(seq, deep, |km| bait.get(km) & bit == 0), floor, ceiling);
    if new == Gauge::Extends {
        (Gauge::Extends, true)
    } else {
        (whole, false)
    }
}

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
/// frozen holdout samples, `find discover`'s organelle-like test passed 979 of 989 mitochondrial
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
    sketch: Option<crate::seedfree::GateSketch>,
    warnings: &mut Vec<String>,
) -> Result<DepthGate> {
    let gate = match (replay.filter(|r| r.ready), sketch) {
        (Some(cache), Some(s)) => {
            DepthGate::long_reads_from_sketch(s, || cache.sequences(), "first pass and recruitment cache")?
        }
        (None, Some(s)) => {
            DepthGate::long_reads_from_sketch(s, || read_batches(inputs), "first pass and inputs")?
        }
        (Some(cache), None) => DepthGate::long_reads_from(|| cache.sequences(), "recruitment cache")?,
        (None, None) => DepthGate::long_reads(inputs)?,
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
    // The same for the genome each target's reads stand for (with the depth sketch).
    let mut previous_size = vec![0u64; n_targets];
    let mut last_size_growth = vec![0u64; n_targets];
    // Depth sketch, built after the first pass when the gate or the downsampling needs it, or
    // the first time a target needs bootstrapping.
    let mut gate: Option<DepthGate> = None;
    let wants_gate = cfg.extend_gate || cfg.max_read_depth.is_some();
    // The gate's first sketch pass, counted during the first pass over the inputs instead of in
    // a pass of its own (129-181 s for the two passes on 13-15 GB of DToL HiFi): whenever a gate
    // may be built, by the condition after the first pass or by a bootstrap.
    let mut gate_sketch = ((wants_gate || cfg.bootstrap)
        && (cfg.iterations.max(1) > 1 || cfg.max_read_depth.is_some()))
    .then(crate::seedfree::GateSketch::default);
    let mut core_depth: Vec<Option<f64>> = vec![None; n_targets];
    let mut warnings: Vec<String> = Vec::new();
    let mut boot_used: Vec<FxHashSet<u64>> = vec![FxHashSet::default(); n_targets];
    // the k-mer gate's table and each target's k-mer depth, from after the first pass
    let mut kmer_table: Option<crate::frontier::KmerCounts> = None;
    let mut kmer_target_depth: Vec<Option<f64>> = vec![None; n_targets];
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
            // only the first pass: without a cache every pass reads the inputs
            if it == 0 { gate_sketch.as_mut() } else { None },
        )?;
        if crate::recruit::profiling() {
            eprintln!("[ovasm] profile iteration {} pre-scan+scan {:.1}s", it + 1, t0.elapsed().as_secs_f64());
        }
        let last = it + 1 == cfg.iterations.max(1);
        if gate.is_none() && wants_gate && (!last || cfg.max_read_depth.is_some()) {
            gate = Some(build_gate(&cfg.inputs, replay.as_ref(), gate_sketch.take(), &mut warnings)?);
        }
        let mut depths = pass_depths(gate.as_ref(), &res);
        if cfg.kmer_gate && kmer_table.is_none() && !last {
            let source = match replay.as_ref().filter(|r| r.ready) {
                Some(cache) => cache.sequences(),
                None => read_batches(&cfg.inputs),
            };
            let table = crate::frontier::KmerCounts::from_source(source, KMER_GATE_SLOTS_LOG2, bait.k)?;
            for t in 0..n_targets {
                let reads = res.all_sample[t].reads();
                kmer_target_depth[t] = kmer_depth(&table, &bait, t, &reads);
            }
            eprintln!(
                "[ovasm] k-mer gate: {} k-mers of {} reads counted in {:.0}s; target depths {:?}",
                bait.k,
                table.reads,
                table.seconds,
                kmer_target_depth
            );
            kmer_table = Some(table);
        }
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
        let mut frontier_reads = vec![0u64; n_targets];
        let mut frontier_steps = vec![0u64; n_targets];
        let mut new_part = vec![0u64; n_targets];
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
                let reads_settled =
                    settled(now, previous[t], last_growth[t], cfg.saturation, cfg.extend_gate);
                // With the depth sketch, the genome the reads stand for must have settled too: a
                // target that took in many deep reads grows by a small share of its reads while
                // its genome still grows (Medicago polymorpha's mitochondrion: 31,212 reads at
                // ~6,500x in its second pass; at its sixth, +1,755 reads, 3.6%, its genome 131 to
                // 169 kb, and 61 kb of it never recruited), and by less than a read length
                // (`genome_settled`).
                let size = depths.implied_size[t].filter(|_| cfg.extend_gate);
                let read_len = res.recruited[t].bases / res.recruited[t].reads.max(1);
                let size_settled = size.map_or(true, |s| {
                    let (before, grew) = (previous_size[t], last_size_growth[t]);
                    genome_settled(s, before, grew, cfg.saturation, read_len)
                });
                if !recount && reads_settled && size_settled {
                    saturated[t] = true;
                    stopped[t] = Some("saturated");
                }
                last_growth[t] = now.saturating_sub(previous[t]);
                previous[t] = now;
                if let Some(s) = size {
                    last_size_growth[t] = s.saturating_sub(previous_size[t]);
                    previous_size[t] = s;
                }
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
                let fraction = extension_fraction(
                    considered.reads,
                    considered.bases,
                    genome,
                    cfg.extend_depth,
                );
                let extend_gate = gate.as_ref().filter(|_| cfg.extend_gate);
                let floor = core_depth[t].map(|c| EXTEND_CORE_FLOOR * c);
                let ceiling = core_depth[t].map(|c| HIGH_COPY_FOLD * c);
                let bounds = (floor, ceiling);
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
                    if let (Some(table), Some(d)) = (kmer_table.as_ref(), kmer_target_depth[t]) {
                        if !kmer_extends(table, &bait, t, &r.seq, d) {
                            gated_out[t] += 1;
                            return;
                        }
                    }
                    let Some(g) = extend_gate else {
                        picked.push(r.seq);
                        return;
                    };
                    match extends(g, &bait, t, &r.seq, extension_deep(g, core_depth[t]), bounds) {
                        (Gauge::Extends, by_new) => {
                            new_part[t] += u64::from(by_new);
                            picked.push(r.seq);
                        }
                        (Gauge::Shallow, _) => shallow_out[t] += 1,
                        (Gauge::Deep, _) => deep_out[t] += 1,
                        (Gauge::Gated, _) => gated_out[t] += 1,
                    }
                })?;
                extended_from[t] = picked.len() as u64;
                added[t] = bait.extend_from_reads(
                    t,
                    picked.iter().map(Vec::as_slice),
                    cfg.extend_min_count,
                );
                used[t].extend(fresh);
                // Within-pass extension from the partial reads at its edges (`FRONTIER_STEPS`),
                // gated as extension is.
                let Some(g) = extend_gate.filter(|_| added[t] > 0) else {
                    continue;
                };
                let deep = extension_deep(g, core_depth[t]);
                let mut pending = res.frontier[t].reads();
                for _ in 0..FRONTIER_STEPS {
                    bait.prepare_lookup();
                    let covered: Vec<bool> = pending
                        .par_iter()
                        .map(|s| {
                            matches!(classify(s, &bait, &cfg.params).0,
                                Call::Organelle(mask) if mask & (1 << t) != 0)
                        })
                        .collect();
                    let (mut ready, mut rest) = (Vec::new(), Vec::new());
                    for (s, c) in pending.into_iter().zip(covered) {
                        if c {
                            ready.push(s);
                        } else {
                            rest.push(s);
                        }
                    }
                    pending = rest;
                    let judged: Vec<(&[u8], bool)> = ready
                        .into_par_iter()
                        .filter_map(|s| match extends(g, &bait, t, s, deep, bounds) {
                            (Gauge::Extends, by_new) => Some((s, by_new)),
                            _ => None,
                        })
                        .collect();
                    new_part[t] += judged.iter().filter(|j| j.1).count() as u64;
                    let picked: Vec<&[u8]> = judged.into_iter().map(|j| j.0).collect();
                    if picked.is_empty() {
                        break;
                    }
                    frontier_reads[t] += picked.len() as u64;
                    frontier_steps[t] += 1;
                    let more = bait.extend_from_reads(t, picked, cfg.extend_min_count);
                    added[t] += more;
                    if more == 0 {
                        break;
                    }
                }
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
                    gate = Some(build_gate(&cfg.inputs, replay.as_ref(), gate_sketch.take(), &mut warnings)?);
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
            tandem_reads: res.tandem,
            bait_kmers_added: added.clone(),
            extended_from_reads: extended_from,
            saturated: saturated.clone(),
            bootstrapped_from_reads: bootstrapped,
            extension_gated_out: gated_out,
            extension_shallow_out: shallow_out,
            extension_deep_out: deep_out,
            frontier_reads,
            frontier_steps,
            extension_new_part: new_part,
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
            None,
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
            if failed.is_ok() && (fraction >= 1.0 || unit_hash(sample_name(&r.id, cfg.mates_together), cfg.salt) < fraction) {
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
        tool: "ovasm find recruit",
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
    fn the_kmer_gate_lets_organelle_sequence_extend_and_keeps_the_nucleus_out() {
        // an organelle at 200x, its bait the first 1,000 bp; the nucleus at 5x, joined to the
        // organelle's 900-1,000 bp by a NUMT
        let organelle = pseudo_random(2000, 41);
        let nucleus = pseudo_random(1000, 42);
        let table = crate::frontier::KmerCounts::with_k(20, 21);
        let count = |seq: &[u8], copies: usize| {
            for _ in 0..copies {
                for (_, km) in CanonicalKmers::new(seq, 21) {
                    table.add(km);
                }
            }
        };
        count(&organelle, 200);
        count(&nucleus, 5);
        let mut bait = BaitSet::new(21);
        let t = bait.target("mitochondrion").unwrap();
        bait.add_sequence(t, &organelle[..1000]);
        bait.prepare_lookup();
        let reads: Vec<&[u8]> = (0..800).step_by(20).map(|a| &organelle[a..a + 150]).collect();
        let depth = kmer_depth(&table, &bait, t, &reads).unwrap();
        assert!((depth - 200.0).abs() < 10.0, "{depth}");
        // 120 bp on the bait and 30 bp further into the organelle: extends
        assert!(kmer_extends(&table, &bait, t, &organelle[880..1030], depth));
        // 120 bp on the bait and 30 bp into the nucleus past the NUMT: does not
        let numt_edge = [&organelle[880..1000], &nucleus[..30]].concat();
        assert!(!kmer_extends(&table, &bait, t, &numt_edge, depth));
    }

    #[test]
    fn mates_are_sampled_together_only_when_asked() {
        let (a, b) = (b"SRR13333747.1 1/1".as_slice(), b"SRR13333747.1 1/2".as_slice());
        assert_eq!(sample_name(a, true), sample_name(b, true));
        assert_eq!(unit_hash(sample_name(a, true), 3), unit_hash(sample_name(b, true), 3));
        assert_eq!(sample_name(a, false), a);
        assert_ne!(unit_hash(sample_name(a, false), 3), unit_hash(sample_name(b, false), 3));
    }

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

    #[test]
    fn reads_skipped_by_their_cover_bound_keep_the_calls_of_a_full_reclassification() {
        // a repeat in every "nuclear" read; bait grows by organelle pieces and by the repeat
        let organelle = pseudo_random(30_000, 701);
        let repeat = pseudo_random(60, 702);
        let mut reads: Vec<Vec<u8>> = Vec::new();
        for i in 0..60u64 {
            let mut r = pseudo_random(3000, 800 + i);
            let at = (i as usize * 37) % 2900;
            r[at..at + 60].copy_from_slice(&repeat);
            if i % 5 == 0 {
                // a NUMT-like read: 600 bp of organelle in nuclear flanks
                let o = (i as usize * 401) % 29_000;
                r[1000..1600].copy_from_slice(&organelle[o..o + 600]);
            }
            reads.push(r);
        }
        for i in 0..20 {
            reads.push(organelle[i * 1400..i * 1400 + 3000].to_vec());
        }
        let p = ClassifyParams {
            min_distinct: 0.0,
            min_cover: 0.8,
            min_hit_frac: 0.2,
            numt_min_block: 300,
            block_gap: 50,
        };
        let mut bait = bait_from(&organelle[..2000], 21);
        bait.target("plastid").unwrap();
        bait.track_changes();
        let (mut calls, mut bounds): (Vec<Call>, Vec<u16>) = reads
            .iter()
            .map(|r| {
                let (c, s) = classify(r, &bait, &p);
                (c, s.cover as u16)
            })
            .unzip();
        let mut skipped = 0;
        let pieces: Vec<(usize, Vec<u8>)> = (0..14)
            .map(|i| (0, organelle[2000 + i * 2000..4000 + i * 2000].to_vec()))
            .chain([(0, repeat[..40].to_vec()), (1, repeat.clone()), (1, organelle[100..900].to_vec())])
            .collect();
        for (target, piece) in pieces {
            bait.add_sequence(target, &piece);
            bait.prepare_lookup();
            let delta = bait.take_changes();
            for (i, r) in reads.iter().enumerate() {
                let fresh = classify(r, &bait, &p);
                if let Some(added) = added_cover(r, &delta, p.block_gap) {
                    match still_other(calls[i], bounds[i], added, r.len(), &p) {
                        Some(b) => {
                            skipped += 1;
                            assert!(b as usize >= fresh.1.cover, "bound {b} below cover {}", fresh.1.cover);
                            bounds[i] = b;
                        }
                        None => {
                            calls[i] = fresh.0;
                            bounds[i] = fresh.1.cover as u16;
                        }
                    }
                }
                assert_eq!(calls[i], fresh.0, "read {i}");
            }
        }
        assert!(skipped > 30, "the bound skipped only {skipped} reclassifications");
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
            min_distinct: 0.0,
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
                min_distinct: 0.0,
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
            min_distinct: 0.0,
            min_cover: 0.8,
            min_hit_frac: 0.05,
            numt_min_block: 5000,
            block_gap: 50,
        };
        assert_eq!(classify(&noisy, &bait, &bridged).0, Call::Organelle(1));
        let strict = ClassifyParams {
            min_distinct: 0.0,
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
            min_distinct: 0.0,
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
                min_distinct: 0.0,
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
            mates_together: false,
            kmer_gate: false,
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
                min_distinct: 0.0,
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
            mates_together: false,
            kmer_gate: false,
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
        assert_eq!(gated.depth_gate.as_ref().unwrap().source, "first pass and inputs");
        assert_eq!(
            cached.depth_gate.as_ref().unwrap().source,
            "first pass and recruitment cache"
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
                min_distinct: 0.0,
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
            mates_together: false,
            kmer_gate: false,
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
                min_distinct: 0.0,
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
            mates_together: false,
            kmer_gate: false,
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
    fn a_tandem_array_the_bait_covers_is_not_recruited() {
        // a satellite: one 1 kb unit along a whole 17 kb read
        let unit = pseudo_random(1000, 61);
        let array: Vec<u8> = unit.iter().cycle().take(17_000).copied().collect();
        assert!(distinct_share(&pseudo_random(17_000, 62), 21) > 0.95);
        let d = distinct_share(&array, 21);
        assert!(d < 0.1, "{d}");
        // two tandem copies of a 5 kb organelle segment in a 15 kb read stay above the cutoff
        let segment = pseudo_random(5000, 63);
        let duplication = [&segment[..], &segment[..], &pseudo_random(5000, 64)[..]].concat();
        let d = distinct_share(&duplication, 21);
        assert!(d > HIFI_MIN_DISTINCT && d < 0.8, "{d}");
        // a bait holding the unit (and its junction) covers the whole array
        let mut bait = BaitSet::new(21);
        let t = bait.target("mitochondrion").unwrap();
        bait.add_sequence(t, &[&unit[..], &unit[..]].concat());
        let mut p = params();
        assert_eq!(classify(&array, &bait, &p).0, Call::Organelle(1));
        p.min_distinct = HIFI_MIN_DISTINCT;
        let (call, score) = classify(&array, &bait, &p);
        assert_eq!(call, Call::Other);
        assert!(score.tandem);
    }

    #[test]
    fn a_target_settles_when_its_reads_and_its_genome_have() {
        // Medicago polymorpha's mitochondrion at its sixth pass: +1,755 reads after +4,382
        // (settled), its genome 131,310 to 169,392 bp after +39,696 (still growing)
        assert!(settled(50_356, 48_601, 4_382, 0.05, true));
        assert!(!settled(169_392, 131_310, 39_696, 0.05, true));
        // a close seed still stops at its second pass: the first pass's reads count as growth
        assert!(settled(10_300, 10_000, 10_000, 0.05, true));
        // no growth is not steady growth; ungated, only the 5% rule
        assert!(settled(10_000, 10_000, 0, 0.05, true));
        assert!(settled(10_300, 10_000, 300, 0.05, false));
        assert!(!settled(10_300, 10_000, 300, 0.05, true));
        // the Lophophora williamsii mitochondrion at its 14th pass with the extension of d06c6fe:
        // its genome 1,691,974 to 1,756,606 bp after +150,009, under half and under 5%, but eight
        // of its 8.1 kb reads; at the 26th of 0.3.2, +7,391 bp after +46,115
        assert!(settled(1_756_606, 1_691_974, 150_009, 0.05, true));
        assert!(!genome_settled(1_756_606, 1_691_974, 150_009, 0.05, 8_137));
        assert!(genome_settled(1_906_388, 1_898_997, 46_115, 0.05, 8_140));
    }

    #[test]
    fn extension_takes_at_least_min_extend_reads() {
        // the Schoenoplectus lacustris plastid's first pass: 289 reads at 5,355x for a "genome"
        // of 214 bp; all of them extend
        assert_eq!(extension_fraction(289, 4_380_000, 214.0, 200.0), 1.0);
        // its sixth: 9,115 reads, 139 Mb over 8,923 bp
        let f = extension_fraction(9_115, 139_255_348, 8_923.0, 200.0);
        assert!((f - 5000.0 / 9115.0).abs() < 1e-9, "{f}");
        // a deep plastid with many reads is still subsampled to about extend_depth
        let f = extension_fraction(391_044, 4_052_425_282, 120_971.0, 200.0);
        assert!((f - 5000.0 / 391_044.0).abs() < 1e-9, "{f}");
        // and to extend_depth when that keeps more than MIN_EXTEND_READS (20,000x of 2 Mb)
        let f = extension_fraction(4_000_000, 40_000_000_000, 2_000_000.0, 200.0);
        assert!((f - 0.01).abs() < 1e-9, "{f}");
        // shallow enough: everything
        assert_eq!(extension_fraction(100_000, 1_000_000_000, 10_000_000.0, 200.0), 1.0);
    }

    #[test]
    fn frontier_reads_extend_a_gated_target_within_a_pass() {
        // A 30 kb mitochondrion at 120x seeded with 4 kb of itself: one pass moves each edge by a
        // fraction of a read; from the pass's partial reads at the edges the gated target moves
        // them further within the pass, and recruits the same reads in far fewer passes.
        let mito = pseudo_random(30_000, 11);
        let nuclear = pseudo_random(300_000, 33);
        let dir = std::env::temp_dir().join(format!("ovasm-frontier-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("reads.fa");
        std::fs::write(&fa, circle_reads(&mito, 25, "m") + &linear_reads(&nuclear, 600, "n")).unwrap();
        let mut bait = BaitSet::new(21);
        let t = bait.target("mitochondrion").unwrap();
        bait.add_sequence(t, &mito[..4000]);
        let plain = run(bait.clone(), test_config(&fa, dir.join("plain"), 60, false)).unwrap();
        let gated = run(bait, test_config(&fa, dir.join("gated"), 60, true)).unwrap();
        let mito_reads = |r: &RecruitReport| count_reads(&r.outputs[0].output, b'm');
        assert_eq!(mito_reads(&plain), 1200);
        assert_eq!(mito_reads(&gated), 1200);
        assert_eq!(count_reads(&gated.outputs[0].output, b'n'), 0);
        assert!(gated.outputs[0].converged && plain.outputs[0].converged);
        let frontier: u64 = gated.iterations.iter().map(|it| it.frontier_reads[0]).sum();
        assert!(frontier > 0);
        assert!(
            3 * gated.iterations.len() < plain.iterations.len(),
            "passes: gated {}, plain {}",
            gated.iterations.len(),
            plain.iterations.len()
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    /// `len` bp reads every `step` bp along `genome` (around it when `circular`), named `tag<i>`.
    fn tiled(genome: &[u8], step: usize, len: usize, circular: bool, tag: &str) -> String {
        let (seq, last) = match circular {
            true => ([genome, &genome[..len]].concat(), genome.len()),
            false => (genome.to_vec(), genome.len() - len),
        };
        (0..last)
            .step_by(step)
            .enumerate()
            .map(|(i, st)| {
                let read = std::str::from_utf8(&seq[st..st + len]).unwrap();
                format!(">{tag}{i}\n{read}\n")
            })
            .collect()
    }

    #[test]
    fn extension_gets_past_a_tandem_array_its_bait_holds() {
        // A 52 kb mitochondrion at 200x in 20 kb reads: 20 kb, a 6 kb tandem array (40 copies of
        // 150 bp), 20 kb, the array again; seeded with the first 20 kb and the arrays, as the
        // Schoenoplectus lacustris plastid's bait held its inverted repeat's arrays. Every read
        // into the second 20 kb that the bait covers by 80% holds a whole array: 30% of its
        // k-mers, 80 times as deep as the rest and few distinct. Only the sequence it adds is
        // organelle-like, and from that the bait closes the circle.
        let unit = pseudo_random(150, 71);
        let array: Vec<u8> = unit.iter().cycle().take(6000).copied().collect();
        let (a, b) = (pseudo_random(20_000, 72), pseudo_random(20_000, 73));
        let mito = [&a[..], &array, &b, &array].concat();
        let nuclear = pseudo_random(300_000, 33);
        let dir = std::env::temp_dir().join(format!("ovasm-array-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("reads.fa");
        let text =
            tiled(&mito, 100, 20_000, true, "m") + &tiled(&nuclear, 4000, 20_000, false, "n");
        std::fs::write(&fa, text).unwrap();
        let mut bait = BaitSet::new(21);
        let t = bait.target("mitochondrion").unwrap();
        bait.add_sequence(t, &mito[..26_000]);
        bait.add_sequence(t, &mito[46_000..]);
        let mut config = test_config(&fa, dir.join("gated"), 60, true);
        config.params.min_distinct = HIFI_MIN_DISTINCT;
        let gated = run(bait, config).unwrap();
        assert_eq!(count_reads(&gated.outputs[0].output, b'm'), 520);
        assert_eq!(count_reads(&gated.outputs[0].output, b'n'), 0);
        let by_new: u64 = gated.iterations.iter().map(|it| it.extension_new_part[0]).sum();
        assert!(by_new > 0);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn extension_gets_through_a_stretch_of_its_target_near_the_organelle_threshold() {
        // A 30 kb mitochondrion at 120x beside a nucleus at 24x (organelle threshold 72x), 3 kb of
        // it missing from half of its molecules and so at 60x, as Lophophora williamsii's had two
        // stretches at 66-90x (threshold 66x); seeded away from that stretch. Its k-mers count as
        // deep from halfway between the nucleus and the target (54x), and the reads across it
        // extend the bait; the nuclear reads stay out.
        let full = pseudo_random(30_000, 74);
        let short = [&full[..15_000], &full[18_000..]].concat();
        let nuclear = pseudo_random(300_000, 33);
        let dir = std::env::temp_dir().join(format!("ovasm-near-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("reads.fa");
        let text = circle_reads(&full, 50, "f")
            + &circle_reads(&short, 50, "s")
            + &linear_reads(&nuclear, 125, "n");
        std::fs::write(&fa, text).unwrap();
        let mut bait = BaitSet::new(21);
        let t = bait.target("mitochondrion").unwrap();
        bait.add_sequence(t, &full[..6000]);
        let gated = run(bait, test_config(&fa, dir.join("gated"), 60, true)).unwrap();
        assert_eq!(count_reads(&gated.outputs[0].output, b'f'), 600);
        assert_eq!(count_reads(&gated.outputs[0].output, b's'), 540);
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
        let mut forward = Sample::<DEPTH_SAMPLE>::default();
        let mut backward = Sample::<DEPTH_SAMPLE>::default();
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
