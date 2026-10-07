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
}

/// Counts of one pass; its recruited reads went to the pass files (see `scan`).
struct ScanResult {
    recruited: Vec<TargetCounts>,
    numt: TargetCounts,
    reads: u64,
    bases: u64,
    multi_target: u64,
}

/// A target's recruited reads of the current pass, in input order.
fn pass_path(out_dir: &Path, target: &str) -> PathBuf {
    out_dir.join(format!(".pass.{target}.fastq"))
}

/// Every record of a FASTQ written by `write_rec`, in order. A pass that recruited nothing
/// leaves an empty file, which the parser rejects, so it is read as no records.
fn for_each_rec(path: &Path, mut f: impl FnMut(Rec)) -> Result<()> {
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

/// One pass over the inputs. Recruited reads are written, in input order, to each target's
/// pass file and NUMT candidates to `numt_path`, not kept: a plastid at 20,000x recruits
/// gigabases (Col-0 HiFi: 3.2 Gb, 14 GB of memory when held), more than a laptop has.
fn scan(
    inputs: &[PathBuf],
    bait: &BaitSet,
    p: &ClassifyParams,
    out_dir: &Path,
    numt_path: &Path,
) -> Result<ScanResult> {
    let (rx, reader) = read_batches(inputs);
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
    let mut out = ScanResult {
        recruited: vec![TargetCounts::default(); bait.targets.len()],
        numt: TargetCounts::default(),
        reads: 0,
        bases: 0,
        multi_target: 0,
    };
    for batch in rx {
        out.reads += batch.len() as u64;
        out.bases += batch.iter().map(|r| r.seq.len() as u64).sum::<u64>();
        let calls: Vec<(Call, Rec)> = batch
            .into_par_iter()
            .filter_map(|rec| {
                let (call, _) = classify(&rec.seq, bait, p);
                (call != Call::Other).then_some((call, rec))
            })
            .collect();
        for (call, rec) in calls {
            let len = rec.seq.len() as u64;
            match call {
                Call::Organelle(mask) => {
                    if mask.count_ones() > 1 {
                        out.multi_target += 1;
                    }
                    for t in (0..8).filter(|t| mask & (1 << t) != 0) {
                        write_rec(&mut writers[t], &rec)?;
                        out.recruited[t].reads += 1;
                        out.recruited[t].bases += len;
                    }
                }
                Call::NumtCandidate(_) => {
                    write_rec(&mut numt_w, &rec)?;
                    out.numt.reads += 1;
                    out.numt.bases += len;
                }
                Call::Other => {}
            }
        }
    }
    reader.join().expect("reader thread panicked")?;
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
    pub sampling_fraction: f64,
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
    pub iterations_requested: usize,
    /// Passes actually run: iteration stops early once no target's bait grows.
    pub iterations_run: usize,
    pub extend_min_count: u32,
    pub extend_depth: f64,
    pub saturation: f64,
    pub target_depth: Option<f64>,
    pub salt: u64,
    pub inputs: Vec<String>,
    pub targets: Vec<crate::bait::Target>,
    pub iterations: Vec<IterationStats>,
    pub outputs: Vec<TargetOutput>,
    pub numt_candidates: u64,
    pub numt_output: String,
    pub bootstrap: bool,
    pub extend_gate: bool,
    /// Depth sketch of the bootstrap or extension gate, when one was built.
    pub bootstrap_nuclear_peak_depth: Option<u32>,
    pub bootstrap_organelle_min_depth: Option<f64>,
    pub elapsed_seconds: f64,
    pub peak_rss_mb: Option<f64>,
}

pub struct RecruitConfig {
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
    pub extend_gate: bool,
}

/// A bootstrap read must not be this much shallower than the target's own reads (when the
/// target has some): a mitochondrial read carrying a plastid insertion is ~10-30x shallower
/// than plastid reads, so it cannot bootstrap the plastid, and the reverse needs the other
/// target to be complete (see `bootstrap_candidates`).
const BOOTSTRAP_DEPTH_FLOOR: f64 = 0.5;
/// Own reads whose depth sets that floor (a sample is enough).
const BOOTSTRAP_DEPTH_SAMPLE: usize = 2000;
/// Fewest own reads for a depth floor.
const BOOTSTRAP_DEPTH_MIN_READS: usize = 10;
/// A bootstrap read must not be covered this much by another target's bait.
const BOOTSTRAP_MAX_OTHER_COVER: f64 = 0.5;

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
/// and not be mostly covered by another target; and it is only used once every other target's
/// bait has stopped growing, so their own reads are recruited, not left as candidates.
fn bootstrap_candidates(
    numt_path: &Path,
    own_path: &Path,
    t: usize,
    bait: &BaitSet,
    p: &ClassifyParams,
    gate: &DepthGate,
    used: &FxHashSet<u64>,
) -> Result<(Vec<Vec<u8>>, Vec<u64>)> {
    let mut own: Vec<f64> = Vec::new();
    for_each_rec(own_path, |r| {
        if own.len() < BOOTSTRAP_DEPTH_SAMPLE {
            own.extend(gate.median_depth(&r.seq));
        }
    })?;
    own.sort_by(f64::total_cmp);
    let floor = if own.len() >= BOOTSTRAP_DEPTH_MIN_READS {
        BOOTSTRAP_DEPTH_FLOOR * own[own.len() / 2]
    } else {
        0.0
    };
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

pub fn run(mut bait: BaitSet, cfg: RecruitConfig) -> Result<RecruitReport> {
    let start = Instant::now();
    std::fs::create_dir_all(&cfg.out_dir)?;
    let mut iterations = Vec::new();
    let mut result = None;
    let n_targets = bait.targets.len();
    let numt_path = cfg.out_dir.join("numt_candidates.fastq");
    // Reads already used to extend each target's bait: later passes count only new reads, so a
    // deep, already complete target (a plastid at 20,000x) is not re-counted every pass.
    // (64-bit hashes of the read names: millions of names, 8 bytes each instead of a string.)
    let mut used: Vec<FxHashSet<u64>> = vec![FxHashSet::default(); n_targets];
    // A target seeded from reads of the sample holds its whole genome already: extending it
    // only adds passes over the input (each a full decompression) for nothing.
    let mut saturated: Vec<bool> = bait.targets.iter().map(|t| t.from_reads).collect();
    let mut previous = vec![0u64; n_targets];
    // Bootstrap depth sketch, built the first time a target needs it.
    let mut gate: Option<DepthGate> = None;
    let mut boot_used: Vec<FxHashSet<u64>> = vec![FxHashSet::default(); n_targets];
    for it in 0..cfg.iterations.max(1) {
        let t0 = Instant::now();
        let bait_kmers = bait.len();
        let res = scan(&cfg.inputs, &bait, &cfg.params, &cfg.out_dir, &numt_path)?;
        let last = it + 1 == cfg.iterations.max(1);
        let mut added = vec![0; n_targets];
        let mut extended_from = vec![0u64; n_targets];
        let mut gated_out = vec![0u64; n_targets];
        if cfg.extend_gate && !last && gate.is_none() {
            gate = Some(DepthGate::long_reads(&cfg.inputs)?);
        }
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
                if !recount
                    && previous[t] > 0
                    && (now.saturating_sub(previous[t]) as f64)
                        < cfg.saturation * previous[t] as f64
                {
                    saturated[t] = true;
                }
                previous[t] = now;
                if saturated[t] {
                    continue;
                }
                // New reads of this pass, from its file: their depth sets the subsample
                // fraction, then the subsample is read back (at most `extend_depth`).
                let path = pass_path(&cfg.out_dir, &bait.targets[t].name);
                let (mut fresh, mut fresh_bases) = (Vec::new(), 0u64);
                for_each_rec(&path, |r| {
                    let key = name_key(&r.id);
                    if recount || !used[t].contains(&key) {
                        fresh.push(key);
                        fresh_bases += r.seq.len() as u64;
                    }
                })?;
                let seed_len = bait.targets[t].seed_length.max(1) as f64;
                let depth = fresh_bases as f64 / seed_len;
                let fraction = if depth > cfg.extend_depth {
                    cfg.extend_depth / depth
                } else {
                    1.0
                };
                let mut picked: Vec<Vec<u8>> = Vec::new();
                let extend_gate = gate.as_ref().filter(|_| cfg.extend_gate);
                for_each_rec(&path, |r| {
                    if (recount || !used[t].contains(&name_key(&r.id)))
                        && (fraction >= 1.0 || unit_hash(&r.id, EXTEND_SALT) < fraction)
                    {
                        if extend_gate.is_some_and(|g| !g.organelle_like(&r.seq)) {
                            gated_out[t] += 1;
                        } else {
                            picked.push(r.seq);
                        }
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
                {
                    continue;
                }
                if gate.is_none() {
                    let g = DepthGate::long_reads(&cfg.inputs)?;
                    eprintln!(
                        "[ovasm] bootstrap depth gate: nuclear peak {}x, organelle threshold {:.0}x (sketch depth)",
                        g.nuclear_peak_depth, g.organelle_min_depth
                    );
                    gate = Some(g);
                }
                let (picked, keys) = bootstrap_candidates(
                    &numt_path,
                    &pass_path(&cfg.out_dir, &bait.targets[t].name),
                    t,
                    &bait,
                    &cfg.params,
                    gate.as_ref().expect("built above"),
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
                // its bait changed: extend it again from what the next pass recruits
                saturated[t] = false;
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
        });
        result = Some(res);
        // Nothing new to bait with: this pass already used the final bait set.
        if !last && added.iter().all(|&a| a == 0) {
            break;
        }
    }
    let res = result.expect("at least one iteration");

    // Final outputs from the last pass's files, downsampled to `target_depth`.
    let mut outputs = Vec::new();
    for (t, counts) in res.recruited.iter().enumerate() {
        let target = &bait.targets[t];
        let depth = if target.seed_length > 0 {
            counts.bases as f64 / target.seed_length as f64
        } else {
            0.0
        };
        let fraction = match cfg.target_depth {
            Some(d) if depth > d => d / depth,
            _ => 1.0,
        };
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
        outputs.push(TargetOutput {
            name: target.name.clone(),
            recruited_reads: counts.reads,
            recruited_bases: counts.bases,
            seed_length: target.seed_length,
            estimated_depth: depth,
            sampling_fraction: fraction,
            kept_reads,
            kept_bases,
            output: path.display().to_string(),
        });
    }

    let report = RecruitReport {
        tool: "ovasm recruit",
        version: env!("CARGO_PKG_VERSION"),
        k: bait.k,
        params: cfg.params.clone(),
        iterations_requested: cfg.iterations,
        iterations_run: iterations.len(),
        extend_min_count: cfg.extend_min_count,
        extend_depth: cfg.extend_depth,
        saturation: cfg.saturation,
        target_depth: cfg.target_depth,
        salt: cfg.salt,
        inputs: cfg.inputs.iter().map(|p| p.display().to_string()).collect(),
        targets: bait.targets.clone(),
        iterations,
        outputs,
        numt_candidates: res.numt.reads,
        numt_output: numt_path.display().to_string(),
        bootstrap: cfg.bootstrap,
        extend_gate: cfg.extend_gate,
        bootstrap_nuclear_peak_depth: gate.as_ref().map(|g| g.nuclear_peak_depth),
        bootstrap_organelle_min_depth: gate.as_ref().map(|g| g.organelle_min_depth),
        elapsed_seconds: start.elapsed().as_secs_f64(),
        peak_rss_mb: peak_rss_mb(),
    };
    std::fs::write(
        cfg.out_dir.join("recruit.json"),
        serde_json::to_string_pretty(&report)?,
    )?;
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
        let boot = run(bait, cfg(true)).unwrap();
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
        let gated = run(bait, cfg(true)).unwrap();
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
        std::fs::write(&reads, fasta).unwrap();
        let mut bait = bait_from(&genome[..5000], 21);
        bait.targets[0].seed_length = genome.len() as u64;
        let cfg = |iterations| RecruitConfig {
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
        };
        let one = run(bait.clone(), cfg(1)).unwrap();
        let many = run(bait, cfg(40)).unwrap();
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
}
