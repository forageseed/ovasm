//! Stage 4 (S4.2): self-correction of noisy organelle long reads (`ovasm correct`).
//!
//! Organelle reads are hundreds to thousands of times deeper than the genome they come from, so
//! even at 85-93% read identity every true k-mer of a small k is seen many times while error
//! k-mers are mostly seen once. Solid k-mers (above the valley after the error peak of the count
//! histogram) form a de Bruijn graph of the organelle. Each read is corrected against it:
//!
//! - maximal runs of consecutive solid k-mers in the read are anchors;
//! - between two anchors the graph is searched for paths from the last k-mer of one to the first
//!   of the next whose length is close to the read's distance between them; the path closest to
//!   the read segment (edit distance) replaces it; an unbridged gap keeps the read's bases;
//! - read ends outside the first and last anchor are trimmed (nothing vouches for them).
//!
//! Rounds with increasing k (e.g. 17, then 25 on the corrected reads) resolve more repeats as the
//! reads get cleaner. Only counts >= 2 are stored: a Bloom filter absorbs first sightings.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Result};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::Serialize;

use crate::recruit::{peak_rss_mb, read_batches, Rec};

#[derive(Clone, Copy)]
pub struct CorrectParams {
    pub k: usize,
    /// Solid-count threshold; `None` picks the valley after the error peak.
    pub min_count: Option<u32>,
    /// Allowed relative length difference between a bridging path and the read gap.
    pub length_tolerance: f64,
    /// Longest read gap (bp) searched for a bridge.
    pub max_gap: usize,
    /// Search budget (graph steps) per gap.
    pub max_steps: usize,
    pub keep_ends: bool,
    /// A k-mer counts as solid only if its count reaches this fraction of its best sibling
    /// (same k-mer with a different first or last base): systematic errors such as homopolymer
    /// length calls recur in many reads but stay a minority beside the true k-mer.
    pub relative: f64,
    /// Auto threshold: at least this fraction of the solid k-mers' peak count. Systematic
    /// errors recur in many reads and pass the histogram valley in later rounds.
    pub peak_fraction: f64,
}

#[inline]
fn code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

const BASES: [u8; 4] = *b"ACGT";

/// Reverse complement of a 2-bit packed k-mer.
#[inline]
fn rc_kmer(mut v: u64, k: usize) -> u64 {
    let mut r = 0u64;
    for _ in 0..k {
        r = (r << 2) | (3 - (v & 3));
        v >>= 2;
    }
    r
}

#[inline]
fn canon(v: u64, k: usize) -> u64 {
    v.min(rc_kmer(v, k))
}

/// Forward 2-bit k-mers of a read with their start positions (windows with non-ACGT skipped).
fn forward_kmers(seq: &[u8], k: usize) -> Vec<(usize, u64)> {
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let mut out = Vec::with_capacity(seq.len());
    let (mut v, mut valid) = (0u64, 0usize);
    for (i, &b) in seq.iter().enumerate() {
        match code(b) {
            Some(c) => {
                v = ((v << 2) | c) & mask;
                valid += 1;
                if valid >= k {
                    out.push((i + 1 - k, v));
                }
            }
            None => valid = 0,
        }
    }
    out
}

#[inline]
fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Bloom filter holding k-mers seen at least once.
pub(crate) struct Bloom {
    bits: Vec<u64>,
    mask: u64,
}

impl Bloom {
    pub(crate) fn new(log2_bits: u32) -> Self {
        Self {
            bits: vec![0; 1 << (log2_bits - 6)],
            mask: (1u64 << log2_bits) - 1,
        }
    }

    /// Insert; returns whether the k-mer was (probably) present before.
    pub(crate) fn insert(&mut self, x: u64) -> bool {
        let h = mix64(x);
        let mut present = true;
        for i in 0..3 {
            let bit = (h.rotate_left(21 * i) ^ (h >> 7)) & self.mask;
            let (w, b) = ((bit >> 6) as usize, bit & 63);
            if self.bits[w] >> b & 1 == 0 {
                present = false;
                self.bits[w] |= 1 << b;
            }
        }
        present
    }
}

/// Solid-count threshold: the first local minimum of the count histogram (from 2 up), the
/// valley between the error peak and the organelle peak.
fn valley(hist: &[u64]) -> u32 {
    for c in 3..hist.len().saturating_sub(1) {
        if hist[c] <= hist[c - 1] && hist[c] < hist[c + 1] {
            return c as u32;
        }
    }
    3
}

/// Count of `v` if it passes both the absolute threshold and the relative one against its
/// siblings (first or last base changed); `None` otherwise.
fn trusted(counts: &FxHashMap<u64, u32>, v: u64, k: usize, min: u32, relative: f64) -> bool {
    let c = match counts.get(&canon(v, k)) {
        Some(&c) if c >= min => c,
        _ => return false,
    };
    if relative <= 0.0 {
        return true;
    }
    let shift = 2 * (k - 1);
    let mut best = c;
    for b in 0..4u64 {
        let last = (v & !3) | b;
        let first = (v & !(3 << shift)) | (b << shift);
        for s in [last, first] {
            if s != v {
                if let Some(&n) = counts.get(&canon(s, k)) {
                    best = best.max(n);
                }
            }
        }
    }
    c as f64 >= relative * best as f64
}

/// Bounded depth-first search from `from` to `to` in the solid graph for paths of `lo..=hi`
/// steps; successors agreeing with the read's next base are tried first. Returns the spelled
/// bases (one per step) of up to `max_paths` paths.
#[allow(clippy::too_many_arguments)]
fn bridges(
    counts: &FxHashMap<u64, u32>,
    min: u32,
    relative: f64,
    k: usize,
    from: u64,
    to: u64,
    lo: usize,
    hi: usize,
    read_gap: &[u8],
    max_steps: usize,
) -> Vec<Vec<u8>> {
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let mut out = Vec::new();
    let mut steps = 0usize;
    // stack of (k-mer, depth, next successor index to try); path holds spelled bases
    let mut stack: Vec<(u64, usize, u8)> = vec![(from, 0, 0)];
    let mut path: Vec<u8> = Vec::new();
    while let Some(top) = stack.last_mut() {
        let (v, depth, tried) = *top;
        if tried == 4 || depth >= hi {
            stack.pop();
            path.pop();
            continue;
        }
        top.2 += 1;
        // try the read's base first, then the others in a fixed order
        let want = read_gap.get(depth).and_then(|&b| code(b)).unwrap_or(0);
        let c = (want + tried as u64) % 4;
        let next = ((v << 2) | c) & mask;
        if !trusted(counts, next, k, min, relative) {
            continue;
        }
        steps += 1;
        if steps > max_steps || out.len() >= 8 {
            break;
        }
        path.push(BASES[c as usize]);
        if next == to && depth + 1 >= lo {
            out.push(path.clone());
        }
        stack.push((next, depth + 1, 0));
    }
    out
}

/// Levenshtein distance (two rows).
fn edit_distance(a: &[u8], b: &[u8]) -> usize {
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    let mut cur = vec![0; b.len() + 1];
    for i in 1..=a.len() {
        cur[0] = i;
        for j in 1..=b.len() {
            let sub = prev[j - 1] + usize::from(a[i - 1] != b[j - 1]);
            cur[j] = sub.min(prev[j] + 1).min(cur[j - 1] + 1);
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    prev[b.len()]
}

#[derive(Default, Clone, Copy, Serialize)]
pub struct ReadStats {
    pub bases_in: u64,
    pub bases_out: u64,
    pub anchored_bases: u64,
    pub gaps: u64,
    pub gaps_bridged: u64,
    pub trimmed_bases: u64,
}

impl std::ops::AddAssign for ReadStats {
    fn add_assign(&mut self, o: Self) {
        self.bases_in += o.bases_in;
        self.bases_out += o.bases_out;
        self.anchored_bases += o.anchored_bases;
        self.gaps += o.gaps;
        self.gaps_bridged += o.gaps_bridged;
        self.trimmed_bases += o.trimmed_bases;
    }
}

/// Correct one read; `None` when it has no solid anchor at all.
fn correct_read(
    seq: &[u8],
    counts: &FxHashMap<u64, u32>,
    min_count: u32,
    p: &CorrectParams,
) -> (Option<Vec<u8>>, ReadStats) {
    let is_solid = |v: u64| trusted(counts, v, p.k, min_count, p.relative);
    // a true k-mer can fall below the solid threshold by chance; a path may pass weak ones
    let weak = (min_count / 3).max(2);
    let k = p.k;
    let mut st = ReadStats {
        bases_in: seq.len() as u64,
        ..Default::default()
    };
    let km = forward_kmers(seq, k);
    // maximal runs of consecutive solid k-mers: (first index, last index) into `km`
    let mut runs: Vec<(usize, usize)> = Vec::new();
    for (i, &(pos, v)) in km.iter().enumerate() {
        if !is_solid(v) {
            continue;
        }
        match runs.last_mut() {
            Some(r) if km[r.1].0 + 1 == pos && r.1 + 1 == i => r.1 = i,
            _ => runs.push((i, i)),
        }
    }
    if runs.is_empty() {
        st.trimmed_bases = seq.len() as u64;
        return (None, st);
    }
    let (first_pos, last_end) = (km[runs[0].0].0, km[runs[runs.len() - 1].1].0 + k);
    let mut out: Vec<u8> = Vec::with_capacity(seq.len());
    if p.keep_ends {
        out.extend_from_slice(&seq[..first_pos]);
    } else {
        st.trimmed_bases += first_pos as u64;
    }
    // first anchor in full
    out.extend_from_slice(&seq[first_pos..km[runs[0].1].0 + k]);
    st.anchored_bases += (km[runs[0].1].0 + k - first_pos) as u64;
    for w in runs.windows(2) {
        let (a, b) = (w[0], w[1]);
        let (apos, _) = km[a.1];
        let (bpos, _) = km[b.0];
        let (bend, _) = km[b.1];
        st.gaps += 1;
        // The k-mers at an anchor's edge can be solid errors (a homopolymer deletion yields the
        // same wrong k-mer in many reads), leaving no path between the anchors: retry with the
        // anchors pulled back so the bridge also replaces their edge.
        let mut bridged: Option<(usize, usize, Vec<u8>)> = None;
        let mut tried: Vec<(usize, usize, u32)> = Vec::new();
        let attempts = [
            (0, 0, min_count),
            (3, 0, min_count),
            (0, 3, min_count),
            (3, 3, min_count),
            (8, 8, min_count),
            (k, k, min_count),
            (0, 0, weak),
            (3, 3, weak),
        ];
        for (sa, sb, min) in attempts {
            let (sa, sb) = (sa.min(a.1 - a.0), sb.min(b.1 - b.0));
            if tried.contains(&(sa, sb, min)) {
                continue;
            }
            tried.push((sa, sb, min));
            let (ia, ib) = (a.1 - sa, b.0 + sb);
            let (pa, va) = km[ia];
            let (pb, vb) = km[ib];
            let dist = pb - pa;
            if dist > p.max_gap {
                continue;
            }
            let gap = &seq[pa + k..pb + k];
            // short gaps: a few clustered indels shift the length by more than a fraction
            let slack = ((dist as f64 * p.length_tolerance).ceil() as usize).max(8);
            let cands = bridges(
                counts,
                min,
                p.relative,
                k,
                va,
                vb,
                dist.saturating_sub(slack).max(1),
                dist + slack,
                gap,
                p.max_steps,
            );
            if let Some(best) = cands.into_iter().min_by_key(|c| edit_distance(c, gap)) {
                bridged = Some((ia, ib, best));
                break;
            }
        }
        match bridged {
            Some((ia, ib, best)) => {
                // drop the part of anchor A after k-mer `ia`, already written
                out.truncate(out.len() - (apos - km[ia].0));
                out.extend_from_slice(&best);
                out.extend_from_slice(&seq[km[ib].0 + k..bend + k]);
                st.gaps_bridged += 1;
            }
            None => {
                out.extend_from_slice(&seq[apos + k..bpos + k]);
                out.extend_from_slice(&seq[bpos + k..bend + k]);
            }
        }
        st.anchored_bases += (bend + k - bpos) as u64;
    }
    if p.keep_ends {
        out.extend_from_slice(&seq[last_end..]);
    } else {
        st.trimmed_bases += (seq.len() - last_end) as u64;
    }
    st.bases_out = out.len() as u64;
    (Some(out), st)
}

#[derive(Serialize)]
pub struct RoundReport {
    pub k: usize,
    pub distinct_kmers_ge2: usize,
    pub min_count: u32,
    pub min_count_auto: bool,
    pub solid_kmers: usize,
    /// Count with the most solid k-mers (the organelle's k-mer depth).
    pub solid_peak: u32,
    /// Solid k-mers came from short reads.
    pub hybrid: bool,
    pub reads_in: u64,
    pub reads_out: u64,
    pub stats: ReadStats,
    pub seconds: f64,
}

#[derive(Serialize)]
pub struct CorrectReport {
    pub rounds: Vec<RoundReport>,
    pub output: String,
    pub elapsed_seconds: f64,
    pub peak_rss_mb: Option<f64>,
}

fn load(inputs: &[PathBuf]) -> Result<Vec<Rec>> {
    let (rx, reader) = read_batches(inputs);
    let mut reads = Vec::new();
    for batch in rx {
        reads.extend(batch);
    }
    reader.join().expect("reader thread panicked")?;
    Ok(reads)
}

/// Canonical k-mers seen at least twice in `reads` (a Bloom filter absorbs first sightings).
fn count_kmers(reads: &[Rec], k: usize) -> FxHashMap<u64, u32> {
    let mut bloom = Bloom::new(30);
    let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
    for chunk in reads.chunks(4096) {
        let kms: Vec<Vec<u64>> = chunk
            .par_iter()
            .map(|r| {
                forward_kmers(&r.seq, k)
                    .into_iter()
                    .map(|(_, v)| canon(v, k))
                    .collect()
            })
            .collect();
        for v in kms.into_iter().flatten() {
            if let Some(c) = counts.get_mut(&v) {
                *c += 1;
            } else if bloom.insert(v) {
                counts.insert(v, 2);
            }
        }
    }
    counts
}

/// One correction round. Solid k-mers come from `short` reads when given (hybrid correction:
/// accurate reads with an unrelated error profile), otherwise from the long reads themselves.
fn round(
    reads: Vec<Rec>,
    short: Option<&[Rec]>,
    p: &CorrectParams,
) -> Result<(Vec<Rec>, RoundReport)> {
    let t0 = Instant::now();
    let k = p.k;
    let counts = count_kmers(short.unwrap_or(&reads), k);
    let mut hist = vec![0u64; 10_001];
    for &c in counts.values() {
        hist[(c as usize).min(10_000)] += 1;
    }
    let low = valley(&hist);
    let solid_peak = (low as usize..hist.len())
        .max_by_key(|&c| hist[c])
        .unwrap_or(0) as u32;
    let min_count = p
        .min_count
        .unwrap_or_else(|| low.max((p.peak_fraction * solid_peak as f64) as u32));
    let distinct = counts.len();
    let solid_kmers = counts.values().filter(|&&c| c >= min_count).count();
    let results: Vec<(Option<Vec<u8>>, ReadStats)> = reads
        .par_iter()
        .map(|r| correct_read(&r.seq, &counts, min_count, p))
        .collect();
    let mut stats = ReadStats::default();
    let reads_in = reads.len() as u64;
    let mut out = Vec::with_capacity(reads.len());
    for (r, (seq, s)) in reads.into_iter().zip(results) {
        stats += s;
        if let Some(seq) = seq.filter(|s| s.len() >= 2 * k) {
            out.push(Rec {
                id: r.id,
                seq,
                qual: None,
            });
        }
    }
    let report = RoundReport {
        k,
        distinct_kmers_ge2: distinct,
        min_count,
        min_count_auto: p.min_count.is_none(),
        solid_kmers,
        solid_peak,
        hybrid: short.is_some(),
        reads_in,
        reads_out: out.len() as u64,
        stats,
        seconds: t0.elapsed().as_secs_f64(),
    };
    Ok((out, report))
}

/// Correct `inputs` in rounds of k `ks`: against the short reads' solid k-mers when
/// `short_inputs` are given (hybrid), else against the reads' own. After hybrid rounds, `self_ks`
/// rounds use the reads' own k-mers: where the short reads leave a stretch uncovered (their
/// recruitment missed it; Nipponbare 5.2-6.3 kb had no Illumina read at 864x elsewhere) the long
/// reads keep their raw errors and the assembly breaks there. Stretches the hybrid rounds
/// corrected are already consistent across reads and stay as they are.
pub fn run(
    inputs: Vec<PathBuf>,
    short_inputs: Vec<PathBuf>,
    ks: &[usize],
    self_ks: &[usize],
    base: CorrectParams,
    out_fasta: &Path,
    out_json: &Path,
) -> Result<CorrectReport> {
    if ks.is_empty() || ks.iter().chain(self_ks).any(|&k| !(11..=31).contains(&k)) {
        bail!("each --k and --self-k must be in 11..=31");
    }
    let t0 = Instant::now();
    let mut reads = load(&inputs)?;
    let short = if short_inputs.is_empty() {
        None
    } else {
        Some(load(&short_inputs)?)
    };
    let mut rounds = Vec::new();
    let plan: Vec<(usize, Option<&[Rec]>)> = ks
        .iter()
        .map(|&k| (k, short.as_deref()))
        .chain(
            self_ks
                .iter()
                .filter(|_| short.is_some())
                .map(|&k| (k, None)),
        )
        .collect();
    for (k, solid_from) in plan {
        let (next, rep) = round(reads, solid_from, &CorrectParams { k, ..base })?;
        eprintln!(
            "[ovasm] correct k={}{}: solid >= {} ({} k-mers, peak {}); {} -> {} reads; {:.1}% anchored, {}/{} gaps bridged, {:.1}% trimmed; {:.1}s",
            k, if rep.hybrid { " (hybrid)" } else { "" }, rep.min_count, rep.solid_kmers, rep.solid_peak, rep.reads_in, rep.reads_out,
            100.0 * rep.stats.anchored_bases as f64 / rep.stats.bases_in.max(1) as f64,
            rep.stats.gaps_bridged, rep.stats.gaps,
            100.0 * rep.stats.trimmed_bases as f64 / rep.stats.bases_in.max(1) as f64,
            rep.seconds
        );
        rounds.push(rep);
        reads = next;
    }
    if let Some(dir) = out_fasta.parent() {
        std::fs::create_dir_all(dir)?;
    }
    let mut w = BufWriter::new(File::create(out_fasta)?);
    for r in &reads {
        w.write_all(b">")?;
        w.write_all(&r.id)?;
        w.write_all(b"\n")?;
        w.write_all(&r.seq)?;
        w.write_all(b"\n")?;
    }
    w.flush()?;
    let report = CorrectReport {
        rounds,
        output: out_fasta.display().to_string(),
        elapsed_seconds: t0.elapsed().as_secs_f64(),
        peak_rss_mb: peak_rss_mb(),
    };
    std::fs::write(out_json, serde_json::to_string_pretty(&report)?)?;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(n: usize, seed: u64) -> Vec<u8> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                BASES[(s >> 62) as usize]
            })
            .collect()
    }

    /// Copy `truth` with substitutions, insertions and deletions at `rate` (deterministic).
    fn noisy(truth: &[u8], rate: f64, seed: u64) -> Vec<u8> {
        let mut s = seed;
        let mut rnd = || {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (s >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut out = Vec::with_capacity(truth.len());
        for &b in truth {
            let x = rnd();
            if x < rate / 3.0 {
                out.push(BASES[((rnd() * 4.0) as usize).min(3)]); // substitution (maybe same)
            } else if x < 2.0 * rate / 3.0 {
                out.push(b);
                out.push(BASES[((rnd() * 4.0) as usize).min(3)]); // insertion
            } else if x < rate {
                // deletion
            } else {
                out.push(b);
            }
        }
        out
    }

    fn revcomp(s: &[u8]) -> Vec<u8> {
        crate::evidence::revcomp(s)
    }

    #[test]
    fn rc_and_canonical_kmers_agree_with_sequence() {
        let s = b"ACGTTGCAAGGCTTAGC";
        let k = s.len();
        let (_, f) = forward_kmers(s, k)[0];
        let (_, r) = forward_kmers(&revcomp(s), k)[0];
        assert_eq!(rc_kmer(f, k), r);
        assert_eq!(canon(f, k), canon(r, k));
    }

    #[test]
    fn valley_is_the_first_local_minimum() {
        let hist = [0, 0, 900, 300, 40, 12, 20, 60, 90, 70];
        assert_eq!(valley(&hist), 5);
    }

    /// Edit distance of `q` against the best-matching substring of `t` (free gaps at both ends of t).
    fn semi_global(q: &[u8], t: &[u8]) -> usize {
        let mut prev = vec![0usize; t.len() + 1];
        let mut cur = vec![0usize; t.len() + 1];
        for i in 1..=q.len() {
            cur[0] = i;
            for j in 1..=t.len() {
                let sub = prev[j - 1] + usize::from(q[i - 1] != t[j - 1]);
                cur[j] = sub.min(prev[j] + 1).min(cur[j - 1] + 1);
            }
            std::mem::swap(&mut prev, &mut cur);
        }
        *prev.iter().min().unwrap()
    }

    #[test]
    fn noisy_reads_come_back_close_to_the_truth() {
        // 8 kb circular genome, 250 reads of 3 kb at 10% error (subs + indels) on both strands.
        let g = seq(8000, 5);
        let circ = [g.clone(), g[..3000].to_vec()].concat();
        let mut reads = Vec::new();
        let mut truths = Vec::new();
        for i in 0..250u64 {
            let st = ((i * 997) % 8000) as usize;
            let truth = circ[st..st + 3000].to_vec();
            let mut r = noisy(&truth, 0.10, 100 + i);
            let mut t = truth;
            if i % 2 == 1 {
                r = revcomp(&r);
                t = revcomp(&t);
            }
            reads.push(Rec {
                id: format!("r{i}").into_bytes(),
                seq: r,
                qual: None,
            });
            truths.push(t);
        }
        let identity = |pairs: &[(Vec<u8>, Vec<u8>)]| {
            let (err, len) = pairs.iter().fold((0usize, 0usize), |(e, l), (q, t)| {
                (e + semi_global(q, t), l + q.len())
            });
            1.0 - err as f64 / len as f64
        };
        let before: Vec<(Vec<u8>, Vec<u8>)> = reads
            .iter()
            .zip(&truths)
            .map(|(r, t)| (r.seq.clone(), t.clone()))
            .collect();
        let p = CorrectParams {
            k: 15,
            min_count: None,
            length_tolerance: 0.15,
            max_gap: 1500,
            max_steps: 20_000,
            keep_ends: false,
            relative: 0.0,
            peak_fraction: 0.3,
        };
        let (out, rep) = round(reads, None, &p).unwrap();
        assert!(rep.stats.gaps_bridged as f64 >= 0.99 * rep.stats.gaps as f64);
        let truth_of: FxHashMap<Vec<u8>, Vec<u8>> = (0..250)
            .map(|i| (format!("r{i}").into_bytes(), truths[i].clone()))
            .collect();
        let after: Vec<(Vec<u8>, Vec<u8>)> = out
            .iter()
            .map(|r| (r.seq.clone(), truth_of[&r.id].clone()))
            .collect();
        let (b, a) = (identity(&before), identity(&after));
        eprintln!("identity {b:.4} -> {a:.5}");
        assert!(b < 0.92, "reads start noisy: {b:.4}");
        assert!(a > 0.999, "corrected identity {a:.5} (from {b:.4})");
        assert!(
            out.len() >= 245,
            "almost every read survives: {}",
            out.len()
        );
    }
}
