//! Native ONT R10 Simplex error correction (`hifiasm --ont` phased overlap algorithm;
//! Cheng, Qu, Li et al., *Nature* 2026).
//!
//! Why R10 needs a distinct error-correction algorithm from R9 / PacBio CLR:
//! - On legacy ONT R9 (`88-92%` raw identity) or quality-less FASTA, random errors dominate every
//!   21-mer; `correct.rs`'s `k = 17, 25` solid-kmer de Bruijn graph bridges noisy gaps rapidly.
//! - On ONT R10.4.1 Simplex (`98-99.3%` raw identity, calibrated Phred FASTQ), random errors are
//!   sparse and the dominant failure modes are:
//!     1. recurrent context-specific / 5mC & 6mA methylation basecalling errors (strongly
//!        strand-biased and enriched at lower Phred Q),
//!     2. homopolymer run-length jitter, and
//!     3. repeat copy / plastid-IR / MTPT homogenization if corrected on an unphased small-k dBG.
//!
//! Following `hifiasm (--ont)`, `run_r10` executes four stages over 2 iterative rounds:
//! 1. Read-level QC (`Q_eff >= 10`), homopolymer-compressed (`HPC`) `(k=15, w=5)` minimizer
//!    indexing, internal/tail chimera clipping at unbridged junctions, and bin-balanced dovetail
//!    overlap selection.
//! 2. Anchor-guided banded Needleman-Wunsch pileups in target-oriented HPC coordinates, with
//!    3-way candidate informative site pre-filtering:
//!    (a) base quality threshold (`Q >= 10`, minor `Q >= 14`),
//!    (b) homopolymer-length invariance in HPC coordinates + long-HP boundary guard, and
//!    (c) symmetric directional strand-bias pruning to strip methylation/systematic errors.
//! 3. $O(n^2)$-time, $O(n)$-space dynamic programming over pairwise four-gamete compatibility to
//!    extract co-segregating haplotype / repeat-copy chains while rejecting unlinked recurrent
//!    errors, inferring each target read's true haplotype across its chain even when the target
//!    read itself carries a single-site error.
//! 4. In-phase peer selection + strand-normalized Phred-weighted HPC base, indel, and homopolymer
//!    run-length consensus.

use std::time::Instant;

use anyhow::Result;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;

use crate::correct::{valley_and_peak, Bloom, ReadStats, RoundReport};
use crate::recruit::Rec;

/// Error-correction chemistry policy for `ovasm run ont` and `ovasm assemble correct`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, clap::ValueEnum, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum OntChemistry {
    /// Automatically use R10 phased overlap correction when FASTQ reads carry calibrated R10-grade
    /// Phred quality scores and high solid 21-mer fraction; keep legacy solid-kmer dBG otherwise.
    #[default]
    Auto,
    /// Force `hifiasm (--ont)`-style phased overlap error correction (ONT R10.4 / Simplex).
    #[value(alias = "r10.4", alias = "r10.4.1", alias = "simplex")]
    R10,
    /// Force the existing `k = 17, 25` solid-kmer de Bruijn graph self-correction (ONT R9 / CLR).
    #[value(alias = "r9", alias = "clr")]
    Legacy,
}

impl OntChemistry {
    pub fn resolve(chemistry: Self, r10_flag: bool) -> Self {
        if r10_flag {
            Self::R10
        } else {
            chemistry
        }
    }

    pub fn as_arg(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::R10 => "r10",
            Self::Legacy => "legacy",
        }
    }
}

const DETECT_K: usize = 21;
const AUTO_MIN_READS: usize = 6;
const AUTO_MIN_QUAL_FRAC: f64 = 0.90;
const AUTO_MIN_EFFECTIVE_Q: f64 = 15.0;
const AUTO_MIN_SOLID_PEAK: u32 = 5;
const AUTO_MIN_SOLID_FRAC: f64 = 0.45;
const AUTO_HIFI_MIN_EFFECTIVE_Q: f64 = 25.0;
const AUTO_HIFI_MIN_SOLID_FRAC: f64 = 0.90;

const K_HPC: usize = 15;
const W_HPC: usize = 5;
const MIN_RAW_LEN: usize = 200;
const MIN_READ_Q_EFF: f64 = 10.0;
const DEFAULT_SYNTH_Q: u8 = 20;

const MIN_SEED_HITS: usize = 6;
const MIN_CHAIN_ANCHORS: usize = 6;
const MAX_CHAIN_GAP: usize = 450;
const MIN_OVLP_HPC: usize = 300;
const MIN_CHIMERA_SPAN_HPC: usize = 360;
const MAX_OVERHANG_HPC: usize = 95;
const BRIDGE_FLANK_HPC: usize = 150;
const END_GRACE_HPC: usize = 32;

const BIN_WIDTH_HPC: usize = 800;
const BIN_PEER_CAP: usize = 32;
const MAX_PEERS: usize = 60;

const MIN_SITE_Q: u8 = 10;
const MIN_MINOR_COUNT: u32 = 3;
const MIN_MINOR_FREQ: f64 = 0.10;
const MAX_CANDIDATE_SITES: usize = 256;
const MIN_PHASE_SEP_HPC: usize = 8;

const DEL_BASE: u8 = 4;
const UNCOV_BASE: u8 = 255;
const BASES: [u8; 4] = *b"ACGT";

#[derive(Clone, Debug, Serialize)]
pub struct ChemistryProfile {
    pub requested: OntChemistry,
    pub selected: &'static str,
    pub use_r10: bool,
    pub reads_examined: usize,
    pub bases_examined: u64,
    pub qual_base_fraction: f64,
    pub effective_phred_q: Option<f64>,
    pub solid_21mer_peak: u32,
    pub solid_21mer_fraction: f64,
    pub reason: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct R10RoundReport {
    pub round: usize,
    pub reads_in: u64,
    pub reads_qc_kept: u64,
    pub reads_chimera_clipped: u64,
    pub bases_chimera_trimmed: u64,
    pub reads_out: u64,
    pub bases_in: u64,
    pub bases_out: u64,
    pub mean_peers: f64,
    pub candidate_sites: u64,
    pub pruned_strand_bias_sites: u64,
    pub pruned_incompatible_sites: u64,
    pub phased_sites_kept: u64,
    pub base_edits: u64,
    pub hp_length_edits: u64,
    pub seconds: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct R10Report {
    pub profile: ChemistryProfile,
    pub rounds: Vec<R10RoundReport>,
}

#[inline]
fn base_code(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

#[inline]
fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[inline]
fn rc_kmer_2bit(mut v: u64, k: usize) -> u64 {
    let mut r = 0u64;
    for _ in 0..k {
        r = (r << 2) | (3 - (v & 3));
        v >>= 2;
    }
    r
}

/// Precomputed Phred error probabilities $10^{-Q/10}$ for $Q \in 0..=60$.
fn phred_err_lut() -> [f64; 61] {
    let mut lut = [0.0f64; 61];
    for (q, slot) in lut.iter_mut().enumerate() {
        *slot = 10.0f64.powf(-(q as f64) / 10.0);
    }
    lut
}

#[inline]
fn decode_phred(ascii_q: u8) -> u8 {
    ascii_q.saturating_sub(33).clamp(2, 50)
}

fn effective_q_of_slice(qual: &[u8], lut: &[f64; 61]) -> f64 {
    if qual.is_empty() {
        return 0.0;
    }
    let sum_err: f64 = qual
        .iter()
        .map(|&q| lut[decode_phred(q) as usize])
        .sum();
    let mean_err = (sum_err / qual.len() as f64).clamp(1e-6, 1.0);
    -10.0 * mean_err.log10()
}

fn forward_21mers(seq: &[u8]) -> Vec<u64> {
    const MASK: u64 = (1u64 << (2 * DETECT_K)) - 1;
    let mut out = Vec::with_capacity(seq.len());
    let (mut v, mut valid) = (0u64, 0usize);
    for &b in seq {
        match base_code(b) {
            Some(c) => {
                v = ((v << 2) | u64::from(c)) & MASK;
                valid += 1;
                if valid >= DETECT_K {
                    let rc = rc_kmer_2bit(v, DETECT_K);
                    out.push(v.min(rc));
                }
            }
            None => valid = 0,
        }
    }
    out
}

fn looks_like_ccs_header(id: &[u8]) -> bool {
    // PacBio CCS/HiFi headers (`m64011_.../123/ccs`, `ccs` token, or SAM tags `rq:f:` / `np:i:`).
    let lower: Vec<u8> = id.iter().take(256).map(u8::to_ascii_lowercase).collect();
    for i in 0..lower.len() {
        if lower[i..].starts_with(b"/ccs")
            && matches!(lower.get(i + 4), None | Some(b' ' | b'\t' | b'/'))
        {
            return true;
        }
        if (lower[i..].starts_with(b"rq:f:") || lower[i..].starts_with(b"np:i:")
            || lower[i..].starts_with(b"rq=") || lower[i..].starts_with(b"np="))
            && (i == 0 || matches!(lower[i - 1], b' ' | b'\t'))
        {
            return true;
        }
    }
    false
}

/// Inspect read Phred calibration and 21-mer solidity to decide between R10 phased overlap
/// correction and the legacy solid-kmer dBG corrector.
pub fn detect_chemistry(reads: &[Rec], mode: OntChemistry) -> ChemistryProfile {
    let lut = phred_err_lut();
    let mut bases_examined = 0u64;
    let mut qual_bases = 0u64;
    let mut err_sum = 0.0f64;
    let mut ccs_headers = 0usize;

    for r in reads {
        bases_examined += r.seq.len() as u64;
        if looks_like_ccs_header(&r.id) {
            ccs_headers += 1;
        }
        if let Some(q) = r.qual.as_ref().filter(|q| q.len() == r.seq.len()) {
            qual_bases += q.len() as u64;
            err_sum += q
                .iter()
                .map(|&b| lut[decode_phred(b) as usize])
                .sum::<f64>();
        }
    }
    let ccs_frac = if !reads.is_empty() {
        ccs_headers as f64 / reads.len() as f64
    } else {
        0.0
    };
    let qual_base_fraction = if bases_examined > 0 {
        qual_bases as f64 / bases_examined as f64
    } else {
        0.0
    };
    let effective_phred_q = (qual_bases > 0).then(|| {
        let mean_err = (err_sum / qual_bases as f64).clamp(1e-6, 1.0);
        -10.0 * mean_err.log10()
    });

    // Fast canonical 21-mer solidity check (skip full table if explicitly forced to Legacy).
    let (solid_21mer_peak, solid_21mer_fraction) = if mode == OntChemistry::Legacy || reads.is_empty()
    {
        (0, 0.0)
    } else {
        let mut bloom = Bloom::new(28);
        let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
        for chunk in reads.chunks(2048) {
            let kms: Vec<Vec<u64>> = chunk.par_iter().map(|r| forward_21mers(&r.seq)).collect();
            for km in kms.into_iter().flatten() {
                if let Some(c) = counts.get_mut(&km) {
                    *c += 1;
                } else if bloom.insert(km) {
                    counts.insert(km, 2);
                }
            }
        }
        let (low, peak) = valley_and_peak(&counts);
        let min_solid = low.max(((0.25 * peak as f64) as u32).max(3));
        let (tot_occ, solid_occ) = reads
            .par_iter()
            .map(|r| {
                let kms = forward_21mers(&r.seq);
                let tot = kms.len() as u64;
                let sol = kms
                    .into_iter()
                    .filter(|km| counts.get(km).copied().unwrap_or(0) >= min_solid)
                    .count() as u64;
                (tot, sol)
            })
            .reduce(|| (0, 0), |a, b| (a.0 + b.0, a.1 + b.1));
        let frac = if tot_occ > 0 {
            solid_occ as f64 / tot_occ as f64
        } else {
            0.0
        };
        (peak, frac)
    };

    let (use_r10, reason) = match mode {
        OntChemistry::R10 => (true, "forced by --chemistry r10 / --r10".to_string()),
        OntChemistry::Legacy => (false, "forced by --chemistry legacy".to_string()),
        OntChemistry::Auto => {
            if reads.len() < AUTO_MIN_READS {
                (
                    false,
                    format!("auto -> legacy: too few reads ({} < {AUTO_MIN_READS})", reads.len()),
                )
            } else if qual_base_fraction < AUTO_MIN_QUAL_FRAC {
                (
                    false,
                    format!(
                        "auto -> legacy: input lacks per-base FASTQ qualities ({:.1}% < {:.0}%)",
                        100.0 * qual_base_fraction,
                        100.0 * AUTO_MIN_QUAL_FRAC
                    ),
                )
            } else if effective_phred_q.unwrap_or(0.0) < AUTO_MIN_EFFECTIVE_Q {
                (
                    false,
                    format!(
                        "auto -> legacy: effective Phred Q {:.2} < {:.1} (R9/CLR noise profile)",
                        effective_phred_q.unwrap_or(0.0),
                        AUTO_MIN_EFFECTIVE_Q
                    ),
                )
            } else if solid_21mer_peak < AUTO_MIN_SOLID_PEAK {
                (
                    false,
                    format!(
                        "auto -> legacy: 21-mer solid peak {} < {AUTO_MIN_SOLID_PEAK}",
                        solid_21mer_peak
                    ),
                )
            } else if solid_21mer_fraction < AUTO_MIN_SOLID_FRAC {
                (
                    false,
                    format!(
                        "auto -> legacy: solid 21-mer fraction {:.3} < {:.2}",
                        solid_21mer_fraction, AUTO_MIN_SOLID_FRAC
                    ),
                )
            } else if ccs_frac >= 0.50
                || (effective_phred_q.unwrap_or(0.0) >= AUTO_HIFI_MIN_EFFECTIVE_Q
                    && solid_21mer_fraction >= AUTO_HIFI_MIN_SOLID_FRAC)
            {
                let msg = format!(
                    "reads look like PacBio HiFi (effective Q {:.2}, solid 21-mer fraction {:.3}{}); keeping legacy path — use `ovasm run hifi` for HiFi or `--chemistry r10` to force R10",
                    effective_phred_q.unwrap_or(0.0),
                    solid_21mer_fraction,
                    if ccs_frac >= 0.50 { ", CCS headers" } else { "" }
                );
                eprintln!("[ovasm] warning: {msg}");
                (false, format!("auto -> legacy: {msg}"))
            } else {
                (
                    true,
                    format!(
                        "auto -> r10: effective Phred Q {:.2}, solid 21-mer fraction {:.3} (peak {})",
                        effective_phred_q.unwrap_or(0.0),
                        solid_21mer_fraction,
                        solid_21mer_peak
                    ),
                )
            }
        }
    };

    ChemistryProfile {
        requested: mode,
        selected: if use_r10 { "r10" } else { "legacy" },
        use_r10,
        reads_examined: reads.len(),
        bases_examined,
        qual_base_fraction,
        effective_phred_q,
        solid_21mer_peak,
        solid_21mer_fraction,
        reason,
    }
}

/// Homopolymer-compressed representation of one read.
#[derive(Clone)]
struct HpcRead {
    id: Vec<u8>,
    /// `0..4` for `A,C,G,T`; `4` for ambiguous `N` boundary.
    base: Vec<u8>,
    /// Homopolymer run length (`>= 1`).
    rlen: Vec<u16>,
    /// Maximum base Phred Q in the run (confidence that this run's base identity exists).
    q_max: Vec<u8>,
    /// Effective Phred Q across the run (confidence in the run's length).
    q_mean: Vec<u8>,
    raw_len: usize,
}

impl HpcRead {
    fn hpc_len(&self) -> usize {
        self.base.len()
    }

    fn slice(&self, start: usize, end: usize) -> Self {
        let raw_len: usize = self.rlen[start..end].iter().map(|&l| l as usize).sum();
        Self {
            id: self.id.clone(),
            base: self.base[start..end].to_vec(),
            rlen: self.rlen[start..end].to_vec(),
            q_max: self.q_max[start..end].to_vec(),
            q_mean: self.q_mean[start..end].to_vec(),
            raw_len,
        }
    }

    fn expand_fasta(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.raw_len);
        for (&b, &len) in self.base.iter().zip(&self.rlen) {
            let ch = if (b as usize) < 4 { BASES[b as usize] } else { b'N' };
            out.resize(out.len() + len as usize, ch);
        }
        out
    }
}

fn encode_hpc_read(
    id: Vec<u8>,
    seq: &[u8],
    qual_phred: &[u8],
    lut: &[f64; 61],
) -> Option<HpcRead> {
    let acgt = seq.iter().filter(|&&b| base_code(b).is_some()).count();
    if acgt < MIN_RAW_LEN {
        return None;
    }
    let mut base: Vec<u8> = Vec::with_capacity(seq.len() * 4 / 5 + 8);
    let mut rlen: Vec<u16> = Vec::with_capacity(base.capacity());
    let mut q_max: Vec<u8> = Vec::with_capacity(base.capacity());
    let mut q_mean: Vec<u8> = Vec::with_capacity(base.capacity());

    let mut i = 0usize;
    while i < seq.len() {
        if let Some(code) = base_code(seq[i]) {
            let mut j = i + 1;
            while j < seq.len() && base_code(seq[j]) == Some(code) && (j - i) < u16::MAX as usize {
                j += 1;
            }
            let len = (j - i) as u16;
            let mut mx = 0u8;
            let mut err_sum = 0.0f64;
            for &qp in &qual_phred[i..j] {
                let q = qp.clamp(2, 50);
                mx = mx.max(q);
                err_sum += lut[q as usize];
            }
            let mean_err = (err_sum / f64::from(len)).clamp(1e-5, 1.0);
            let qm = (-10.0 * mean_err.log10()).round().clamp(2.0, 50.0) as u8;
            base.push(code);
            rlen.push(len);
            q_max.push(mx);
            q_mean.push(qm);
            i = j;
        } else {
            base.push(4);
            rlen.push(1);
            q_max.push(2);
            q_mean.push(2);
            i += 1;
        }
    }
    if base.len() < K_HPC * 2 {
        return None;
    }
    Some(HpcRead {
        id,
        base,
        rlen,
        q_max,
        q_mean,
        raw_len: seq.len(),
    })
}

/// Extract `(pos, canonical_hash, canonical_strand)` `(k=15, w=5)` minimizers from HPC bases.
fn hpc_minimizers(hpc_base: &[u8]) -> Vec<(u32, u64, u8)> {
    const MASK: u64 = (1u64 << (2 * K_HPC)) - 1;
    if hpc_base.len() < K_HPC {
        return Vec::new();
    }
    let mut kmers: Vec<(u32, u64, u8)> = Vec::with_capacity(hpc_base.len());
    let (mut fwd, mut rev, mut valid) = (0u64, 0u64, 0usize);
    let rc_shift = 2 * (K_HPC - 1);
    for (i, &b) in hpc_base.iter().enumerate() {
        if b < 4 {
            let c = u64::from(b);
            fwd = ((fwd << 2) | c) & MASK;
            rev = (rev >> 2) | ((3 - c) << rc_shift);
            valid += 1;
            if valid >= K_HPC {
                let pos = (i + 1 - K_HPC) as u32;
                // K_HPC = 15 is odd, so `fwd != rev` always holds.
                let (canon, strand) = if fwd < rev { (fwd, 0u8) } else { (rev, 1u8) };
                kmers.push((pos, mix64(canon), strand));
            }
        } else {
            valid = 0;
            fwd = 0;
            rev = 0;
        }
    }
    if kmers.is_empty() {
        return Vec::new();
    }
    if kmers.len() <= W_HPC {
        let best = kmers.iter().copied().min_by_key(|x| x.1).unwrap();
        return vec![best];
    }
    let mut out: Vec<(u32, u64, u8)> = Vec::with_capacity(kmers.len() * 2 / (W_HPC + 1) + 4);
    let mut last_idx = usize::MAX;
    for w_start in 0..=kmers.len() - W_HPC {
        let mut min_idx = w_start;
        for j in (w_start + 1)..(w_start + W_HPC) {
            if kmers[j].1 < kmers[min_idx].1 {
                min_idx = j;
            }
        }
        if min_idx != last_idx {
            out.push(kmers[min_idx]);
            last_idx = min_idx;
        }
    }
    out
}

#[derive(Clone, Copy)]
struct IndexHit {
    read_id: u32,
    pos_and_strand: u32, // `(pos << 1) | strand`
}

struct HpcIndex {
    /// `hash -> (offset, len)` into `hits`.
    table: FxHashMap<u64, (u32, u32)>,
    hits: Vec<IndexHit>,
    read_mins: Vec<Vec<(u32, u64, u8)>>,
}

impl HpcIndex {
    fn build(reads: &[HpcRead], max_occ: usize) -> Self {
        let read_mins: Vec<Vec<(u32, u64, u8)>> =
            reads.par_iter().map(|r| hpc_minimizers(&r.base)).collect();

        let mut counts: FxHashMap<u64, u32> =
            FxHashMap::with_capacity_and_hasher(reads.len() * 64, Default::default());
        for mins in &read_mins {
            for &(_, h, _) in mins {
                *counts.entry(h).or_insert(0) += 1;
            }
        }

        let mut table: FxHashMap<u64, (u32, u32)> =
            FxHashMap::with_capacity_and_hasher(counts.len(), Default::default());
        let mut total = 0u32;
        for (h, c) in counts {
            if c >= 2 && (c as usize) <= max_occ {
                table.insert(h, (total, 0));
                total = total.saturating_add(c);
            }
        }

        let mut hits = vec![
            IndexHit {
                read_id: 0,
                pos_and_strand: 0,
            };
            total as usize
        ];
        for (rid, mins) in read_mins.iter().enumerate() {
            let rid = rid as u32;
            for &(pos, h, strand) in mins {
                if let Some(slot) = table.get_mut(&h) {
                    let idx = (slot.0 + slot.1) as usize;
                    hits[idx] = IndexHit {
                        read_id: rid,
                        pos_and_strand: (pos << 1) | u32::from(strand),
                    };
                    slot.1 += 1;
                }
            }
        }

        Self {
            table,
            hits,
            read_mins,
        }
    }

    #[inline]
    fn lookup(&self, h: u64) -> &[IndexHit] {
        match self.table.get(&h) {
            Some(&(off, len)) => &self.hits[off as usize..(off + len) as usize],
            None => &[],
        }
    }
}

#[derive(Clone)]
struct OverlapChain {
    peer_id: u32,
    is_rev: bool,
    t_start: u32,
    t_end: u32,
    q_start: u32,
    q_end: u32,
    overhang: u32,
    score: i32,
    /// Anchor `(t_pos, q_oriented_pos)` pairs of exact `K_HPC`-mers in `t`'s forward frame.
    anchors: Vec<(u32, u32)>,
}

fn collect_chains_for_read(
    t_id: usize,
    reads: &[HpcRead],
    index: &HpcIndex,
    min_span: usize,
    require_dovetail: bool,
) -> Vec<OverlapChain> {
    let t_len = reads[t_id].hpc_len();
    let t_mins = &index.read_mins[t_id];
    if t_mins.len() < MIN_SEED_HITS {
        return Vec::new();
    }

    // Pass 1: count shared minimizers and coarse `[t_min, t_max]` / `[q_min, q_max]` bounds per
    // candidate `(peer_id, is_rev)` so 300x-600x organelles do not sort/chain thousands of
    // redundant or non-dovetail peers per read.
    let mut cand_summary: FxHashMap<u32, (u16, u32, u32, u32, u32)> =
        FxHashMap::with_capacity_and_hasher(256, Default::default());
    for &(t_pos, h, t_strand) in t_mins {
        for hit in index.lookup(h) {
            let q_id = hit.read_id as usize;
            if q_id == t_id {
                continue;
            }
            let q_len = reads[q_id].hpc_len() as u32;
            if q_len < (K_HPC * 2) as u32 {
                continue;
            }
            let q_pos_raw = hit.pos_and_strand >> 1;
            let q_strand = (hit.pos_and_strand & 1) as u8;
            let is_rev = t_strand != q_strand;
            let q_pos = if is_rev {
                q_len - (K_HPC as u32) - q_pos_raw
            } else {
                q_pos_raw
            };
            let peer_key = (hit.read_id << 1) | u32::from(is_rev);
            cand_summary
                .entry(peer_key)
                .and_modify(|e| {
                    e.0 = e.0.saturating_add(1);
                    e.1 = e.1.min(t_pos);
                    e.2 = e.2.max(t_pos);
                    e.3 = e.3.min(q_pos);
                    e.4 = e.4.max(q_pos);
                })
                .or_insert((1, t_pos, t_pos, q_pos, q_pos));
        }
    }

    let mut viable: Vec<(u32, u16, u32, u32)> = Vec::with_capacity(cand_summary.len().min(256));
    let oh_coarse = (MAX_OVERHANG_HPC + 48) as u32;
    for (peer_key, (cnt, t_lo, t_hi, q_lo, q_hi)) in cand_summary {
        if (cnt as usize) < MIN_SEED_HITS {
            continue;
        }
        let q_len = reads[(peer_key >> 1) as usize].hpc_len() as u32;
        let needed_span = min_span.min(t_len * 3 / 4).min(q_len as usize * 3 / 4);
        let t_span_est = (t_hi + K_HPC as u32).saturating_sub(t_lo) as usize;
        if t_span_est < needed_span {
            continue;
        }
        if require_dovetail {
            let left_oh = t_lo.min(q_lo);
            let right_oh =
                (t_len as u32).saturating_sub(t_hi + K_HPC as u32).min(q_len.saturating_sub(q_hi + K_HPC as u32));
            if left_oh > oh_coarse || right_oh > oh_coarse {
                continue;
            }
        }
        viable.push((peer_key, cnt, t_lo, t_hi + K_HPC as u32));
    }
    if viable.is_empty() {
        return Vec::new();
    }

    const MAX_CHAIN_CANDS: usize = 128;
    let keep_peers: FxHashSet<u32> = if viable.len() <= MAX_CHAIN_CANDS {
        viable.into_iter().map(|x| x.0).collect()
    } else {
        viable.sort_unstable_by_key(|x| std::cmp::Reverse(x.1));
        let n_bins = (t_len / BIN_WIDTH_HPC) + 1;
        let mut bin_cov = vec![0u16; n_bins];
        let mut set = FxHashSet::with_capacity_and_hasher(MAX_CHAIN_CANDS, Default::default());
        for (pk, _, ts, te) in viable {
            if set.len() >= MAX_CHAIN_CANDS {
                break;
            }
            let b0 = (ts as usize / BIN_WIDTH_HPC).min(n_bins - 1);
            let b1 = (te.saturating_sub(1) as usize / BIN_WIDTH_HPC).min(n_bins - 1);
            if (b0..=b1).any(|b| bin_cov[b] < 48) {
                for slot in &mut bin_cov[b0..=b1] {
                    *slot = slot.saturating_add(1);
                }
                set.insert(pk);
            }
        }
        set
    };

    // Pass 2: gather `(peer_key, t_pos, q_oriented_pos)` only for viable peers.
    let mut raw_hits: Vec<(u32, u32, u32)> = Vec::with_capacity(keep_peers.len() * 64);
    for &(t_pos, h, t_strand) in t_mins {
        for hit in index.lookup(h) {
            let q_id = hit.read_id as usize;
            if q_id == t_id {
                continue;
            }
            let q_pos_raw = hit.pos_and_strand >> 1;
            let q_strand = (hit.pos_and_strand & 1) as u8;
            let is_rev = t_strand != q_strand;
            let peer_key = (hit.read_id << 1) | u32::from(is_rev);
            if !keep_peers.contains(&peer_key) {
                continue;
            }
            let q_len = reads[q_id].hpc_len() as u32;
            let q_pos = if is_rev {
                q_len - (K_HPC as u32) - q_pos_raw
            } else {
                q_pos_raw
            };
            raw_hits.push((peer_key, t_pos, q_pos));
        }
    }
    if raw_hits.len() < MIN_SEED_HITS {
        return Vec::new();
    }
    raw_hits.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));

    let mut chains = Vec::new();
    let mut dp_score: Vec<i32> = Vec::new();
    let mut dp_prev: Vec<u32> = Vec::new();

    let mut i = 0usize;
    while i < raw_hits.len() {
        let key = raw_hits[i].0;
        let mut j = i + 1;
        while j < raw_hits.len() && raw_hits[j].0 == key {
            j += 1;
        }
        let slice = &raw_hits[i..j];
        i = j;
        if slice.len() < MIN_SEED_HITS {
            continue;
        }

        // Deduplicate identical `(t_pos, q_pos)` pairs.
        let mut pts: Vec<(u32, u32)> = Vec::with_capacity(slice.len());
        for &(_, tp, qp) in slice {
            if pts.last().copied() != Some((tp, qp)) {
                pts.push((tp, qp));
            }
        }
        if pts.len() < MIN_SEED_HITS {
            continue;
        }

        let m = pts.len();
        dp_score.clear();
        dp_score.resize(m, K_HPC as i32);
        dp_prev.clear();
        dp_prev.resize(m, u32::MAX);

        let mut best_end = 0usize;
        for cur in 1..m {
            let (tc, qc) = pts[cur];
            let mut best_sc = K_HPC as i32;
            let mut best_pr = u32::MAX;
            let mut back = cur;
            let mut steps = 0usize;
            while back > 0 && steps < 28 {
                back -= 1;
                let (tp, qp) = pts[back];
                if tc <= tp {
                    continue;
                }
                let dt = (tc - tp) as usize;
                if dt > MAX_CHAIN_GAP {
                    break;
                }
                steps += 1;
                if qc <= qp {
                    continue;
                }
                let dq = (qc - qp) as usize;
                if dq > MAX_CHAIN_GAP {
                    continue;
                }
                let diff = dt.abs_diff(dq);
                let max_skew = (dt.max(dq) / 8) + 24;
                if diff > max_skew {
                    continue;
                }
                let bonus = dt.min(dq).min(K_HPC) as i32;
                let gap_pen = (diff as i32) + ((dt.max(dq) / 32) as i32);
                let cand = dp_score[back] + bonus - gap_pen;
                if cand > best_sc {
                    best_sc = cand;
                    best_pr = back as u32;
                }
            }
            dp_score[cur] = best_sc;
            dp_prev[cur] = best_pr;
            if best_sc > dp_score[best_end] {
                best_end = cur;
            }
        }

        let mut anchors: Vec<(u32, u32)> = Vec::new();
        let mut cur = best_end as u32;
        while cur != u32::MAX {
            anchors.push(pts[cur as usize]);
            cur = dp_prev[cur as usize];
        }
        anchors.reverse();
        if anchors.len() < MIN_CHAIN_ANCHORS {
            continue;
        }

        let peer_id = key >> 1;
        let is_rev = (key & 1) != 0;
        let q_len = reads[peer_id as usize].hpc_len() as u32;
        let t_start = anchors[0].0;
        let q_start = anchors[0].1;
        let t_end = anchors.last().unwrap().0 + K_HPC as u32;
        let q_end = anchors.last().unwrap().1 + K_HPC as u32;
        let t_span = (t_end - t_start) as usize;
        let needed_span = min_span.min(t_len * 3 / 4).min(q_len as usize * 3 / 4);
        if t_span < needed_span {
            continue;
        }
        // Require reasonable anchor density (>= 1 anchor per 35 HPC bases) so noisy non-homologous
        // stretches cannot chain.
        if anchors.len() * 35 < t_span {
            continue;
        }

        let left_oh = t_start.min(q_start);
        let right_oh = (t_len as u32 - t_end).min(q_len - q_end);
        let overhang = left_oh + right_oh;
        if require_dovetail
            && (left_oh as usize > MAX_OVERHANG_HPC || right_oh as usize > MAX_OVERHANG_HPC)
        {
            continue;
        }

        chains.push(OverlapChain {
            peer_id,
            is_rev,
            t_start,
            t_end,
            q_start,
            q_end,
            overhang,
            score: dp_score[best_end],
            anchors,
        });
    }

    chains
}

/// Select up to `MAX_PEERS` dovetail peers with uniform bin depth across ultra-long reads.
fn select_balanced_peers(mut chains: Vec<OverlapChain>, t_len: usize) -> Vec<OverlapChain> {
    if chains.len() <= MAX_PEERS {
        chains.sort_unstable_by(|a, b| b.score.cmp(&a.score).then(a.overhang.cmp(&b.overhang)));
        return chains;
    }
    chains.sort_unstable_by(|a, b| b.score.cmp(&a.score).then(a.overhang.cmp(&b.overhang)));
    let n_bins = (t_len / BIN_WIDTH_HPC) + 1;
    let mut bin_cov = vec![0u16; n_bins];
    let mut selected = Vec::with_capacity(MAX_PEERS);

    for ch in chains {
        if selected.len() >= MAX_PEERS {
            break;
        }
        let b0 = (ch.t_start as usize / BIN_WIDTH_HPC).min(n_bins - 1);
        let b1 = (ch.t_end.saturating_sub(1) as usize / BIN_WIDTH_HPC).min(n_bins - 1);
        let needs_cov = (b0..=b1).any(|b| (bin_cov[b] as usize) < BIN_PEER_CAP);
        if !needs_cov {
            continue;
        }
        for slot in &mut bin_cov[b0..=b1] {
            *slot = slot.saturating_add(1);
        }
        selected.push(ch);
    }
    selected
}

/// Stage 1A: clip unsupported chimeric tails and internal unbridged chimeric joins without eroding
/// genuine read ends.
fn clip_chimeras(reads: Vec<HpcRead>, max_occ: usize) -> (Vec<HpcRead>, u64, u64) {
    let index = HpcIndex::build(&reads, max_occ);
    let results: Vec<(Option<HpcRead>, bool, u64)> = (0..reads.len())
        .into_par_iter()
        .map(|t_id| {
            let r = &reads[t_id];
            let t_len = r.hpc_len();
            let chains =
                collect_chains_for_read(t_id, &reads, &index, MIN_CHIMERA_SPAN_HPC, false);
            if chains.len() < 2 {
                // Keep sparse reads intact rather than over-clipping low-depth regions.
                return (Some(r.clone()), false, 0);
            }
            let mut cov_diff = vec![0i32; t_len + 1];
            let mut br_diff = vec![0i32; t_len + 1];
            for ch in &chains {
                let q_len = reads[ch.peer_id as usize].hpc_len();
                let mut s = ch.t_start as usize;
                let mut e = (ch.t_end as usize).min(t_len);
                // Minimizer windows `(k=15, w=5)` stop up to ~19 HPC bases inside true read ends;
                // extend to `0` / `t_len` when the alignment reaches the end of `t` (or `q` overhangs `t`).
                if s <= END_GRACE_HPC && (ch.q_start >= ch.t_start || ch.q_start as usize <= END_GRACE_HPC) {
                    s = 0;
                }
                if e + END_GRACE_HPC >= t_len
                    && ((q_len.saturating_sub(ch.q_end as usize)) >= (t_len - e)
                        || (ch.q_end as usize) + END_GRACE_HPC >= q_len)
                {
                    e = t_len;
                }
                if s < e {
                    cov_diff[s] += 1;
                    cov_diff[e] -= 1;
                }
                let bs = s + BRIDGE_FLANK_HPC;
                let be = e.saturating_sub(BRIDGE_FLANK_HPC);
                if bs < be {
                    br_diff[bs] += 1;
                    br_diff[be] -= 1;
                }
            }
            let mut cov = vec![0i32; t_len];
            let mut br = vec![0i32; t_len];
            let (mut c_acc, mut b_acc) = (0i32, 0i32);
            for p in 0..t_len {
                c_acc += cov_diff[p];
                b_acc += br_diff[p];
                cov[p] = c_acc;
                br[p] = b_acc;
            }

            // Identify maximal segments with `cov >= 2` and no internal unbridged break.
            // Internal breaks are tested at least `BRIDGE_FLANK_HPC + END_GRACE_HPC + 32` inside
            // `t` so read-end minimizer margins are never mistaken for chimeras.
            let internal_margin = BRIDGE_FLANK_HPC + END_GRACE_HPC + 32;
            let mut best_seg: Option<(usize, usize)> = None;
            let mut p = 0usize;
            while p < t_len {
                if cov[p] < 2 {
                    p += 1;
                    continue;
                }
                let seg_start = p;
                p += 1;
                while p < t_len && cov[p] >= 2 {
                    let is_break = if p >= internal_margin && p + internal_margin < t_len {
                        let lc = cov[p - 90];
                        let rc = cov[p + 90];
                        (lc >= 6 && rc >= 6 && br[p] <= 1) || (lc >= 3 && rc >= 3 && br[p] == 0)
                    } else {
                        false
                    };
                    if is_break {
                        break;
                    }
                    p += 1;
                }
                let seg_end = p;
                if best_seg.is_none_or(|(bs, be)| seg_end - seg_start > be - bs) {
                    best_seg = Some((seg_start, seg_end));
                }
            }

            let Some((mut s, mut e)) = best_seg else {
                return (None, true, r.raw_len as u64);
            };
            if s <= END_GRACE_HPC {
                s = 0;
            }
            if e + END_GRACE_HPC >= t_len {
                e = t_len;
            }
            if s == 0 && e == t_len {
                return (Some(r.clone()), false, 0);
            }
            let sliced = r.slice(s, e);
            let trimmed = r.raw_len.saturating_sub(sliced.raw_len) as u64;
            if sliced.raw_len < MIN_RAW_LEN || sliced.hpc_len() < K_HPC * 2 {
                (None, true, r.raw_len as u64)
            } else {
                (Some(sliced), true, trimmed)
            }
        })
        .collect();

    let mut out = Vec::with_capacity(reads.len());
    let mut clipped_reads = 0u64;
    let mut trimmed_bases = 0u64;
    for (opt_r, clipped, trim) in results {
        if clipped {
            clipped_reads += 1;
            trimmed_bases += trim;
        }
        if let Some(r) = opt_r {
            out.push(r);
        }
    }
    (out, clipped_reads, trimmed_bases)
}

#[derive(Clone, Copy)]
struct InsRun {
    base: u8,
    rlen: u16,
    q_base: u8,
    q_len: u8,
}

const MAX_POCKET_INS: usize = 4;

#[derive(Clone, Copy)]
struct SmallIns {
    len: u8,
    runs: [InsRun; MAX_POCKET_INS],
}

impl SmallIns {
    fn empty() -> Self {
        Self {
            len: 0,
            runs: [InsRun {
                base: 0,
                rlen: 1,
                q_base: 0,
                q_len: 0,
            }; MAX_POCKET_INS],
        }
    }

    fn push(&mut self, run: InsRun) {
        if (self.len as usize) < MAX_POCKET_INS {
            self.runs[self.len as usize] = run;
            self.len += 1;
        }
    }

    fn as_slice(&self) -> &[InsRun] {
        &self.runs[..self.len as usize]
    }
}

/// Compact aligned view of one peer onto target read `t`.
struct PeerAligned {
    is_rev: bool,
    t_start: u32,
    t_end: u32,
    col_base: Vec<u8>,
    col_rlen: Vec<u16>,
    col_q_base: Vec<u8>,
    col_q_len: Vec<u8>,
    /// Sorted `(t_pos, SmallIns)` for insertions immediately after `t_pos` (`u32::MAX` = before 0).
    ins_after: Vec<(u32, SmallIns)>,
}

impl PeerAligned {
    #[inline]
    fn covers(&self, t_pos: usize) -> bool {
        (t_pos as u32) >= self.t_start && (t_pos as u32) < self.t_end
    }

    #[inline]
    fn spans_after(&self, t_pos: usize) -> bool {
        (t_pos as u32) >= self.t_start && ((t_pos + 1) as u32) < self.t_end
    }

    #[inline]
    fn call_at(&self, t_pos: usize) -> Option<(u8, u16, u8, u8)> {
        if !self.covers(t_pos) {
            return None;
        }
        let idx = t_pos - self.t_start as usize;
        let b = self.col_base[idx];
        if b == UNCOV_BASE {
            None
        } else {
            Some((
                b,
                self.col_rlen[idx],
                self.col_q_base[idx],
                self.col_q_len[idx],
            ))
        }
    }

    #[inline]
    fn ins_at(&self, t_pos: u32) -> Option<&SmallIns> {
        self.ins_after
            .binary_search_by_key(&t_pos, |x| x.0)
            .ok()
            .map(|idx| &self.ins_after[idx].1)
    }
}

struct OrientedPeer<'a> {
    src: &'a HpcRead,
    is_rev: bool,
    len: usize,
}

impl<'a> OrientedPeer<'a> {
    fn new(src: &'a HpcRead, is_rev: bool) -> Self {
        Self {
            src,
            is_rev,
            len: src.hpc_len(),
        }
    }

    #[inline]
    fn base(&self, q_pos: usize) -> u8 {
        if self.is_rev {
            let b = self.src.base[self.len - 1 - q_pos];
            if b < 4 { 3 - b } else { 4 }
        } else {
            self.src.base[q_pos]
        }
    }

    #[inline]
    fn rlen(&self, q_pos: usize) -> u16 {
        if self.is_rev {
            self.src.rlen[self.len - 1 - q_pos]
        } else {
            self.src.rlen[q_pos]
        }
    }

    #[inline]
    fn q_max(&self, q_pos: usize) -> u8 {
        if self.is_rev {
            self.src.q_max[self.len - 1 - q_pos]
        } else {
            self.src.q_max[q_pos]
        }
    }

    #[inline]
    fn q_mean(&self, q_pos: usize) -> u8 {
        if self.is_rev {
            self.src.q_mean[self.len - 1 - q_pos]
        } else {
            self.src.q_mean[q_pos]
        }
    }
}

/// Banded global Needleman-Wunsch on a short inter-anchor HPC pocket `(t[ts..te], q[qs..qe])`.
/// Because `q` is already oriented into `t`'s strand, identical sequences on `+` and `-` peers
/// receive identical tracebacks (zero strand gap-placement skew).
#[allow(clippy::too_many_arguments)]
fn align_pocket(
    t: &HpcRead,
    ts: usize,
    te: usize,
    q: &OrientedPeer<'_>,
    qs: usize,
    qe: usize,
    out_off: usize,
    col_base: &mut [u8],
    col_rlen: &mut [u16],
    col_q_base: &mut [u8],
    col_q_len: &mut [u8],
    ins_after: &mut Vec<(u32, SmallIns)>,
) {
    let n = te - ts;
    let m = qe - qs;
    if n == 0 && m == 0 {
        return;
    }
    if n == 0 {
        if m <= MAX_POCKET_INS && ts > 0 {
            let mut si = SmallIns::empty();
            for j in 0..m {
                let qp = qs + j;
                let qb = q.base(qp);
                if qb < 4 {
                    si.push(InsRun {
                        base: qb,
                        rlen: q.rlen(qp),
                        q_base: q.q_max(qp),
                        q_len: q.q_mean(qp),
                    });
                }
            }
            if si.len > 0 {
                ins_after.push(((ts - 1) as u32, si));
            }
        }
        return;
    }
    if m == 0 {
        let del_q = if qs > 0 && qs < q.len {
            q.q_max(qs - 1).min(q.q_max(qs))
        } else if qs > 0 {
            q.q_max(qs - 1)
        } else {
            q.q_max(0)
        };
        for i in 0..n {
            let idx = (ts + i) - out_off;
            col_base[idx] = DEL_BASE;
            col_rlen[idx] = 0;
            col_q_base[idx] = del_q;
            col_q_len[idx] = del_q;
        }
        return;
    }

    const NEG_INF: i32 = i32::MIN / 4;
    const GAP_PEN: i32 = 4;
    let band = (n.abs_diff(m) + 14).max(18).min(n.max(m) + 1);
    let width = m + 1;
    let mut dp = vec![NEG_INF; (n + 1) * width];
    // 0 = diag, 1 = up (DEL on t), 2 = left (INS from q)
    let mut bt = vec![0u8; (n + 1) * width];

    dp[0] = 0;
    for j in 1..=m.min(band) {
        dp[j] = -(j as i32) * GAP_PEN;
        bt[j] = 2;
    }
    for i in 1..=n {
        let center = (i * m) / n.max(1);
        let j_lo = center.saturating_sub(band);
        let j_hi = (center + band).min(m);
        if j_lo == 0 {
            dp[i * width] = -(i as i32) * GAP_PEN;
            bt[i * width] = 1;
        }
        let tb = t.base[ts + i - 1];
        let tl = t.rlen[ts + i - 1];
        for j in j_lo.max(1)..=j_hi {
            let qb = q.base(qs + j - 1);
            let ql = q.rlen(qs + j - 1);
            let match_sc = if tb == qb && tb < 4 {
                let dlen = tl.abs_diff(ql);
                6 + if dlen <= 1 { 2 } else if dlen <= 2 { 1 } else { 0 }
            } else {
                -5
            };
            let diag = dp[(i - 1) * width + (j - 1)].saturating_add(match_sc);
            let up = dp[(i - 1) * width + j].saturating_sub(GAP_PEN);
            let left = dp[i * width + (j - 1)].saturating_sub(GAP_PEN);

            // Deterministic tie-break: prefer earlier match (`up`/`left` whose match already
            // occurred upstream, or `diag` over starting a new gap).
            let (best, dir) = if diag >= up && diag >= left {
                (diag, 0u8)
            } else if up >= left {
                (up, 1u8)
            } else {
                (left, 2u8)
            };
            dp[i * width + j] = best;
            bt[i * width + j] = dir;
        }
    }

    if dp[n * width + m] <= NEG_INF / 2 {
        return;
    }

    // Traceback from `(n, m)` to `(0, 0)`.
    let mut ops: Vec<(u8, usize, usize)> = Vec::with_capacity(n + m);
    let (mut i, mut j) = (n, m);
    while i > 0 || j > 0 {
        let dir = if i == 0 {
            2
        } else if j == 0 {
            1
        } else {
            bt[i * width + j]
        };
        match dir {
            0 => {
                i -= 1;
                j -= 1;
                ops.push((0, i, j));
            }
            1 => {
                i -= 1;
                ops.push((1, i, j));
            }
            _ => {
                j -= 1;
                ops.push((2, i, j));
            }
        }
    }
    ops.reverse();

    let mut cur_ins_pos: Option<u32> = None;
    let mut cur_ins = SmallIns::empty();
    let flush_ins = |pos: &mut Option<u32>, ins: &mut SmallIns, out: &mut Vec<(u32, SmallIns)>| {
        if let Some(p) = pos.take() {
            if ins.len > 0 {
                out.push((p, *ins));
            }
            *ins = SmallIns::empty();
        }
    };

    for (op, ti, qj) in ops {
        match op {
            0 => {
                flush_ins(&mut cur_ins_pos, &mut cur_ins, ins_after);
                let t_abs = ts + ti;
                let q_abs = qs + qj;
                let qb = q.base(q_abs);
                if qb < 4 {
                    let idx = t_abs - out_off;
                    col_base[idx] = qb;
                    col_rlen[idx] = q.rlen(q_abs);
                    col_q_base[idx] = q.q_max(q_abs);
                    col_q_len[idx] = q.q_mean(q_abs);
                }
            }
            1 => {
                flush_ins(&mut cur_ins_pos, &mut cur_ins, ins_after);
                let t_abs = ts + ti;
                let q_abs = (qs + qj).min(q.len.saturating_sub(1));
                let del_q = if qs + qj > 0 && qs + qj < q.len {
                    q.q_max(qs + qj - 1).min(q.q_max(qs + qj))
                } else {
                    q.q_max(q_abs)
                };
                let idx = t_abs - out_off;
                col_base[idx] = DEL_BASE;
                col_rlen[idx] = 0;
                col_q_base[idx] = del_q;
                col_q_len[idx] = del_q;
            }
            _ => {
                // Insertion right after `ts + ti - 1` (skip leading insertions before `0`).
                if ts + ti > 0 {
                    let p = (ts + ti - 1) as u32;
                    if cur_ins_pos != Some(p) {
                        flush_ins(&mut cur_ins_pos, &mut cur_ins, ins_after);
                        cur_ins_pos = Some(p);
                    }
                    let q_abs = qs + qj;
                    let qb = q.base(q_abs);
                    if qb < 4 {
                        cur_ins.push(InsRun {
                            base: qb,
                            rlen: q.rlen(q_abs),
                            q_base: q.q_max(q_abs),
                            q_len: q.q_mean(q_abs),
                        });
                    }
                }
            }
        }
    }
    flush_ins(&mut cur_ins_pos, &mut cur_ins, ins_after);
}

/// Align one selected dovetail chain onto `t` using exact `K_HPC`-mer anchor blocks + banded NW
/// on inter-anchor pockets and small dovetail ends.
fn align_peer_chain(t: &HpcRead, peer: &HpcRead, ch: &OverlapChain) -> PeerAligned {
    let q = OrientedPeer::new(peer, ch.is_rev);
    let t_len = t.hpc_len();

    // Merge overlapping exact `K_HPC`-mer anchors into contiguous exact blocks `(ts, qs, len)`.
    let mut blocks: Vec<(usize, usize, usize)> = Vec::with_capacity(ch.anchors.len());
    for &(tp, qp) in &ch.anchors {
        let (tp, qp) = (tp as usize, qp as usize);
        if let Some(last) = blocks.last_mut() {
            let dt = tp.saturating_sub(last.0);
            let dq = qp.saturating_sub(last.1);
            if dt == dq && tp <= last.0 + last.2 {
                last.2 = last.2.max((tp - last.0) + K_HPC);
                continue;
            }
        }
        blocks.push((tp, qp, K_HPC));
    }

    let first = blocks[0];
    let last = *blocks.last().unwrap();
    let ext_left = first.0.min(first.1).min(MAX_OVERHANG_HPC);
    let ext_right = (t_len - (last.0 + last.2))
        .min(q.len - (last.1 + last.2))
        .min(MAX_OVERHANG_HPC);

    let cov_start = first.0 - ext_left;
    let cov_end = last.0 + last.2 + ext_right;
    let span = cov_end - cov_start;

    let mut col_base = vec![UNCOV_BASE; span];
    let mut col_rlen = vec![0u16; span];
    let mut col_q_base = vec![0u8; span];
    let mut col_q_len = vec![0u8; span];
    let mut ins_after: Vec<(u32, SmallIns)> = Vec::new();

    // 1. Optional short left dovetail extension before the first exact block.
    if ext_left > 0 {
        align_pocket(
            t,
            first.0 - ext_left,
            first.0,
            &q,
            first.1 - ext_left,
            first.1,
            cov_start,
            &mut col_base,
            &mut col_rlen,
            &mut col_q_base,
            &mut col_q_len,
            &mut ins_after,
        );
    }

    // 2. Copy exact anchor blocks and align inter-anchor pockets.
    for (b_idx, &(bt, bq, blen)) in blocks.iter().enumerate() {
        for d in 0..blen {
            let idx = (bt + d) - cov_start;
            let q_abs = bq + d;
            col_base[idx] = q.base(q_abs);
            col_rlen[idx] = q.rlen(q_abs);
            col_q_base[idx] = q.q_max(q_abs);
            col_q_len[idx] = q.q_mean(q_abs);
        }
        if b_idx + 1 < blocks.len() {
            let next = blocks[b_idx + 1];
            let t_gap_s = bt + blen;
            let q_gap_s = bq + blen;
            if t_gap_s <= next.0 && q_gap_s <= next.1 {
                align_pocket(
                    t,
                    t_gap_s,
                    next.0,
                    &q,
                    q_gap_s,
                    next.1,
                    cov_start,
                    &mut col_base,
                    &mut col_rlen,
                    &mut col_q_base,
                    &mut col_q_len,
                    &mut ins_after,
                );
            } else {
                // Rare anchor overlap with indel shift: pull back up to `K_HPC` bases so pocket is monotonic.
                let pull = (t_gap_s.saturating_sub(next.0))
                    .max(q_gap_s.saturating_sub(next.1))
                    .add_one_clamped(blen.min(K_HPC));
                let ts2 = (bt + blen).saturating_sub(pull);
                let qs2 = (bq + blen).saturating_sub(pull);
                if ts2 <= next.0 && qs2 <= next.1 {
                    align_pocket(
                        t,
                        ts2,
                        next.0,
                        &q,
                        qs2,
                        next.1,
                        cov_start,
                        &mut col_base,
                        &mut col_rlen,
                        &mut col_q_base,
                        &mut col_q_len,
                        &mut ins_after,
                    );
                }
            }
        }
    }

    // 3. Optional short right dovetail extension after the last exact block.
    if ext_right > 0 {
        let te0 = last.0 + last.2;
        let qe0 = last.1 + last.2;
        align_pocket(
            t,
            te0,
            te0 + ext_right,
            &q,
            qe0,
            qe0 + ext_right,
            cov_start,
            &mut col_base,
            &mut col_rlen,
            &mut col_q_base,
            &mut col_q_len,
            &mut ins_after,
        );
    }

    ins_after.sort_unstable_by_key(|x| x.0);
    ins_after.dedup_by_key(|x| x.0);

    PeerAligned {
        is_rev: ch.is_rev,
        t_start: cov_start as u32,
        t_end: cov_end as u32,
        col_base,
        col_rlen,
        col_q_base,
        col_q_len,
        ins_after,
    }
}

trait AddOneClamped {
    fn add_one_clamped(self, cap: usize) -> usize;
}

impl AddOneClamped for usize {
    #[inline]
    fn add_one_clamped(self, cap: usize) -> usize {
        (self + 2).min(cap)
    }
}

/// Candidate bi-allelic site (`peer 0` = target read `t`, `1..=P` = aligned peers).
struct CandidateSite {
    pos: usize,
    weight: i32,
    self_allele: u8,
    self_q: u8,
    self_strong: bool,
    /// `0` = major `a0`, `1` = minor `a1`, `255` = missing/other; index `0` is target `t`.
    state: [u8; MAX_PEERS + 1],
}

/// Symmetric directional strand-bias test for recurrent methylation / systematic ONT errors.
fn is_strand_biased(f0: u32, r0: u32, f1: u32, r1: u32) -> bool {
    let tot_f = f0 + f1;
    let tot_r = r0 + r1;
    if tot_f < 2 || tot_r < 2 {
        return false;
    }
    let tot = (tot_f + tot_r) as f64;
    let c0 = (f0 + r0) as f64;
    let c1 = (f1 + r1) as f64;

    // One-sided minor or major allele when binomial expectation on the empty strand >= 2.2.
    if (f1 == 0 && c1 * (tot_f as f64) / tot >= 2.2)
        || (r1 == 0 && c1 * (tot_r as f64) / tot >= 2.2)
        || (f0 == 0 && c0 * (tot_f as f64) / tot >= 2.2)
        || (r0 == 0 && c0 * (tot_r as f64) / tot >= 2.2)
    {
        return true;
    }

    // Strong (>5x) strand allele-frequency skew.
    let rate_f = (f1 as f64 + 0.25) / (tot_f as f64 + 0.5);
    let rate_r = (r1 as f64 + 0.25) / (tot_r as f64 + 0.5);
    if (f1 == 0 || r1 == 0 || (f1 + r1) >= 6) && (rate_f > 5.0 * rate_r || rate_r > 5.0 * rate_f)
    {
        return true;
    }
    false
}

/// Pairwise four-gamete compatibility test between two candidate sites `u` and `v`.
/// Returns `Some((edge_bonus, same_index_coupling))` when `u` and `v` co-segregate cleanly.
fn sites_compatible(
    u: &CandidateSite,
    v: &CandidateSite,
    n_actors: usize,
) -> Option<(i32, bool)> {
    if v.pos < u.pos + MIN_PHASE_SEP_HPC {
        return None;
    }
    let mut g = [[0u32; 2]; 2];
    for p in 0..n_actors {
        let su = u.state[p];
        let sv = v.state[p];
        if su <= 1 && sv <= 1 {
            g[su as usize][sv as usize] += 1;
        }
    }
    let n_span = g[0][0] + g[0][1] + g[1][0] + g[1][1];
    if n_span < 6 {
        return None;
    }
    let in_phase = g[0][0] + g[1][1];
    let anti_phase = g[0][1] + g[1][0];
    let coupling = in_phase >= anti_phase;
    let (h_a, h_b, disc) = if coupling {
        (g[0][0], g[1][1], anti_phase)
    } else {
        (g[0][1], g[1][0], in_phase)
    };
    let min_hap = h_a.min(h_b);
    let minor_u = (g[1][0] + g[1][1]).min(g[0][0] + g[0][1]);
    let minor_v = (g[0][1] + g[1][1]).min(g[0][0] + g[1][0]);
    let min_minor = minor_u.min(minor_v);

    if min_hap < 2
        || (min_hap as f64) < 0.60 * (min_minor as f64)
        || (disc as f64) > (0.12 * (n_span as f64)).max(1.0)
    {
        return None;
    }
    let bonus = ((min_hap * 8).saturating_sub(disc * 10) as i32).max(8);
    Some((bonus, coupling))
}

/// Whether candidate site `v` strongly conflicts (`disc > 20%` over `>= 8` spanning reads) with an
/// already-accepted multi-site haplotype chain.
fn conflicts_with_accepted(
    v: &CandidateSite,
    accepted: &[usize],
    sites: &[CandidateSite],
    n_actors: usize,
) -> bool {
    for &idx in accepted {
        let u = &sites[idx];
        if u.pos.abs_diff(v.pos) < MIN_PHASE_SEP_HPC {
            continue;
        }
        let mut g = [[0u32; 2]; 2];
        for p in 0..n_actors {
            let su = u.state[p];
            let sv = v.state[p];
            if su <= 1 && sv <= 1 {
                g[su as usize][sv as usize] += 1;
            }
        }
        let n_span = g[0][0] + g[0][1] + g[1][0] + g[1][1];
        if n_span >= 8 {
            let disc = (g[0][1] + g[1][0]).min(g[0][0] + g[1][1]);
            if (disc as f64) > 0.20 * (n_span as f64) {
                return true;
            }
        }
    }
    false
}

struct ReadCorrectionOut {
    read: Option<HpcRead>,
    peers_used: usize,
    cand_sites: u64,
    pruned_strand_bias: u64,
    pruned_incompatible: u64,
    phased_kept: u64,
    base_edits: u64,
    hp_len_edits: u64,
}

/// Steps 2-4 of `hifiasm (--ont)` on a single target read `t`.
fn correct_single_read(
    t_id: usize,
    reads: &[HpcRead],
    index: &HpcIndex,
    lut: &[f64; 61],
) -> ReadCorrectionOut {
    let t = &reads[t_id];
    let t_len = t.hpc_len();
    let chains = select_balanced_peers(
        collect_chains_for_read(t_id, reads, index, MIN_OVLP_HPC, true),
        t_len,
    );
    if chains.is_empty() {
        return ReadCorrectionOut {
            read: Some(t.clone()),
            peers_used: 0,
            cand_sites: 0,
            pruned_strand_bias: 0,
            pruned_incompatible: 0,
            phased_kept: 0,
            base_edits: 0,
            hp_len_edits: 0,
        };
    }

    let peers: Vec<PeerAligned> = chains
        .iter()
        .map(|ch| align_peer_chain(t, &reads[ch.peer_id as usize], ch))
        .collect();
    let n_peers = peers.len();
    let n_actors = n_peers + 1;

    // -------------------------------------------------------------------------
    // Step 2: Candidate informative site discovery + 3-way pre-filtering
    // -------------------------------------------------------------------------
    let mut cand_sites_raw = 0u64;
    let mut pruned_strand_bias = 0u64;
    let mut sites: Vec<CandidateSite> = Vec::new();

    for i in 0..t_len {
        let near_long_hp = t.rlen[i] >= 7
            || (i > 0 && t.rlen[i - 1] >= 7)
            || (i + 1 < t_len && t.rlen[i + 1] >= 7);

        let mut cnt_f = [0u32; 5];
        let mut cnt_r = [0u32; 5];
        let mut q14_cnt = [0u32; 5];
        let mut q18_f = [0u32; 5];
        let mut q18_r = [0u32; 5];
        let mut q_sum = [0u32; 5];

        let tb = t.base[i];
        let tq = t.q_max[i];
        if tb < 4 && tq >= MIN_SITE_Q {
            let a = tb as usize;
            cnt_f[a] += 1;
            q_sum[a] += u32::from(tq);
            if tq >= 14 {
                q14_cnt[a] += 1;
            }
            if tq >= 18 {
                q18_f[a] += 1;
            }
        }

        for p in &peers {
            if let Some((b, _, qb, _)) = p.call_at(i) {
                if qb < MIN_SITE_Q || b as usize >= 5 {
                    continue;
                }
                if near_long_hp && b == DEL_BASE && qb < 15 {
                    continue;
                }
                let a = b as usize;
                if p.is_rev {
                    cnt_r[a] += 1;
                    if qb >= 18 {
                        q18_r[a] += 1;
                    }
                } else {
                    cnt_f[a] += 1;
                    if qb >= 18 {
                        q18_f[a] += 1;
                    }
                }
                q_sum[a] += u32::from(qb);
                if qb >= 14 {
                    q14_cnt[a] += 1;
                }
            }
        }

        // Top two alleles `(a0, a1)` at column `i`.
        let mut order = [0usize, 1, 2, 3, 4];
        order.sort_unstable_by_key(|&a| std::cmp::Reverse((cnt_f[a] + cnt_r[a], q_sum[a])));
        let (a0, a1) = (order[0], order[1]);
        let c0 = cnt_f[a0] + cnt_r[a0];
        let c1 = cnt_f[a1] + cnt_r[a1];
        let cov = c0 + c1;
        if cov < 6 || c1 < MIN_MINOR_COUNT || (c1 as f64) < MIN_MINOR_FREQ * (cov as f64) {
            continue;
        }
        if q14_cnt[a1] < 2 {
            continue;
        }
        // Long-homopolymer boundary deletion artifact guard.
        if near_long_hp && (a1 == DEL_BASE as usize || a0 == DEL_BASE as usize) {
            let del_a = if a1 == DEL_BASE as usize { a1 } else { a0 };
            if q18_f[del_a] < 2 || q18_r[del_a] < 2 {
                continue;
            }
        }

        cand_sites_raw += 1;
        if is_strand_biased(cnt_f[a0], cnt_r[a0], cnt_f[a1], cnt_r[a1]) {
            pruned_strand_bias += 1;
            continue;
        }

        let mut state = [255u8; MAX_PEERS + 1];
        let self_allele = if tb as usize == a0 && tq >= MIN_SITE_Q {
            state[0] = 0;
            0
        } else if tb as usize == a1 && tq >= MIN_SITE_Q {
            state[0] = 1;
            1
        } else {
            255
        };
        for (idx, p) in peers.iter().enumerate() {
            if let Some((b, _, qb, _)) = p.call_at(i) {
                if qb >= MIN_SITE_Q {
                    if b as usize == a0 {
                        state[idx + 1] = 0;
                    } else if b as usize == a1 {
                        state[idx + 1] = 1;
                    }
                }
            }
        }
        let self_strong = c1 >= 4
            && q18_f[a0] >= 2
            && q18_r[a0] >= 2
            && q18_f[a1] >= 2
            && q18_r[a1] >= 2;
        sites.push(CandidateSite {
            pos: i * 2,
            weight: (q_sum[a1] as i32 / 8).max(4),
            self_allele,
            self_q: tq,
            self_strong,
            state,
        });

        // Also check non-HP insertion pocket immediately after column `i` (`NO_INS` vs `INS`).
        let mut ins_f = 0u32;
        let mut ins_r = 0u32;
        let mut no_f = 1u32; // target `t` has `NO_INS`
        let mut no_r = 0u32;
        let mut ins_q14 = 0u32;
        let mut ins_q18_f = 0u32;
        let mut ins_q18_r = 0u32;
        let mut ins_q_sum = 0u32;
        let mut ins_state = [255u8; MAX_PEERS + 1];
        ins_state[0] = 0;

        for (idx, p) in peers.iter().enumerate() {
            if !p.spans_after(i) {
                continue;
            }
            if let Some(si) = p.ins_at(i as u32) {
                let qb = si.as_slice()[0].q_base;
                if qb < MIN_SITE_Q || (near_long_hp && qb < 15) {
                    continue;
                }
                ins_state[idx + 1] = 1;
                ins_q_sum += u32::from(qb);
                if qb >= 14 {
                    ins_q14 += 1;
                }
                if p.is_rev {
                    ins_r += 1;
                    if qb >= 18 {
                        ins_q18_r += 1;
                    }
                } else {
                    ins_f += 1;
                    if qb >= 18 {
                        ins_q18_f += 1;
                    }
                }
            } else {
                ins_state[idx + 1] = 0;
                if p.is_rev {
                    no_r += 1;
                } else {
                    no_f += 1;
                }
            }
        }
        let c_no = no_f + no_r;
        let c_ins = ins_f + ins_r;
        let c_minor = c_no.min(c_ins);
        if c_no + c_ins >= 6
            && c_minor >= MIN_MINOR_COUNT
            && (c_minor as f64) >= MIN_MINOR_FREQ * ((c_no + c_ins) as f64)
            && ins_q14 >= 2
            && (!near_long_hp || (ins_q18_f >= 2 && ins_q18_r >= 2))
        {
            cand_sites_raw += 1;
            if is_strand_biased(no_f, no_r, ins_f, ins_r) {
                pruned_strand_bias += 1;
            } else {
                let strong = c_minor >= 4
                    && no_f >= 2
                    && no_r >= 2
                    && ins_q18_f >= 2
                    && ins_q18_r >= 2;
                sites.push(CandidateSite {
                    pos: i * 2 + 1,
                    weight: (ins_q_sum as i32 / 8).max(4),
                    self_allele: 0,
                    self_q: tq,
                    self_strong: strong,
                    state: ins_state,
                });
            }
        }
    }

    if sites.len() > MAX_CANDIDATE_SITES {
        sites.sort_unstable_by_key(|s| std::cmp::Reverse(s.weight));
        sites.truncate(MAX_CANDIDATE_SITES);
        sites.sort_unstable_by_key(|s| s.pos);
    }

    // -------------------------------------------------------------------------
    // Step 3: O(n^2)-time, O(n)-space DP compatibility phasing across sites
    // -------------------------------------------------------------------------
    let n_sites = sites.len();
    let mut active = vec![true; n_sites];
    let mut target_hap = vec![255u8; n_sites];
    let mut multi_chain_sites: Vec<usize> = Vec::new();
    let mut dp_score = vec![0i32; n_sites];
    let mut dp_prev = vec![usize::MAX; n_sites];
    let mut dp_coup = vec![true; n_sites];

    loop {
        let mut best_v = usize::MAX;
        let mut best_sc = i32::MIN;
        for v in 0..n_sites {
            if !active[v] {
                continue;
            }
            let mut sc = sites[v].weight;
            let mut pr = usize::MAX;
            let mut cp = true;
            for u in 0..v {
                if !active[u] {
                    continue;
                }
                if let Some((bonus, coupling)) = sites_compatible(&sites[u], &sites[v], n_actors) {
                    let cand = dp_score[u] + bonus + sites[v].weight;
                    if cand > sc {
                        sc = cand;
                        pr = u;
                        cp = coupling;
                    }
                }
            }
            dp_score[v] = sc;
            dp_prev[v] = pr;
            dp_coup[v] = cp;
            if sc > best_sc {
                best_sc = sc;
                best_v = v;
            }
        }
        if best_v == usize::MAX {
            break;
        }

        let mut chain: Vec<(usize, bool)> = Vec::new();
        let mut cur = best_v;
        while cur != usize::MAX {
            chain.push((cur, dp_coup[cur]));
            cur = dp_prev[cur];
        }
        chain.reverse();

        if chain.len() >= 2 {
            // Reconstruct haplotype H0 (`rel_allele[r]`) along the chain and vote target `t`
            // between H0 and H1 = 1 - H0 so a single error on `t` at one phased SNP never flips
            // `t`'s haplotype at that site!
            let mut h0 = vec![0u8; chain.len()];
            for r in 1..chain.len() {
                let coup_from_prev = chain[r].1;
                h0[r] = if coup_from_prev { h0[r - 1] } else { 1 - h0[r - 1] };
            }
            let (mut vote_h0, mut vote_h1) = (0u32, 0u32);
            for (r, &(s_idx, _)) in chain.iter().enumerate() {
                active[s_idx] = false;
                let sa = sites[s_idx].self_allele;
                if sa <= 1 {
                    let w = u32::from(sites[s_idx].self_q.max(10));
                    if sa == h0[r] {
                        vote_h0 += w;
                    } else {
                        vote_h1 += w;
                    }
                }
            }
            let pick_h1 = vote_h1 > vote_h0;
            for (r, &(s_idx, _)) in chain.iter().enumerate() {
                target_hap[s_idx] = if pick_h1 { 1 - h0[r] } else { h0[r] };
                multi_chain_sites.push(s_idx);
            }
        } else {
            // Single isolated site: keep only if `self_strong` and not discordant with an accepted
            // multi-site chain; otherwise reject as an unlinked error site.
            let s_idx = chain[0].0;
            active[s_idx] = false;
            if sites[s_idx].self_strong
                && sites[s_idx].self_allele <= 1
                && !conflicts_with_accepted(&sites[s_idx], &multi_chain_sites, &sites, n_actors)
            {
                target_hap[s_idx] = sites[s_idx].self_allele;
            }
        }
    }

    let phased_kept = target_hap.iter().filter(|&&h| h <= 1).count() as u64;
    let pruned_incompatible = (n_sites as u64).saturating_sub(phased_kept);

    // -------------------------------------------------------------------------
    // Step 4: In-phase peer selection + strand-balanced Phred-weighted consensus
    // -------------------------------------------------------------------------
    let mut peer_in_phase = vec![true; n_peers];
    let mut peer_diff_cnt = vec![0u16; n_peers];
    if phased_kept > 0 {
        for p_idx in 0..n_peers {
            let (mut same, mut diff) = (0u32, 0u32);
            for (s_idx, site) in sites.iter().enumerate() {
                let th = target_hap[s_idx];
                if th > 1 {
                    continue;
                }
                let ps = site.state[p_idx + 1];
                if ps <= 1 {
                    if ps == th {
                        same += 1;
                    } else {
                        diff += 1;
                    }
                }
            }
            peer_diff_cnt[p_idx] = diff.min(u16::MAX as u32) as u16;
            let span = same + diff;
            peer_in_phase[p_idx] = if span == 0 {
                true
            } else if span == 1 {
                diff == 0
            } else {
                diff == 0 || (diff == 1 && same >= 4) || (diff as f64) <= 0.12 * (span as f64)
            };
        }
    }

    let mut raw_seq: Vec<u8> = Vec::with_capacity(t.raw_len + 64);
    let mut raw_qual: Vec<u8> = Vec::with_capacity(t.raw_len + 64);
    let mut base_edits = 0u64;
    let mut hp_len_edits = 0u64;
    let mut rlen_obs: Vec<(u16, f64, bool)> = Vec::with_capacity(n_actors);

    for i in 0..t_len {
        let mut in_phase_cov = 0usize;
        for (p_idx, p) in peers.iter().enumerate() {
            if peer_in_phase[p_idx] && p.call_at(i).is_some() {
                in_phase_cov += 1;
            }
        }
        let relax = in_phase_cov < 2;

        let mut wf = [0.0f64; 5];
        let mut wr = [0.0f64; 5];
        let mut nf = 0u32;
        let mut nr = 0u32;

        let tb = t.base[i];
        let tq = t.q_max[i];
        if tb < 4 {
            let self_scale = if in_phase_cov >= 4 { 0.75 } else { 1.0 };
            wf[tb as usize] += f64::from(tq.min(40)) * self_scale;
            nf += 1;
        }

        for (p_idx, p) in peers.iter().enumerate() {
            if !peer_in_phase[p_idx] && !(relax && peer_diff_cnt[p_idx] <= 1) {
                continue;
            }
            if let Some((b, _, qb, _)) = p.call_at(i) {
                if (b as usize) < 5 {
                    let w = f64::from(qb.min(40));
                    if p.is_rev {
                        wr[b as usize] += w;
                        nr += 1;
                    } else {
                        wf[b as usize] += w;
                        nf += 1;
                    }
                }
            }
        }

        if nf + nr == 0 {
            // Preserve ambiguous `N` run if completely uncovered.
            raw_seq.push(b'N');
            raw_qual.push(2);
            continue;
        }

        let both_strands = nf >= 2 && nr >= 2;
        let mut best_a = if tb < 4 { tb as usize } else { 0 };
        let mut best_sc = -1.0f64;
        let mut runner_sc = 0.0f64;
        for a in 0..5 {
            let sc = if both_strands {
                wf[a] / f64::from(nf) + wr[a] / f64::from(nr)
            } else {
                wf[a] + wr[a]
            };
            if sc > best_sc {
                runner_sc = best_sc.max(0.0);
                best_sc = sc;
                best_a = a;
            } else if sc > runner_sc {
                runner_sc = sc;
            }
        }

        if best_a == DEL_BASE as usize {
            base_edits += 1;
        } else {
            if tb as usize != best_a {
                base_edits += 1;
            }
            rlen_obs.clear();
            if tb as usize == best_a {
                let self_scale = if in_phase_cov >= 4 { 0.75 } else { 1.0 };
                rlen_obs.push((t.rlen[i], f64::from(t.q_mean[i].min(35)) * self_scale, false));
            }
            for (p_idx, p) in peers.iter().enumerate() {
                if !peer_in_phase[p_idx] && !(relax && peer_diff_cnt[p_idx] <= 1) {
                    continue;
                }
                if let Some((b, rl, _, ql)) = p.call_at(i) {
                    if b as usize == best_a && rl >= 1 {
                        rlen_obs.push((rl, f64::from(ql.min(35)), p.is_rev));
                    }
                }
            }
            let best_rlen = call_homopolymer_length(&rlen_obs, t.rlen[i]);
            if tb as usize == best_a && best_rlen != t.rlen[i] {
                hp_len_edits += best_rlen.abs_diff(t.rlen[i]) as u64;
            }
            let conf_q = (((best_sc - runner_sc) / (best_sc + 1e-6)) * 38.0)
                .round()
                .clamp(12.0, 40.0) as u8;
            raw_seq.resize(raw_seq.len() + best_rlen as usize, BASES[best_a]);
            raw_qual.resize(raw_qual.len() + best_rlen as usize, conf_q);
        }

        // Consensus check for an HPC insertion immediately after column `i`.
        let (mut no_wf, mut no_wr, mut ins_wf, mut ins_wr) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let (mut sf, mut sr, mut ins_cnt) = (0u32, 0u32, 0u32);
        if i + 1 < t_len {
            let self_scale = if in_phase_cov >= 4 { 0.75 } else { 1.0 };
            no_wf += f64::from(tq.min(t.q_max[i + 1]).min(35)) * self_scale;
            sf += 1;
        }
        let mut ins_candidates: Vec<(&SmallIns, f64)> = Vec::new();
        for (p_idx, p) in peers.iter().enumerate() {
            if !peer_in_phase[p_idx] && !(relax && peer_diff_cnt[p_idx] <= 1) {
                continue;
            }
            if !p.spans_after(i) {
                continue;
            }
            if let Some(si) = p.ins_at(i as u32) {
                let w = f64::from(si.as_slice()[0].q_base.min(40));
                ins_cnt += 1;
                ins_candidates.push((si, w));
                if p.is_rev {
                    ins_wr += w;
                    sr += 1;
                } else {
                    ins_wf += w;
                    sf += 1;
                }
            } else {
                let w = 22.0;
                if p.is_rev {
                    no_wr += w;
                    sr += 1;
                } else {
                    no_wf += w;
                    sf += 1;
                }
            }
        }
        if ins_cnt >= 2 && sf + sr >= 3 {
            let (sc_ins, sc_no) = if sf >= 2 && sr >= 2 {
                (
                    ins_wf / f64::from(sf) + ins_wr / f64::from(sr),
                    no_wf / f64::from(sf) + no_wr / f64::from(sr),
                )
            } else {
                (ins_wf + ins_wr, no_wf + no_wr)
            };
            if sc_ins > sc_no {
                // Emit the highest-weight modal inserted HPC run sequence.
                let chosen = choose_modal_insertion(&ins_candidates);
                base_edits += u64::from(chosen.len);
                for run in chosen.as_slice() {
                    if run.base < 4 {
                        raw_seq.resize(raw_seq.len() + run.rlen as usize, BASES[run.base as usize]);
                        raw_qual.resize(raw_qual.len() + run.rlen as usize, run.q_base.clamp(14, 36));
                    }
                }
            }
        }
    }

    let next_read = encode_hpc_read(t.id.clone(), &raw_seq, &raw_qual, lut);
    let in_phase_total = peer_in_phase.iter().filter(|&&k| k).count();
    ReadCorrectionOut {
        read: next_read,
        peers_used: in_phase_total,
        cand_sites: cand_sites_raw,
        pruned_strand_bias,
        pruned_incompatible,
        phased_kept,
        base_edits,
        hp_len_edits,
    }
}

/// Strand-balanced Phred-weighted mode for homopolymer run lengths, breaking close ties toward the
/// weighted median so ONT deletion skew on long runs does not bias short.
fn call_homopolymer_length(obs: &[(u16, f64, bool)], fallback: u16) -> u16 {
    if obs.is_empty() {
        return fallback.max(1);
    }
    if obs.len() == 1 {
        return obs[0].0.max(1);
    }
    let max_l = obs.iter().map(|x| x.0).max().unwrap_or(1).min(64) as usize;
    if max_l == 1 {
        return 1;
    }
    let mut nf = 0u32;
    let mut nr = 0u32;
    for &(_, _, is_rev) in obs {
        if is_rev {
            nr += 1;
        } else {
            nf += 1;
        }
    }
    let both = nf >= 2 && nr >= 2;
    let mut wf = vec![0.0f64; max_l + 1];
    let mut wr = vec![0.0f64; max_l + 1];
    let mut tot_w = vec![0.0f64; max_l + 1];

    for &(rl, w, is_rev) in obs {
        let l = (rl as usize).min(max_l);
        if is_rev {
            wr[l] += w;
        } else {
            wf[l] += w;
        }
    }
    let mut sum_all = 0.0f64;
    for l in 1..=max_l {
        tot_w[l] = if both {
            wf[l] / f64::from(nf) + wr[l] / f64::from(nr)
        } else {
            wf[l] + wr[l]
        };
        sum_all += tot_w[l];
    }
    let mut med = 1usize;
    let mut acc = 0.0f64;
    for (l, &w) in tot_w.iter().enumerate().skip(1) {
        acc += w;
        if acc >= 0.5 * sum_all {
            med = l;
            break;
        }
    }
    let mut best_l = med;
    let mut best_w = tot_w[med];
    for (l, &w) in tot_w.iter().enumerate().skip(1) {
        if w > best_w * 1.03 || (w >= best_w * 0.97 && l.abs_diff(med) < best_l.abs_diff(med)) {
            best_w = w;
            best_l = l;
        }
    }
    best_l as u16
}

fn choose_modal_insertion(cands: &[(&SmallIns, f64)]) -> SmallIns {
    let mut best = *cands[0].0;
    let mut best_w = 0.0f64;
    for &(si, _) in cands {
        let s_slice = si.as_slice();
        let w_base: f64 = cands
            .iter()
            .filter(|&&(other, _)| {
                other.len == si.len
                    && other
                        .as_slice()
                        .iter()
                        .zip(s_slice)
                        .all(|(a, b)| a.base == b.base)
            })
            .map(|&(_, w)| w)
            .sum();
        if w_base > best_w {
            best_w = w_base;
            best = *si;
        }
    }
    // Refine run lengths across matching-base insertions using each run's length quality `q_len`.
    for pos in 0..(best.len as usize) {
        let b_want = best.runs[pos].base;
        let obs: Vec<(u16, f64, bool)> = cands
            .iter()
            .filter(|&&(other, _)| other.len == best.len && other.runs[pos].base == b_want)
            .map(|&(other, w)| {
                let ql = f64::from(other.runs[pos].q_len.max(8).min(35));
                (other.runs[pos].rlen, 0.5 * (w + ql), false)
            })
            .collect();
        best.runs[pos].rlen = call_homopolymer_length(&obs, best.runs[pos].rlen);
    }
    best
}

/// Execute 2 iterative rounds of `hifiasm (--ont)` phased overlap error correction.
pub fn run_r10(
    reads_in: Vec<Rec>,
    profile: ChemistryProfile,
) -> Result<(Vec<Rec>, Vec<RoundReport>, R10Report)> {
    let lut = phred_err_lut();
    let total_reads_in = reads_in.len() as u64;
    let total_bases_in: u64 = reads_in.iter().map(|r| r.seq.len() as u64).sum();
    let has_real_qual = profile.qual_base_fraction >= AUTO_MIN_QUAL_FRAC;

    // Initial read QC + HPC encoding.
    let mut hpc_reads: Vec<HpcRead> = reads_in
        .into_par_iter()
        .filter_map(|r| {
            if r.seq.len() < MIN_RAW_LEN {
                return None;
            }
            let qp: Vec<u8> = match r.qual.as_deref() {
                Some(q) if q.len() == r.seq.len() => {
                    if has_real_qual && effective_q_of_slice(q, &lut) < MIN_READ_Q_EFF {
                        return None;
                    }
                    q.iter().map(|&b| decode_phred(b)).collect()
                }
                _ => vec![DEFAULT_SYNTH_Q; r.seq.len()],
            };
            encode_hpc_read(r.id, &r.seq, &qp, &lut)
        })
        .collect();

    let qc_kept = hpc_reads.len() as u64;
    let max_occ = ((profile.solid_21mer_peak as usize).saturating_mul(8)).clamp(400, 2400);

    let mut r10_rounds = Vec::with_capacity(2);
    let mut compat_rounds = Vec::with_capacity(2);

    for round_idx in 1..=2 {
        let t0 = Instant::now();
        let r_in_cnt = if round_idx == 1 {
            total_reads_in
        } else {
            hpc_reads.len() as u64
        };
        let b_in_cnt: u64 = if round_idx == 1 {
            total_bases_in
        } else {
            hpc_reads.iter().map(|r| r.raw_len as u64).sum()
        };

        let (clipped_reads, chimera_cnt, chimera_trimmed) = if round_idx == 1 {
            clip_chimeras(hpc_reads, max_occ)
        } else {
            (hpc_reads, 0, 0)
        };

        let index = HpcIndex::build(&clipped_reads, max_occ);
        let outs: Vec<ReadCorrectionOut> = (0..clipped_reads.len())
            .into_par_iter()
            .map(|t_id| correct_single_read(t_id, &clipped_reads, &index, &lut))
            .collect();

        let mut next_reads = Vec::with_capacity(clipped_reads.len());
        let mut sum_peers = 0usize;
        let mut cand_sites = 0u64;
        let mut pruned_sb = 0u64;
        let mut pruned_incomp = 0u64;
        let mut phased_kept = 0u64;
        let mut base_edits = 0u64;
        let mut hp_edits = 0u64;

        for o in outs {
            sum_peers += o.peers_used;
            cand_sites += o.cand_sites;
            pruned_sb += o.pruned_strand_bias;
            pruned_incomp += o.pruned_incompatible;
            phased_kept += o.phased_kept;
            base_edits += o.base_edits;
            hp_edits += o.hp_len_edits;
            if let Some(r) = o.read.filter(|r| r.raw_len >= MIN_RAW_LEN) {
                next_reads.push(r);
            }
        }

        let b_out_cnt: u64 = next_reads.iter().map(|r| r.raw_len as u64).sum();
        let mean_peers = if !clipped_reads.is_empty() {
            sum_peers as f64 / clipped_reads.len() as f64
        } else {
            0.0
        };
        let secs = t0.elapsed().as_secs_f64();

        eprintln!(
            "[ovasm] correct r10 round {round_idx}: {} -> {} reads ({:.1} Mb), mean peers {:.1}; sites {} -> {} phased ({} strand-bias, {} incompatible pruned); {} base edits, {} hp edits, {} chimeras clipped; {:.1}s",
            r_in_cnt,
            next_reads.len(),
            b_out_cnt as f64 / 1e6,
            mean_peers,
            cand_sites,
            phased_kept,
            pruned_sb,
            pruned_incomp,
            base_edits,
            hp_edits,
            chimera_cnt,
            secs
        );

        r10_rounds.push(R10RoundReport {
            round: round_idx,
            reads_in: r_in_cnt,
            reads_qc_kept: if round_idx == 1 { qc_kept } else { r_in_cnt },
            reads_chimera_clipped: chimera_cnt,
            bases_chimera_trimmed: chimera_trimmed,
            reads_out: next_reads.len() as u64,
            bases_in: b_in_cnt,
            bases_out: b_out_cnt,
            mean_peers,
            candidate_sites: cand_sites,
            pruned_strand_bias_sites: pruned_sb,
            pruned_incompatible_sites: pruned_incomp,
            phased_sites_kept: phased_kept,
            base_edits,
            hp_length_edits: hp_edits,
            seconds: secs,
        });

        compat_rounds.push(RoundReport {
            k: K_HPC,
            distinct_kmers_ge2: index.table.len(),
            min_count: 2,
            min_count_auto: true,
            solid_kmers: index.table.len(),
            solid_peak: profile.solid_21mer_peak,
            hybrid: false,
            reads_in: r_in_cnt,
            reads_out: next_reads.len() as u64,
            stats: ReadStats {
                bases_in: b_in_cnt,
                bases_out: b_out_cnt,
                anchored_bases: b_out_cnt,
                gaps: base_edits + hp_edits,
                gaps_bridged: base_edits + hp_edits,
                trimmed_bases: chimera_trimmed,
            },
            seconds: secs,
        });

        hpc_reads = next_reads;
    }

    let out_recs: Vec<Rec> = hpc_reads
        .into_par_iter()
        .map(|r| {
            let seq = r.expand_fasta();
            Rec {
                id: r.id,
                seq,
                qual: None,
            }
        })
        .collect();

    Ok((
        out_recs,
        compat_rounds,
        R10Report {
            profile,
            rounds: r10_rounds,
        },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 11
        }
        fn unit(&mut self) -> f64 {
            self.next() as f64 / (1u64 << 53) as f64
        }
        fn base(&mut self) -> u8 {
            BASES[(self.next() >> 51) as usize]
        }
    }

    fn random_genome(len: usize, seed: u64) -> Vec<u8> {
        let mut r = Rng(seed);
        (0..len).map(|_| r.base()).collect()
    }

    fn revcomp_qual(seq: &[u8], qual: &[u8]) -> (Vec<u8>, Vec<u8>) {
        let rc_s = crate::evidence::revcomp(seq);
        let mut rc_q = qual.to_vec();
        rc_q.reverse();
        (rc_s, rc_q)
    }

    /// Semi-global edit distance of `q` against circular `genome` (both strands).
    fn best_edit_distance(q: &[u8], genome: &[u8]) -> usize {
        let circ = [genome, &genome[..q.len().min(genome.len())]].concat();
        let rc = crate::evidence::revcomp(&circ);
        let semi = |query: &[u8], target: &[u8]| -> usize {
            let mut prev = vec![0usize; target.len() + 1];
            let mut cur = vec![0usize; target.len() + 1];
            for &qb in query {
                cur[0] = prev[0] + 1;
                for (j, &tb) in target.iter().enumerate() {
                    let sub = prev[j] + usize::from(qb != tb);
                    cur[j + 1] = sub.min(prev[j + 1] + 1).min(cur[j] + 1);
                }
                std::mem::swap(&mut prev, &mut cur);
            }
            *prev.iter().min().unwrap()
        };
        semi(q, &circ).min(semi(q, &rc))
    }

    fn sim_r10_read(truth: &[u8], err_rate: f64, rng: &mut Rng) -> (Vec<u8>, Vec<u8>) {
        let mut seq = Vec::with_capacity(truth.len() + 32);
        let mut qual = Vec::with_capacity(truth.len() + 32);
        for &b in truth {
            let u = rng.unit();
            if u < err_rate * 0.45 {
                // Substitution error with lower Q (Q9-Q13).
                let mut alt = rng.base();
                if alt == b {
                    alt = BASES[((base_code(b).unwrap() + 1) % 4) as usize];
                }
                seq.push(alt);
                qual.push(33 + 10 + (rng.next() % 4) as u8);
            } else if u < err_rate * 0.75 {
                // 1-bp insertion (often homopolymer slip) at Q10-Q14.
                seq.push(b);
                qual.push(33 + 24);
                seq.push(if rng.unit() < 0.7 { b } else { rng.base() });
                qual.push(33 + 11);
            } else if u < err_rate {
                // 1-bp deletion.
                continue;
            } else {
                seq.push(b);
                qual.push(33 + 22 + (rng.next() % 8) as u8);
            }
        }
        (seq, qual)
    }

    #[test]
    fn chemistry_auto_distinguishes_r10_fastq_from_r9_dorado_q90_and_fasta() {
        let g = random_genome(6_000, 17);
        let circ = [g.as_slice(), &g[..2_500]].concat();
        let mut rng = Rng(101);

        // 1. R10.4.1 FASTQ (~1.2% error, calibrated Q22-Q29).
        let r10_reads: Vec<Rec> = (0..120)
            .map(|i| {
                let st = (i * 50) % 6_000;
                let (s, q) = sim_r10_read(&circ[st..st + 2_200], 0.012, &mut rng);
                let (s, q) = if i % 2 == 1 { revcomp_qual(&s, &q) } else { (s, q) };
                Rec {
                    id: format!("r10_{i}").into_bytes(),
                    seq: s,
                    qual: Some(q),
                }
            })
            .collect();
        let p_r10 = detect_chemistry(&r10_reads, OntChemistry::Auto);
        assert!(p_r10.use_r10, "expected r10: {p_r10:?}");
        assert_eq!(p_r10.selected, "r10");

        // 2. Quality-less FASTA of identical sequences -> must keep Legacy.
        let fasta_reads: Vec<Rec> = r10_reads
            .iter()
            .map(|r| Rec {
                id: r.id.clone(),
                seq: r.seq.clone(),
                qual: None,
            })
            .collect();
        let p_fa = detect_chemistry(&fasta_reads, OntChemistry::Auto);
        assert!(!p_fa.use_r10, "FASTA must stay legacy: {p_fa:?}");

        // 3. Dorado R9 FASTQ (~9% error with Q90 `{` spikes alongside Q3 `$` noisy bases):
        // arithmetic mean Q would exceed 25, but Q_eff and solid 21-mer fraction stay R9!
        let r9_reads: Vec<Rec> = (0..120)
            .map(|i| {
                let st = (i * 50) % 6_000;
                let mut seq = Vec::with_capacity(2_200);
                let mut qual = Vec::with_capacity(2_200);
                for &b in &circ[st..st + 2_200] {
                    if rng.unit() < 0.09 {
                        let alt = BASES[((base_code(b).unwrap() + 1) % 4) as usize];
                        seq.push(alt);
                        qual.push(b'$'); // Q3
                    } else {
                        seq.push(b);
                        qual.push(b'{'); // Q90 spike!
                    }
                }
                Rec {
                    id: format!("r9_{i}").into_bytes(),
                    seq,
                    qual: Some(qual),
                }
            })
            .collect();
        let p_r9 = detect_chemistry(&r9_reads, OntChemistry::Auto);
        assert!(!p_r9.use_r10, "R9 Dorado Q90 must stay legacy: {p_r9:?}");

        // 4. PacBio HiFi profile (Q ~27-30, ~0.15% error => solid 21-mer ~0.97), even with SRA-renamed
        // headers lacking `/ccs`, must NOT trigger R10 under Auto (`--chemistry r10` still forces it).
        let hifi_reads: Vec<Rec> = (0..80)
            .map(|i| {
                let st = (i * 75) % 6_000;
                let mut seq = Vec::with_capacity(2_200);
                let mut qual = Vec::with_capacity(2_200);
                for &b in &circ[st..st + 2_200] {
                    if rng.unit() < 0.0015 {
                        let alt = BASES[((base_code(b).unwrap() + 1) % 4) as usize];
                        seq.push(alt);
                        qual.push(33 + 18);
                    } else {
                        seq.push(b);
                        qual.push(33 + 30);
                    }
                }
                Rec {
                    id: format!("ERR12779266.{}", i + 1).into_bytes(),
                    seq,
                    qual: Some(qual),
                }
            })
            .collect();
        let p_hifi = detect_chemistry(&hifi_reads, OntChemistry::Auto);
        assert!(!p_hifi.use_r10, "HiFi reads passed to Auto must stay legacy: {p_hifi:?}");
        assert!(p_hifi.effective_phred_q.unwrap() >= 25.0 && p_hifi.solid_21mer_fraction >= 0.95);
        let p_hifi_forced = detect_chemistry(&hifi_reads, OntChemistry::R10);
        assert!(p_hifi_forced.use_r10, "explicit --chemistry r10 must still force R10");

        // 5. PacBio `/ccs` headers at lower Q22 also stay legacy under Auto.
        let ccs_reads: Vec<Rec> = r10_reads
            .iter()
            .enumerate()
            .map(|(i, r)| Rec {
                id: format!("m64011_220101_000000/{}/ccs", i + 1).into_bytes(),
                seq: r.seq.clone(),
                qual: r.qual.clone(),
            })
            .collect();
        let p_ccs = detect_chemistry(&ccs_reads, OntChemistry::Auto);
        assert!(!p_ccs.use_r10, "CCS-header reads must stay legacy under Auto: {p_ccs:?}");
    }

    #[test]
    fn strand_biased_methylation_errors_are_pruned_and_corrected_under_skewed_strands() {
        // 5 kb genome; at 5 CpG loci, 75% of forward (+) reads carry a systematic C->T error at
        // Q16 under a 3:1 forward:reverse strand skew (so raw T outnumbers C overall: 45% vs 55%,
        // or even >50% on forward!). Reverse (-) reads have normal 1% random noise.
        let mut g = random_genome(5_000, 303);
        let meth_sites = [800usize, 1_500, 2_200, 3_100, 3_900];
        for &pos in &meth_sites {
            g[pos - 1] = b'A';
            g[pos] = b'C';
            g[pos + 1] = b'G';
            g[pos + 2] = b'T';
        }
        let circ = [g.as_slice(), &g[..2_500]].concat();
        let mut rng = Rng(404);
        let mut reads = Vec::new();
        for i in 0..160usize {
            let st = (i * 31) % 5_000;
            let mut truth = circ[st..st + 2_200].to_vec();
            let is_fwd = i % 4 != 0; // 3:1 forward:reverse strand skew!
            if is_fwd {
                for &mp in &meth_sites {
                    for off in [mp, mp + 5_000] {
                        if off >= st && off < st + 2_200 && rng.unit() < 0.75 {
                            truth[off - st] = b'T';
                        }
                    }
                }
            }
            let (mut s, mut q) = sim_r10_read(&truth, 0.01, &mut rng);
            if !is_fwd {
                (s, q) = revcomp_qual(&s, &q);
            }
            reads.push(Rec {
                id: format!("m{i}").into_bytes(),
                seq: s,
                qual: Some(q),
            });
        }
        let prof = detect_chemistry(&reads, OntChemistry::R10);
        let (out, _, rep) = run_r10(reads, prof).unwrap();
        assert!(
            rep.rounds[0].pruned_strand_bias_sites > 0,
            "expected strand-biased methylation sites to be pruned: {:?}",
            rep.rounds[0]
        );
        let total_bp: usize = out.iter().map(|r| r.seq.len()).sum();
        let total_err: usize = out.iter().map(|r| best_edit_distance(&r.seq, &g)).sum();
        let identity = 1.0 - total_err as f64 / total_bp as f64;
        assert!(
            identity > 0.9995,
            "expected >99.95% identity despite 75% +strand methylation skew, got {identity:.6} ({total_err} / {total_bp})"
        );
    }

    #[test]
    fn two_copy_diverged_repeat_is_phased_without_homogenization() {
        // 10 kb circular genome holding two 2.4 kb repeat copies (`1000..3400` and `6000..8400`)
        // that differ at 6 SNPs, 1 non-HP base insertion, and 1 homopolymer run length (`A6` vs `A3`).
        // Reads of 2.1 kb completely internal to the repeat must retain their own copy's alleles!
        let mut g = random_genome(10_000, 777);
        let mut copy1 = random_genome(2_400, 888);
        // Place clean `A6` homopolymer in copy1 at 900..906 flanked by C/G.
        copy1[899] = b'C';
        copy1[900..906].fill(b'A');
        copy1[906] = b'G';
        let mut copy2 = copy1.clone();
        for &snp in &[250usize, 550, 1_150, 1_450, 1_750, 2_050] {
            copy2[snp] = BASES[((base_code(copy2[snp]).unwrap() + 2) % 4) as usize];
        }
        // Shorten `A6` -> `A3` in copy2 and compensate length near end so copy2 is also 2,400 bp.
        copy2.drain(903..906);
        copy2.extend_from_slice(b"CTG");

        g[1_000..3_400].copy_from_slice(&copy1);
        g[6_000..8_400].copy_from_slice(&copy2);

        let circ = [g.as_slice(), &g[..2_500]].concat();
        let mut rng = Rng(999);
        let mut reads = Vec::new();
        for i in 0..250usize {
            let st = (i * 40) % 10_000;
            let (s, q) = sim_r10_read(&circ[st..st + 2_100], 0.012, &mut rng);
            let (s, q) = if i % 2 == 1 { revcomp_qual(&s, &q) } else { (s, q) };
            let tag = if (1_050..1_250).contains(&st) {
                "inner_copy1"
            } else if (6_050..6_250).contains(&st) {
                "inner_copy2"
            } else {
                "flank"
            };
            reads.push(Rec {
                id: format!("{tag}_{i}").into_bytes(),
                seq: s,
                qual: Some(q),
            });
        }

        let prof = detect_chemistry(&reads, OntChemistry::Auto);
        assert!(prof.use_r10);
        let (out, _, rep) = run_r10(reads, prof).unwrap();
        assert!(rep.rounds[0].phased_sites_kept > 0);

        let mut inner1_checked = 0usize;
        let mut inner2_checked = 0usize;
        for r in &out {
            let name = std::str::from_utf8(&r.id).unwrap();
            if name.starts_with("inner_copy1_") {
                let d1 = best_edit_distance(&r.seq, &copy1);
                let d2 = best_edit_distance(&r.seq, &copy2);
                assert!(d1 <= 1 && d2 >= 6, "{name}: homogenized! d1={d1}, d2={d2}");
                inner1_checked += 1;
            } else if name.starts_with("inner_copy2_") {
                let d1 = best_edit_distance(&r.seq, &copy1);
                let d2 = best_edit_distance(&r.seq, &copy2);
                assert!(d2 <= 1 && d1 >= 6, "{name}: homogenized! d1={d1}, d2={d2}");
                inner2_checked += 1;
            }
        }
        assert!(inner1_checked >= 4 && inner2_checked >= 4);
    }

    #[test]
    fn chimeric_reads_are_clipped_at_the_unbridged_junction() {
        let g = random_genome(6_000, 1234);
        let nuclear = random_genome(20_000, 5678);
        let circ = [g.as_slice(), &g[..2_500]].concat();
        let mut rng = Rng(42);
        let mut reads: Vec<Rec> = (0..120)
            .map(|i| {
                let st = (i * 50) % 6_000;
                let (s, q) = sim_r10_read(&circ[st..st + 2_200], 0.01, &mut rng);
                let (s, q) = if i % 2 == 1 { revcomp_qual(&s, &q) } else { (s, q) };
                Rec {
                    id: format!("clean_{i}").into_bytes(),
                    seq: s,
                    qual: Some(q),
                }
            })
            .collect();

        // Add 4 NUMT chimeric reads `[mito 2.0 kb | nuclear 2.0 kb]` and 2 artificial mito+mito
        // pore-ligation chimeras `[0..1800 | 3800..5600]`.
        for c in 0..4usize {
            let mut raw = g[500..2_500].to_vec();
            raw.extend_from_slice(&nuclear[c * 3_000..c * 3_000 + 2_000]);
            let (s, q) = sim_r10_read(&raw, 0.01, &mut rng);
            reads.push(Rec {
                id: format!("numt_chim_{c}").into_bytes(),
                seq: s,
                qual: Some(q),
            });
        }
        for c in 0..2usize {
            let mut raw = g[100..1_900].to_vec();
            raw.extend_from_slice(&g[3_700..5_500]);
            let (s, q) = sim_r10_read(&raw, 0.01, &mut rng);
            reads.push(Rec {
                id: format!("self_chim_{c}").into_bytes(),
                seq: s,
                qual: Some(q),
            });
        }

        let prof = detect_chemistry(&reads, OntChemistry::R10);
        let (out, _, rep) = run_r10(reads, prof).unwrap();
        assert!(
            rep.rounds[0].reads_chimera_clipped >= 6,
            "expected >=6 chimeras clipped: {:?}",
            rep.rounds[0]
        );
        for r in &out {
            let name = std::str::from_utf8(&r.id).unwrap();
            let err = best_edit_distance(&r.seq, &g);
            assert!(
                err <= 2,
                "{name} (len {}) still carries chimeric sequence: edit_dist={err}",
                r.seq.len()
            );
            if name.starts_with("clean_") {
                assert!(r.seq.len() >= 2_150, "clean read eroded: {}", r.seq.len());
            }
        }
    }
}
