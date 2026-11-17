//! Native pure-Rust BAM/modBAM (`MM`/`ML`) reader, plant-aware modification caller (`CG`, `CHG`,
//! `CHH` 5mC and `6mA`), adaptive NUMT/NUPT read separator, and repeat/IR-preserving organelle
//! methylation profiler (`ovasm check methylation` and `ovasm run ... --methylation`).

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;
use std::time::Instant;

use anyhow::{bail, ensure, Context, Result};
use clap::{Args, ValueEnum};
use flate2::read::MultiGzDecoder;
use flate2::{Decompress, FlushDecompress};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::Serialize;
use serde_json::{json, Value};

use crate::recruit::Rec;

pub const SCHEMA: &str = "ovasm.methylation.v1";
const CIRCULAR_WRAP_MAX: usize = 25_000;
const MAX_KMER_OCC: usize = 32;
const NUCLEAR_SAMPLE_MAX: usize = 25_000;

/// BAM nibble lookup (`=ACMGRSVTWYHKDBN`).
const BAM_SEQ_LUT: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

#[derive(Clone, Copy, Debug, Default)]
pub struct BamIngestFilter {
    pub max_cg_methylation: Option<f64>,
    pub min_cpgs: usize,
}

static ACTIVE_BAM_FILTER: Mutex<BamIngestFilter> = Mutex::new(BamIngestFilter {
    max_cg_methylation: None,
    min_cpgs: 5,
});
static BAM_FILTERED_READS: AtomicU64 = AtomicU64::new(0);
static BAM_FILTERED_BASES: AtomicU64 = AtomicU64::new(0);

/// Scope guard for recruitment-time BAM methylation filtering.
pub struct BamFilterGuard;

impl Drop for BamFilterGuard {
    fn drop(&mut self) {
        if let Ok(mut g) = ACTIVE_BAM_FILTER.lock() {
            *g = BamIngestFilter::default();
        }
    }
}

pub fn activate_bam_ingest_filter(max_cg_methylation: Option<f64>, min_cpgs: usize) -> BamFilterGuard {
    BAM_FILTERED_READS.store(0, Ordering::Relaxed);
    BAM_FILTERED_BASES.store(0, Ordering::Relaxed);
    if let Ok(mut g) = ACTIVE_BAM_FILTER.lock() {
        *g = BamIngestFilter {
            max_cg_methylation,
            min_cpgs: min_cpgs.max(1),
        };
    }
    BamFilterGuard
}

pub fn bam_ingest_filtered_counts() -> (u64, u64) {
    (
        BAM_FILTERED_READS.load(Ordering::Relaxed),
        BAM_FILTERED_BASES.load(Ordering::Relaxed),
    )
}

fn current_bam_filter() -> BamIngestFilter {
    ACTIVE_BAM_FILTER.lock().map(|g| *g).unwrap_or_default()
}

/// True when `path` has a `.bam` or `.ubam` extension.
pub fn is_bam_path(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("bam") || e.eq_ignore_ascii_case("ubam"))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum CircularMode {
    /// Treat contigs as circular unless FASTA header says `linear` or `circular=false`.
    Auto,
    /// Treat all target contigs as circular.
    All,
    /// Treat all target contigs as linear.
    None,
}

#[derive(Args, Clone, Debug)]
pub struct MethylationArgs {
    /// Target organelle FASTA (`molecules.fasta` or reference FASTA). Optional when `--gfa` is given.
    #[arg(long)]
    pub target: Option<PathBuf>,
    /// Optional organelle GFA (`organelle.gfa` or `assembly.gfa`) to annotate per-segment methylation tags.
    #[arg(long)]
    pub gfa: Option<PathBuf>,
    /// Input BAM/modBAM file(s) carrying `MM`/`ML` base-modification tags (aligned or unaligned).
    #[arg(short = 'b', long = "bam", required = true, num_args = 1..)]
    pub bam: Vec<PathBuf>,
    /// Output directory (`methylation.bed`, `reads.methylation.tsv`, `methylation.json`).
    #[arg(short = 'o', long = "out")]
    pub out: PathBuf,
    /// Minimum probability in [0, 1] to count a modified base call as methylated.
    #[arg(long, default_value_t = 0.5)]
    pub prob_threshold: f64,
    /// Max per-read CG 5mC fraction to keep a read as organelle-derived: `auto`, `none`, or float in [0, 1].
    #[arg(long, default_value = "auto")]
    pub max_cg_methylation: String,
    /// Minimum called CG sites on a read before CG methylation filtering applies.
    #[arg(long, default_value_t = 5)]
    pub min_read_cpgs: usize,
    /// Maximum unaligned terminal flank (bp) outside all target organelle sequences before calling chimeric NUMT/NUPT.
    #[arg(long, default_value_t = 300)]
    pub max_overhang: usize,
    /// Anchor k-mer length (15..=31) for colinear repeat-aware placement.
    #[arg(short = 'k', long, default_value_t = 31)]
    pub kmer: usize,
    /// Circular boundary wrap mode for target sequences.
    #[arg(long, value_enum, default_value_t = CircularMode::Auto)]
    pub circular: CircularMode,
    /// Lineage epigenetic model: `landplants` (`CG`/`CHG`/`CHH` + `6mA`), `metazoa` (`animal`:
    /// `CG`-only NUMT separation with vertebrate global-CpG vs invertebrate mosaic/unmethylated
    /// regime auto-detection and multi-lap ~16 kb circular chaining), or `fungi`.
    #[arg(short = 'x', long, value_enum, default_value_t = crate::seeddb::Lineage::Landplants)]
    pub lineage: crate::seeddb::Lineage,
    /// Worker threads for parallel BGZF decompression and read placement.
    #[arg(short = 't', long)]
    pub threads: Option<usize>,
}

// =============================================================================================
// Pure-Rust BGZF & BAM Streaming Decoder
// =============================================================================================

#[inline]
fn comp_base(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

fn revcomp_in_place(seq: &mut [u8]) {
    seq.reverse();
    for b in seq.iter_mut() {
        *b = comp_base(*b);
    }
}

/// Try to parse the total compressed block length (`BSIZE + 1`) and header length of one BGZF block.
fn parse_bgzf_header(buf: &[u8]) -> Option<(usize, usize)> {
    if buf.len() < 18 {
        return None;
    }
    if buf[0] != 0x1f || buf[1] != 0x8b || buf[2] != 0x08 || (buf[3] & 0x04) == 0 {
        return None;
    }
    let xlen = u16::from_le_bytes([buf[10], buf[11]]) as usize;
    if buf.len() < 12 + xlen + 8 {
        return None;
    }
    let mut off = 12;
    let end = 12 + xlen;
    let mut bsize: Option<usize> = None;
    while off + 4 <= end {
        let si1 = buf[off];
        let si2 = buf[off + 1];
        let slen = u16::from_le_bytes([buf[off + 2], buf[off + 3]]) as usize;
        off += 4;
        if off + slen > end {
            return None;
        }
        if si1 == b'B' && si2 == b'C' && slen == 2 {
            bsize = Some(u16::from_le_bytes([buf[off], buf[off + 1]]) as usize + 1);
        }
        off += slen;
    }
    let total = bsize?;
    let mut hdr_end = end;
    let flg = buf[3];
    if flg & 0x08 != 0 {
        while hdr_end < buf.len() && buf[hdr_end] != 0 {
            hdr_end += 1;
        }
        hdr_end += 1;
    }
    if flg & 0x10 != 0 {
        while hdr_end < buf.len() && buf[hdr_end] != 0 {
            hdr_end += 1;
        }
        hdr_end += 1;
    }
    if flg & 0x02 != 0 {
        hdr_end += 2;
    }
    if total < hdr_end + 8 {
        return None;
    }
    Some((total, hdr_end))
}

fn inflate_bgzf_block(block: &[u8], hdr_len: usize) -> Result<Vec<u8>> {
    let total = block.len();
    let isize_uncomp = u32::from_le_bytes([
        block[total - 4],
        block[total - 3],
        block[total - 2],
        block[total - 1],
    ]) as usize;
    if isize_uncomp == 0 {
        return Ok(Vec::new());
    }
    let mut out = vec![0u8; isize_uncomp];
    let mut dec = Decompress::new(false);
    let payload = &block[hdr_len..total - 8];
    let status = dec
        .decompress(payload, &mut out, FlushDecompress::Finish)
        .context("corrupt BGZF deflate block")?;
    if status != flate2::Status::StreamEnd && dec.total_out() as usize != isize_uncomp {
        bail!("incomplete BGZF block decompression");
    }
    out.truncate(dec.total_out() as usize);
    Ok(out)
}

/// Decompress a BAM stream in parallel BGZF batches (with automatic fallback to `MultiGzDecoder`).
struct BamDecompressor {
    file: File,
    comp_buf: Vec<u8>,
    comp_len: usize,
    eof: bool,
    fallback: Option<MultiGzDecoder<BufReader<File>>>,
    path: PathBuf,
}

impl BamDecompressor {
    fn open(path: &Path) -> Result<Self> {
        let mut f = File::open(path).with_context(|| format!("cannot open BAM {}", path.display()))?;
        let mut probe = [0u8; 64];
        let n = f.read(&mut probe)?;
        let f = File::open(path).with_context(|| format!("cannot reopen BAM {}", path.display()))?;
        if n >= 28 && parse_bgzf_header(&probe[..n]).is_some() {
            Ok(Self {
                file: f,
                comp_buf: vec![0u8; 16 << 20],
                comp_len: 0,
                eof: false,
                fallback: None,
                path: path.to_path_buf(),
            })
        } else {
            Ok(Self {
                file: File::open(path)?,
                comp_buf: Vec::new(),
                comp_len: 0,
                eof: false,
                fallback: Some(MultiGzDecoder::new(BufReader::with_capacity(1 << 20, f))),
                path: path.to_path_buf(),
            })
        }
    }

    fn next_uncompressed_chunk(&mut self, out: &mut Vec<u8>) -> Result<bool> {
        if let Some(gz) = &mut self.fallback {
            let mut tmp = [0u8; 256 * 1024];
            let n = gz
                .read(&mut tmp)
                .with_context(|| format!("gzip read error in {}", self.path.display()))?;
            if n == 0 {
                return Ok(false);
            }
            out.extend_from_slice(&tmp[..n]);
            return Ok(true);
        }
        while !self.eof && self.comp_len < self.comp_buf.len() {
            let n = self
                .file
                .read(&mut self.comp_buf[self.comp_len..])
                .with_context(|| format!("read error in {}", self.path.display()))?;
            if n == 0 {
                self.eof = true;
                break;
            }
            self.comp_len += n;
            if self.comp_len >= (8 << 20) {
                break;
            }
        }
        if self.comp_len == 0 {
            return Ok(false);
        }
        let mut specs: Vec<(usize, usize, usize)> = Vec::with_capacity(256);
        let mut pos = 0usize;
        while pos < self.comp_len {
            if self.comp_len - pos < 28 && !self.eof {
                break;
            }
            match parse_bgzf_header(&self.comp_buf[pos..self.comp_len]) {
                Some((total, hdr_len)) => {
                    if pos + total > self.comp_len {
                        if self.eof {
                            bail!("truncated BGZF block at end of {}", self.path.display());
                        }
                        break;
                    }
                    specs.push((pos, total, hdr_len));
                    pos += total;
                }
                None => {
                    bail!("invalid BGZF block header in {}", self.path.display());
                }
            }
        }
        if specs.is_empty() {
            return Ok(false);
        }
        let slice = &self.comp_buf[..pos];
        let blocks: Result<Vec<Vec<u8>>> = specs
            .par_iter()
            .map(|&(off, total, hdr_len)| inflate_bgzf_block(&slice[off..off + total], hdr_len))
            .collect();
        let blocks = blocks?;
        let add_bytes: usize = blocks.iter().map(Vec::len).sum();
        out.reserve(add_bytes);
        for b in blocks {
            out.extend_from_slice(&b);
        }
        self.comp_buf.copy_within(pos..self.comp_len, 0);
        self.comp_len -= pos;
        Ok(true)
    }
}

/// Decoded primary BAM record in original sequenced 5'->3' orientation.
#[derive(Clone, Debug)]
pub struct BamRead {
    pub name: Vec<u8>,
    #[allow(dead_code)]
    pub flag: u16,
    /// Original sequenced 5'->3' bases (already reverse-complemented back if `flag & 0x10 != 0`).
    pub seq: Vec<u8>,
    /// Phred+33 ASCII qualities in original sequenced 5'->3' order (`None` if missing `0xff`).
    pub qual: Option<Vec<u8>>,
    /// Raw `MM`/`Mm` tag string bytes (without trailing NUL), if present.
    pub mm: Option<Vec<u8>>,
    /// Raw `ML`/`Ml` `B:C` probability array (`0..=255`), if present.
    pub ml: Option<Vec<u8>>,
}

fn parse_aux_mm_ml(mut aux: &[u8]) -> (Option<Vec<u8>>, Option<Vec<u8>>) {
    let mut mm = None;
    let mut ml = None;
    while aux.len() >= 3 {
        let t0 = aux[0];
        let t1 = aux[1];
        let vtype = aux[2];
        aux = &aux[3..];
        let is_mm = t0 == b'M' && (t1 == b'M' || t1 == b'm');
        let is_ml = t0 == b'M' && (t1 == b'L' || t1 == b'l');
        match vtype {
            b'A' | b'c' | b'C' => {
                if aux.is_empty() {
                    break;
                }
                aux = &aux[1..];
            }
            b's' | b'S' => {
                if aux.len() < 2 {
                    break;
                }
                aux = &aux[2..];
            }
            b'i' | b'I' | b'f' => {
                if aux.len() < 4 {
                    break;
                }
                aux = &aux[4..];
            }
            b'd' => {
                if aux.len() < 8 {
                    break;
                }
                aux = &aux[8..];
            }
            b'Z' | b'H' => {
                let Some(nul) = aux.iter().position(|&b| b == 0) else {
                    break;
                };
                if is_mm && vtype == b'Z' {
                    mm = Some(aux[..nul].to_vec());
                }
                aux = &aux[nul + 1..];
            }
            b'B' => {
                if aux.len() < 5 {
                    break;
                }
                let sub = aux[0];
                let cnt = u32::from_le_bytes([aux[1], aux[2], aux[3], aux[4]]) as usize;
                aux = &aux[5..];
                let elem = match sub {
                    b'c' | b'C' => 1usize,
                    b's' | b'S' => 2usize,
                    b'i' | b'I' | b'f' => 4usize,
                    _ => break,
                };
                let bytes = cnt.saturating_mul(elem);
                if aux.len() < bytes {
                    break;
                }
                if is_ml && sub == b'C' {
                    ml = Some(aux[..bytes].to_vec());
                }
                aux = &aux[bytes..];
            }
            _ => break,
        }
        if mm.is_some() && ml.is_some() {
            break;
        }
    }
    (mm, ml)
}

/// Stream all primary non-QC-fail records from `bam_path`.
pub fn for_each_bam_read(
    bam_path: &Path,
    need_mods: bool,
    mut cb: impl FnMut(BamRead) -> Result<bool>,
) -> Result<()> {
    let mut dec = BamDecompressor::open(bam_path)?;
    let mut buf: Vec<u8> = Vec::with_capacity(32 << 20);
    let mut cur = 0usize;

    let ensure_bytes = |dec: &mut BamDecompressor, buf: &mut Vec<u8>, cur: &mut usize, need: usize| -> Result<bool> {
        if *cur > (16 << 20) {
            buf.drain(..*cur);
            *cur = 0;
        }
        while buf.len() - *cur < need {
            if !dec.next_uncompressed_chunk(buf)? {
                return Ok(false);
            }
        }
        Ok(true)
    };

    if !ensure_bytes(&mut dec, &mut buf, &mut cur, 8)? {
        bail!("empty BAM file: {}", bam_path.display());
    }
    ensure!(
        &buf[cur..cur + 4] == b"BAM\x01",
        "{} is not a BAM file (missing BAM\\1 magic)",
        bam_path.display()
    );
    let l_text = i32::from_le_bytes([buf[cur + 4], buf[cur + 5], buf[cur + 6], buf[cur + 7]]) as usize;
    cur += 8;
    ensure!(
        ensure_bytes(&mut dec, &mut buf, &mut cur, l_text + 4)?,
        "truncated BAM header text in {}",
        bam_path.display()
    );
    cur += l_text;
    let n_ref = i32::from_le_bytes([buf[cur], buf[cur + 1], buf[cur + 2], buf[cur + 3]]) as usize;
    cur += 4;
    for _ in 0..n_ref {
        ensure!(
            ensure_bytes(&mut dec, &mut buf, &mut cur, 4)?,
            "truncated BAM ref list in {}",
            bam_path.display()
        );
        let l_name = i32::from_le_bytes([buf[cur], buf[cur + 1], buf[cur + 2], buf[cur + 3]]) as usize;
        cur += 4;
        ensure!(
            ensure_bytes(&mut dec, &mut buf, &mut cur, l_name + 4)?,
            "truncated BAM ref entry in {}",
            bam_path.display()
        );
        cur += l_name + 4;
    }

    loop {
        if !ensure_bytes(&mut dec, &mut buf, &mut cur, 4)? {
            ensure!(
                buf.len() == cur,
                "truncated BAM record size at EOF in {}",
                bam_path.display()
            );
            break;
        }
        let block_size = u32::from_le_bytes([buf[cur], buf[cur + 1], buf[cur + 2], buf[cur + 3]]) as usize;
        cur += 4;
        ensure!(block_size >= 32, "malformed BAM block_size {block_size}");
        ensure!(
            ensure_bytes(&mut dec, &mut buf, &mut cur, block_size)?,
            "truncated BAM record payload in {}",
            bam_path.display()
        );
        let rec = &buf[cur..cur + block_size];
        cur += block_size;

        let l_read_name = rec[8] as usize;
        let n_cigar_op = u16::from_le_bytes([rec[12], rec[13]]) as usize;
        let flag = u16::from_le_bytes([rec[14], rec[15]]);
        if flag & 0x0f00 != 0 {
            continue;
        }
        let l_seq = i32::from_le_bytes([rec[16], rec[17], rec[18], rec[19]]) as usize;
        if l_seq == 0 {
            continue;
        }
        let name_start = 32usize;
        let cigar_start = name_start + l_read_name;
        let seq_bytes = (l_seq + 1) / 2;
        let seq_start = cigar_start + n_cigar_op * 4;
        let qual_start = seq_start + seq_bytes;
        let aux_start = qual_start + l_seq;
        if aux_start > block_size {
            bail!("corrupt BAM record offsets in {}", bam_path.display());
        }
        let raw_name = &rec[name_start..cigar_start];
        let name_end = raw_name
            .iter()
            .position(|&b| b == 0)
            .unwrap_or(raw_name.len());
        let name = raw_name[..name_end].to_vec();

        let mut seq = Vec::with_capacity(l_seq);
        let packed = &rec[seq_start..qual_start];
        for i in 0..l_seq {
            let byte = packed[i >> 1];
            let nib = if i & 1 == 0 { byte >> 4 } else { byte & 0x0f } as usize;
            seq.push(BAM_SEQ_LUT[nib]);
        }
        let raw_qual = &rec[qual_start..aux_start];
        let qual = if raw_qual.first() == Some(&0xff) {
            None
        } else {
            let mut q: Vec<u8> = raw_qual.iter().map(|&b| b.saturating_add(33).min(126)).collect();
            if flag & 0x10 != 0 {
                q.reverse();
            }
            Some(q)
        };
        if flag & 0x10 != 0 {
            revcomp_in_place(&mut seq);
        }
        let (mm, ml) = if need_mods {
            parse_aux_mm_ml(&rec[aux_start..block_size])
        } else {
            (None, None)
        };
        if !cb(BamRead {
            name,
            flag,
            seq,
            qual,
            mm,
            ml,
        })? {
            break;
        }
    }
    Ok(())
}

/// First `limit` primary read lengths in a BAM file (for `assemble::median_read_length`).
pub fn first_bam_read_lengths(path: &Path, limit: usize) -> Result<Vec<usize>> {
    let mut lens = Vec::with_capacity(limit.min(10_000));
    for_each_bam_read(path, false, |r| {
        lens.push(r.seq.len());
        Ok(lens.len() < limit)
    })?;
    Ok(lens)
}

/// Pre-scan up to `sample_limit` primary BAM records carrying `MM`/`ML` CpG calls to infer
/// a lineage-aware pre-recruitment `--max-read-methylation` cutoff when `--methylation` is
/// enabled without an explicit `--max-read-methylation` value.
///
/// Returns `Some((max_cg, regime_name))` when hypermethylated nuclear reads are detected so
/// nuclear NUMTs/NUPTs are filtered before recruitment and graph construction, or `None` when
/// BAMs lack `MM`/`ML` tags or belong to an unmethylated clade (e.g. Diptera, Nematoda) where
/// pre-recruitment CpG pruning should stay off to avoid false read loss.
pub fn auto_ingest_max_cg_methylation(
    paths: &[PathBuf],
    lineage: crate::seeddb::Lineage,
) -> Option<(f64, &'static str)> {
    const SAMPLE_LIMIT: usize = 2000;
    const MIN_CPGS: usize = 5;
    let mut cgs: Vec<f64> = Vec::new();
    for path in paths.iter().filter(|p| is_bam_path(p)) {
        if cgs.len() >= SAMPLE_LIMIT {
            break;
        }
        let _ = for_each_bam_read(path, true, |r| {
            if let (Some(mm), Some(ml)) = (r.mm.as_deref(), r.ml.as_deref()) {
                let s = summarize_read_mods(&r.seq, mm, ml, 0.5);
                if s.cg_calls >= MIN_CPGS {
                    cgs.push(s.cg_frac());
                }
            }
            Ok(cgs.len() < SAMPLE_LIMIT)
        });
    }
    if cgs.len() < 16 {
        return None;
    }
    cgs.sort_by(|a, b| a.total_cmp(b));
    let med = quantile_sorted(&cgs, 0.50);
    let p90 = quantile_sorted(&cgs, 0.90);
    let high: Vec<f64> = cgs.iter().copied().filter(|&x| x >= 0.40).collect();
    let high_share = high.len() as f64 / cgs.len() as f64;
    let high_med = quantile_sorted(&high, 0.50);
    match lineage {
        crate::seeddb::Lineage::Landplants => {
            if high_share >= 0.02 || p90 >= 0.35 {
                Some((0.40, EpigeneticRegime::PlantCgChgRddm.as_str()))
            } else {
                None
            }
        }
        crate::seeddb::Lineage::Metazoa => {
            if med >= 0.40 || (high_share >= 0.04 && high_med >= 0.50) {
                Some((0.30, EpigeneticRegime::VertebrateGlobalCpg.as_str()))
            } else if high_share >= 0.01 || p90 >= 0.25 {
                Some((0.38, EpigeneticRegime::InvertebrateMosaicCpg.as_str()))
            } else {
                // DNMT-deficient / unmethylated animal clade (e.g. Diptera, Nematoda):
                // keep all reads during recruitment and let `check methylation` profile cleanly.
                None
            }
        }
        crate::seeddb::Lineage::Fungi => {
            if high_share >= 0.02 || p90 >= 0.35 {
                Some((0.40, EpigeneticRegime::FungiVariable5mC.as_str()))
            } else {
                None
            }
        }
    }
}

/// Convert BAM primary records into `recruit::Rec` records, applying any active BAM methylation filter.
pub fn for_each_bam_rec(path: &Path, mut f: impl FnMut(Rec)) -> Result<()> {
    let filt = current_bam_filter();
    let need_mods = filt.max_cg_methylation.is_some();
    for_each_bam_read(path, need_mods, |r| {
        if let Some(max_cg) = filt.max_cg_methylation {
            if let (Some(mm), Some(ml)) = (r.mm.as_deref(), r.ml.as_deref()) {
                let summary = summarize_read_mods(&r.seq, mm, ml, 0.5);
                if summary.cg_calls >= filt.min_cpgs && summary.cg_frac() > max_cg {
                    BAM_FILTERED_READS.fetch_add(1, Ordering::Relaxed);
                    BAM_FILTERED_BASES.fetch_add(r.seq.len() as u64, Ordering::Relaxed);
                    return Ok(true);
                }
            }
        }
        f(Rec {
            id: r.name,
            seq: r.seq,
            qual: r.qual,
        });
        Ok(true)
    })
}

/// Stream BAM records in `BATCH_BASES` batches into `tx` for `recruit::read_file`.
pub fn stream_bam_file(
    path: &Path,
    batch_bases: usize,
    tx: &crossbeam_channel::Sender<Vec<Rec>>,
) -> Result<()> {
    let filt = current_bam_filter();
    let need_mods = filt.max_cg_methylation.is_some();
    let (mut batch, mut bases) = (Vec::new(), 0usize);
    for_each_bam_read(path, need_mods, |r| {
        if let Some(max_cg) = filt.max_cg_methylation {
            if let (Some(mm), Some(ml)) = (r.mm.as_deref(), r.ml.as_deref()) {
                let summary = summarize_read_mods(&r.seq, mm, ml, 0.5);
                if summary.cg_calls >= filt.min_cpgs && summary.cg_frac() > max_cg {
                    BAM_FILTERED_READS.fetch_add(1, Ordering::Relaxed);
                    BAM_FILTERED_BASES.fetch_add(r.seq.len() as u64, Ordering::Relaxed);
                    return Ok(true);
                }
            }
        }
        bases += r.seq.len();
        batch.push(Rec {
            id: r.name,
            seq: r.seq,
            qual: r.qual,
        });
        if bases >= batch_bases {
            if tx.send(std::mem::take(&mut batch)).is_err() {
                return Ok(false);
            }
            bases = 0;
        }
        Ok(true)
    })?;
    if !batch.is_empty() {
        let _ = tx.send(batch);
    }
    Ok(())
}

// =============================================================================================
// SAMtags §1.7 MM/ML Decoder & Plant Context Classifier (CG / CHG / CHH / 6mA)
// =============================================================================================

#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ModContext {
    Cg = 0,
    Chg = 1,
    Chh = 2,
    M6a = 3,
}

impl ModContext {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Cg => "CG",
            Self::Chg => "CHG",
            Self::Chh => "CHH",
            Self::M6a => "6mA",
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ReadModCall {
    pub read_pos: u32,
    pub on_orig_strand: bool,
    pub context: ModContext,
    pub qual: u8,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct ReadModSummary {
    pub cg_calls: usize,
    pub cg_mod: usize,
    pub chg_calls: usize,
    pub chg_mod: usize,
    pub chh_calls: usize,
    pub chh_mod: usize,
    pub m6a_calls: usize,
    pub m6a_mod: usize,
}

impl ReadModSummary {
    pub fn cg_frac(&self) -> f64 {
        if self.cg_calls == 0 {
            0.0
        } else {
            self.cg_mod as f64 / self.cg_calls as f64
        }
    }
    pub fn chg_frac(&self) -> f64 {
        if self.chg_calls == 0 {
            0.0
        } else {
            self.chg_mod as f64 / self.chg_calls as f64
        }
    }
    pub fn chh_frac(&self) -> f64 {
        if self.chh_calls == 0 {
            0.0
        } else {
            self.chh_mod as f64 / self.chh_calls as f64
        }
    }
    pub fn m6a_frac(&self) -> f64 {
        if self.m6a_calls == 0 {
            0.0
        } else {
            self.m6a_mod as f64 / self.m6a_calls as f64
        }
    }
}

#[inline]
fn is_h_base(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'T' | b'a' | b'c' | b't')
}

#[inline]
fn classify_c_context_plus(seq: &[u8], pos: usize) -> Option<ModContext> {
    if pos + 1 >= seq.len() {
        return None;
    }
    let b1 = seq[pos + 1].to_ascii_uppercase();
    if b1 == b'G' {
        return Some(ModContext::Cg);
    }
    if is_h_base(b1) && pos + 2 < seq.len() {
        let b2 = seq[pos + 2].to_ascii_uppercase();
        if b2 == b'G' {
            return Some(ModContext::Chg);
        }
        if is_h_base(b2) {
            return Some(ModContext::Chh);
        }
    }
    None
}

#[inline]
fn classify_c_context_minus(seq: &[u8], pos: usize) -> Option<ModContext> {
    if pos == 0 {
        return None;
    }
    let b1 = comp_base(seq[pos - 1]);
    if b1 == b'G' {
        return Some(ModContext::Cg);
    }
    if is_h_base(b1) && pos >= 2 {
        let b2 = comp_base(seq[pos - 2]);
        if b2 == b'G' {
            return Some(ModContext::Chg);
        }
        if is_h_base(b2) {
            return Some(ModContext::Chh);
        }
    }
    None
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum ParsedCode {
    Cytosine5mC,
    Cytosine5hmC,
    Adenine6mA,
    Other,
}

fn parse_mod_codes(raw: &[u8]) -> (Vec<ParsedCode>, Option<u8>) {
    if raw.is_empty() {
        return (Vec::new(), None);
    }
    let (core, mode) = match raw.last().copied() {
        Some(b'?') | Some(b'.') => (&raw[..raw.len() - 1], raw.last().copied()),
        _ => (raw, None),
    };
    if core.is_empty() {
        return (Vec::new(), mode);
    }
    if core.iter().all(|b| b.is_ascii_digit()) {
        let code = match core {
            b"21839" => ParsedCode::Cytosine5mC,
            b"76792" => ParsedCode::Cytosine5hmC,
            b"28871" => ParsedCode::Adenine6mA,
            _ => ParsedCode::Other,
        };
        return (vec![code], mode);
    }
    let codes = core
        .iter()
        .map(|&c| match c {
            b'm' => ParsedCode::Cytosine5mC,
            b'h' => ParsedCode::Cytosine5hmC,
            b'a' => ParsedCode::Adenine6mA,
            _ => ParsedCode::Other,
        })
        .collect();
    (codes, mode)
}

/// Decode all `5mC`/`5hmC` (`CG`/`CHG`/`CHH`) and `6mA` calls on `orig_seq` (5'->3').
pub fn decode_read_mods(orig_seq: &[u8], mm: &[u8], ml: &[u8]) -> Vec<ReadModCall> {
    if mm.is_empty() || orig_seq.is_empty() {
        return Vec::new();
    }
    let mut ml_idx = 0usize;
    let mut c_plus: FxHashMap<u32, u8> = FxHashMap::default();
    let mut c_minus: FxHashMap<u32, u8> = FxHashMap::default();
    let mut a_plus: FxHashMap<u32, u8> = FxHashMap::default();
    let mut a_minus: FxHashMap<u32, u8> = FxHashMap::default();
    let mut c_plus_implicit_zero = false;
    let mut a_plus_implicit_zero = false;

    for series in mm.split(|&b| b == b';') {
        if series.is_empty() {
            continue;
        }
        let mut parts = series.split(|&b| b == b',');
        let Some(header) = parts.next() else {
            continue;
        };
        if header.len() < 3 {
            continue;
        }
        let fund_base = header[0].to_ascii_uppercase();
        let strand = header[1];
        if strand != b'+' && strand != b'-' {
            continue;
        }
        let (codes, mode) = parse_mod_codes(&header[2..]);
        let n_codes = codes.len();
        if n_codes == 0 {
            continue;
        }
        let deltas: Vec<usize> = parts
            .filter(|s| !s.is_empty())
            .filter_map(|s| std::str::from_utf8(s).ok()?.parse::<usize>().ok())
            .collect();

        let has_5mc_or_hmc = codes
            .iter()
            .any(|&c| matches!(c, ParsedCode::Cytosine5mC | ParsedCode::Cytosine5hmC));
        let has_6ma = codes.contains(&ParsedCode::Adenine6mA);

        if strand == b'+' && mode == Some(b'.') {
            if fund_base == b'C' && has_5mc_or_hmc {
                c_plus_implicit_zero = true;
            }
            if fund_base == b'A' && has_6ma {
                a_plus_implicit_zero = true;
            }
        }

        if !has_5mc_or_hmc && !has_6ma {
            ml_idx = ml_idx.saturating_add(deltas.len() * n_codes);
            continue;
        }

        let target_base = if strand == b'+' {
            fund_base
        } else {
            comp_base(fund_base)
        };
        let mut positions: Vec<u32> = Vec::new();
        if strand == b'+' {
            for (i, &b) in orig_seq.iter().enumerate() {
                if fund_base == b'N' || b.to_ascii_uppercase() == target_base {
                    positions.push(i as u32);
                }
            }
        } else {
            for (i, &b) in orig_seq.iter().enumerate().rev() {
                if fund_base == b'N' || b.to_ascii_uppercase() == target_base {
                    positions.push(i as u32);
                }
            }
        }

        let mut occ_cursor = 0usize;
        for skip in deltas {
            occ_cursor += skip;
            if occ_cursor < positions.len() && ml_idx + n_codes <= ml.len() {
                let pos = positions[occ_cursor];
                for (ci, &code) in codes.iter().enumerate() {
                    let q = ml[ml_idx + ci];
                    match code {
                        ParsedCode::Cytosine5mC | ParsedCode::Cytosine5hmC => {
                            let map = if strand == b'+' { &mut c_plus } else { &mut c_minus };
                            let entry = map.entry(pos).or_insert(0);
                            *entry = entry.saturating_add(q);
                        }
                        ParsedCode::Adenine6mA => {
                            let map = if strand == b'+' { &mut a_plus } else { &mut a_minus };
                            let entry = map.entry(pos).or_insert(0);
                            *entry = (*entry).max(q);
                        }
                        ParsedCode::Other => {}
                    }
                }
            }
            ml_idx += n_codes;
            occ_cursor += 1;
        }
    }

    let mut out = Vec::with_capacity(c_plus.len() + c_minus.len() + a_plus.len());

    if c_plus_implicit_zero {
        for i in 0..orig_seq.len() {
            if orig_seq[i].eq_ignore_ascii_case(&b'C') {
                if let Some(ctx) = classify_c_context_plus(orig_seq, i) {
                    let q = c_plus.get(&(i as u32)).copied().unwrap_or(0);
                    out.push(ReadModCall {
                        read_pos: i as u32,
                        on_orig_strand: true,
                        context: ctx,
                        qual: q,
                    });
                }
            }
        }
    } else {
        let mut keys: Vec<u32> = c_plus.keys().copied().collect();
        keys.sort_unstable();
        for pos in keys {
            if let Some(ctx) = classify_c_context_plus(orig_seq, pos as usize) {
                out.push(ReadModCall {
                    read_pos: pos,
                    on_orig_strand: true,
                    context: ctx,
                    qual: c_plus[&pos],
                });
            }
        }
    }

    if !c_minus.is_empty() {
        let mut keys: Vec<u32> = c_minus.keys().copied().collect();
        keys.sort_unstable();
        for pos in keys {
            if let Some(ctx) = classify_c_context_minus(orig_seq, pos as usize) {
                out.push(ReadModCall {
                    read_pos: pos,
                    on_orig_strand: false,
                    context: ctx,
                    qual: c_minus[&pos],
                });
            }
        }
    }

    if a_plus_implicit_zero {
        for i in 0..orig_seq.len() {
            if orig_seq[i].eq_ignore_ascii_case(&b'A') {
                let q = a_plus.get(&(i as u32)).copied().unwrap_or(0);
                out.push(ReadModCall {
                    read_pos: i as u32,
                    on_orig_strand: true,
                    context: ModContext::M6a,
                    qual: q,
                });
            }
        }
    } else if !a_plus.is_empty() {
        let mut keys: Vec<u32> = a_plus.keys().copied().collect();
        keys.sort_unstable();
        for pos in keys {
            out.push(ReadModCall {
                read_pos: pos,
                on_orig_strand: true,
                context: ModContext::M6a,
                qual: a_plus[&pos],
            });
        }
    }

    if !a_minus.is_empty() {
        let mut keys: Vec<u32> = a_minus.keys().copied().collect();
        keys.sort_unstable();
        for pos in keys {
            out.push(ReadModCall {
                read_pos: pos,
                on_orig_strand: false,
                context: ModContext::M6a,
                qual: a_minus[&pos],
            });
        }
    }

    out
}

pub fn summarize_calls(calls: &[ReadModCall], prob_threshold: f64) -> ReadModSummary {
    let thr = (prob_threshold.clamp(0.0, 1.0) * 255.0).round() as u8;
    let mut s = ReadModSummary::default();
    for c in calls {
        let is_mod = usize::from(c.qual >= thr && c.qual > 0);
        match c.context {
            ModContext::Cg => {
                s.cg_calls += 1;
                s.cg_mod += is_mod;
            }
            ModContext::Chg => {
                s.chg_calls += 1;
                s.chg_mod += is_mod;
            }
            ModContext::Chh => {
                s.chh_calls += 1;
                s.chh_mod += is_mod;
            }
            ModContext::M6a => {
                s.m6a_calls += 1;
                s.m6a_mod += is_mod;
            }
        }
    }
    s
}

pub fn summarize_read_mods(orig_seq: &[u8], mm: &[u8], ml: &[u8], prob_threshold: f64) -> ReadModSummary {
    let calls = decode_read_mods(orig_seq, mm, ml);
    summarize_calls(&calls, prob_threshold)
}

// =============================================================================================
// Repeat/IR-Aware Colinear Anchor Index & Base-to-Base Alignment
// =============================================================================================

#[derive(Clone, Debug)]
struct TargetContig {
    name: String,
    seq: Vec<u8>,
    ext_seq: Vec<u8>,
    circular: bool,
    true_cpg_dyads: usize,
    true_chg_sites: usize,
    true_chh_sites: usize,
    true_a_sites: usize,
}

#[derive(Clone, Copy, Debug)]
struct KmerHit {
    contig: u16,
    pos: u32,
    is_rc: bool,
}

struct TargetIndex {
    k: usize,
    contigs: Vec<TargetContig>,
    bloom: Vec<u64>,
    bloom_mask: usize,
    map: FxHashMap<u64, Vec<KmerHit>>,
}

#[inline]
fn base_2bit(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

#[inline]
fn hash64(mut x: u64) -> u64 {
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}

fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    hash64(h)
}

fn scan_canonical_kmers(seq: &[u8], k: usize, step: usize, mut f: impl FnMut(usize, u64, bool)) {
    if k == 0 || k > 31 || seq.len() < k {
        return;
    }
    let mask: u64 = (1u64 << (2 * k)) - 1;
    let rc_shift = 2 * (k - 1);
    let mut fwd: u64 = 0;
    let mut rev: u64 = 0;
    let mut valid = 0usize;
    let step = step.max(1);
    for (i, &b) in seq.iter().enumerate() {
        if let Some(v) = base_2bit(b) {
            fwd = ((fwd << 2) | v) & mask;
            rev = (rev >> 2) | ((3 - v) << rc_shift);
            valid += 1;
            if valid >= k {
                let start = i + 1 - k;
                if start % step == 0 || i + 1 == seq.len() {
                    if fwd <= rev {
                        f(start, fwd, false);
                    } else {
                        f(start, rev, true);
                    }
                }
            }
        } else {
            valid = 0;
            fwd = 0;
            rev = 0;
        }
    }
}

#[inline]
fn contig_c_context(seq: &[u8], circular: bool, pos: usize, plus_strand: bool) -> Option<ModContext> {
    let n = seq.len();
    if pos >= n {
        return None;
    }
    if plus_strand {
        if seq[pos] != b'C' {
            return None;
        }
        let b1 = if pos + 1 < n {
            seq[pos + 1]
        } else if circular {
            seq[(pos + 1) % n]
        } else {
            return None;
        };
        if b1 == b'G' {
            return Some(ModContext::Cg);
        }
        if is_h_base(b1) {
            let b2 = if pos + 2 < n {
                seq[pos + 2]
            } else if circular {
                seq[(pos + 2) % n]
            } else {
                return None;
            };
            if b2 == b'G' {
                return Some(ModContext::Chg);
            }
            if is_h_base(b2) {
                return Some(ModContext::Chh);
            }
        }
        None
    } else {
        if seq[pos] != b'G' {
            return None;
        }
        let p1 = if pos >= 1 {
            pos - 1
        } else if circular {
            n - 1
        } else {
            return None;
        };
        let b1 = comp_base(seq[p1]);
        if b1 == b'G' {
            return Some(ModContext::Cg);
        }
        if is_h_base(b1) {
            let p2 = if pos >= 2 {
                pos - 2
            } else if circular && n >= 2 {
                (pos + n) - 2
            } else {
                return None;
            };
            let b2 = comp_base(seq[p2]);
            if b2 == b'G' {
                return Some(ModContext::Chg);
            }
            if is_h_base(b2) {
                return Some(ModContext::Chh);
            }
        }
        None
    }
}

fn count_contig_motifs(seq: &[u8], circular: bool) -> (usize, usize, usize, usize) {
    let n = seq.len();
    if n == 0 {
        return (0, 0, 0, 0);
    }
    let (mut cpg, mut chg, mut chh, mut a_sites) = (0usize, 0usize, 0usize, 0usize);
    for i in 0..n {
        let b0 = seq[i];
        if b0 == b'A' || b0 == b'T' {
            a_sites += 1;
        }
        match contig_c_context(seq, circular, i, true) {
            Some(ModContext::Cg) => cpg += 1,
            Some(ModContext::Chg) => chg += 1,
            Some(ModContext::Chh) => chh += 1,
            _ => {}
        }
        match contig_c_context(seq, circular, i, false) {
            Some(ModContext::Chg) => chg += 1,
            Some(ModContext::Chh) => chh += 1,
            _ => {}
        }
    }
    (cpg, chg, chh, a_sites)
}

fn build_circular_ext_seq(seq: &[u8], circular: bool) -> Vec<u8> {
    if !circular || seq.len() <= 64 {
        return seq.to_vec();
    }
    // Large plant organelles (120-2,000 kb) only need one 25 kb wrap across the origin, whereas
    // compact animal/fungal mitogenomes (~14-20 kb) are frequently spanned >1x by single 20-50 kb
    // ONT/HiFi rolling-circle or origin-crossing reads. Wrapping up to `seq.len() + 2 * 25 kb`
    // (capped at 6 laps, well below MAX_KMER_OCC = 32) prevents multi-lap animal mtDNA reads from
    // hitting the end of `ext_seq` and being misclassified as chimeric NUMTs.
    let target_len = (seq.len() + 2 * CIRCULAR_WRAP_MAX).min(seq.len() * 6);
    let mut ext = Vec::with_capacity(target_len);
    ext.extend_from_slice(seq);
    while ext.len() < target_len {
        let take = (target_len - ext.len()).min(seq.len());
        ext.extend_from_slice(&seq[..take]);
    }
    ext
}

fn load_targets_from_fasta(path: &Path, mode: CircularMode) -> Result<Vec<TargetContig>> {
    let mut reader = needletail::parse_fastx_file(path)
        .with_context(|| format!("cannot read target FASTA {}", path.display()))?;
    let mut out = Vec::new();
    while let Some(rec) = reader.next() {
        let rec = rec?;
        let full_id = String::from_utf8_lossy(rec.id()).into_owned();
        let name = full_id
            .split_whitespace()
            .next()
            .unwrap_or(&full_id)
            .to_string();
        let hdr_lc = full_id.to_ascii_lowercase();
        let circular = match mode {
            CircularMode::All => true,
            CircularMode::None => false,
            CircularMode::Auto => {
                !(hdr_lc.contains("circular=false") || hdr_lc.contains("topology=linear") || hdr_lc == "linear")
            }
        };
        let seq: Vec<u8> = rec.seq().iter().map(|b| b.to_ascii_uppercase()).collect();
        if seq.is_empty() {
            continue;
        }
        let ext_seq = build_circular_ext_seq(&seq, circular);
        let (true_cpg_dyads, true_chg_sites, true_chh_sites, true_a_sites) =
            count_contig_motifs(&seq, circular);
        out.push(TargetContig {
            name,
            seq,
            ext_seq,
            circular,
            true_cpg_dyads,
            true_chg_sites,
            true_chh_sites,
            true_a_sites,
        });
    }
    ensure!(!out.is_empty(), "no sequences in {}", path.display());
    Ok(out)
}

fn load_targets_from_gfa(path: &Path) -> Result<Vec<TargetContig>> {
    let f = BufReader::new(File::open(path).with_context(|| format!("cannot read GFA {}", path.display()))?);
    let mut segs: Vec<(String, Vec<u8>)> = Vec::new();
    let mut self_loops: rustc_hash::FxHashSet<String> = rustc_hash::FxHashSet::default();
    for line in f.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix("S\t") {
            let mut cols = rest.split('\t');
            let Some(name) = cols.next() else { continue };
            let Some(seq_str) = cols.next() else { continue };
            if seq_str == "*" || seq_str.is_empty() {
                continue;
            }
            segs.push((name.to_string(), seq_str.bytes().map(|b| b.to_ascii_uppercase()).collect()));
        } else if let Some(rest) = line.strip_prefix("L\t") {
            let cols: Vec<&str> = rest.split('\t').collect();
            if cols.len() >= 4 && cols[0] == cols[2] && cols[1] == cols[3] {
                self_loops.insert(cols[0].to_string());
            }
        }
    }
    ensure!(!segs.is_empty(), "no sequence segments in GFA {}", path.display());
    let mut out = Vec::with_capacity(segs.len());
    for (name, seq) in segs {
        let circular = self_loops.contains(&name);
        let ext_seq = build_circular_ext_seq(&seq, circular);
        let (true_cpg_dyads, true_chg_sites, true_chh_sites, true_a_sites) =
            count_contig_motifs(&seq, circular);
        out.push(TargetContig {
            name,
            seq,
            ext_seq,
            circular,
            true_cpg_dyads,
            true_chg_sites,
            true_chh_sites,
            true_a_sites,
        });
    }
    Ok(out)
}

impl TargetIndex {
    fn build(contigs: Vec<TargetContig>, k: usize) -> Self {
        let bits_pow2 = 1usize << 25;
        let words = bits_pow2 >> 6;
        let bloom_mask = bits_pow2 - 1;
        let mut bloom = vec![0u64; words];
        let mut map: FxHashMap<u64, Vec<KmerHit>> = FxHashMap::default();

        for (cid, c) in contigs.iter().enumerate() {
            scan_canonical_kmers(&c.ext_seq, k, 1, |pos, canon, is_rc| {
                let h = hash64(canon);
                let b1 = (h as usize) & bloom_mask;
                let b2 = ((h >> 25) as usize) & bloom_mask;
                bloom[b1 >> 6] |= 1u64 << (b1 & 63);
                bloom[b2 >> 6] |= 1u64 << (b2 & 63);
                let vec = map.entry(canon).or_default();
                if vec.len() <= MAX_KMER_OCC {
                    vec.push(KmerHit {
                        contig: cid as u16,
                        pos: pos as u32,
                        is_rc,
                    });
                }
            });
        }
        map.retain(|_, v| v.len() <= MAX_KMER_OCC);
        Self {
            k,
            contigs,
            bloom,
            bloom_mask,
            map,
        }
    }

    #[inline]
    fn may_contain(&self, canon: u64) -> bool {
        let h = hash64(canon);
        let b1 = (h as usize) & self.bloom_mask;
        let b2 = ((h >> 25) as usize) & self.bloom_mask;
        ((self.bloom[b1 >> 6] >> (b1 & 63)) & 1) != 0
            && ((self.bloom[b2 >> 6] >> (b2 & 63)) & 1) != 0
    }

    fn quick_screen(&self, seq: &[u8]) -> bool {
        let mut hits = 0usize;
        let mut done = false;
        scan_canonical_kmers(seq, self.k, 16, |_, canon, _| {
            if !done && self.may_contain(canon) && self.map.contains_key(&canon) {
                hits += 1;
                if hits >= 2 {
                    done = true;
                }
            }
        });
        done
    }
}

#[derive(Clone, Debug)]
struct AnchorChain {
    contig: u16,
    is_rev: bool,
    anchors: Vec<(u32, u32)>,
    read_start: u32,
    read_end: u32,
    covered_bases: u32,
}

fn chain_read(idx: &TargetIndex, orig_seq: &[u8]) -> Vec<AnchorChain> {
    let k = idx.k as u32;
    let mut bins: FxHashMap<(u16, bool, i32), Vec<(u32, u32)>> = FxHashMap::default();
    scan_canonical_kmers(orig_seq, idx.k, 4, |r_pos, canon, r_is_rc| {
        if !idx.may_contain(canon) {
            return;
        }
        let Some(hits) = idx.map.get(&canon) else {
            return;
        };
        let rp = r_pos as u32;
        for h in hits {
            let is_rev = r_is_rc != h.is_rc;
            let diag: i64 = if !is_rev {
                h.pos as i64 - rp as i64
            } else {
                h.pos as i64 + rp as i64
            };
            let bin = (diag >> 8) as i32;
            bins.entry((h.contig, is_rev, bin)).or_default().push((rp, h.pos));
        }
    });

    if bins.is_empty() {
        return Vec::new();
    }

    let mut groups: FxHashMap<(u16, bool), Vec<(i32, Vec<(u32, u32)>)>> = FxHashMap::default();
    for ((cid, is_rev, bin), pts) in bins {
        if pts.len() >= 3 {
            groups.entry((cid, is_rev)).or_default().push((bin, pts));
        }
    }

    let mut chains = Vec::new();
    for ((cid, is_rev), mut bvec) in groups {
        bvec.sort_by_key(|(b, _)| *b);
        let mut merged: Vec<Vec<(u32, u32)>> = Vec::new();
        for (_, mut pts) in bvec {
            if let Some(prev_pts) = merged.last_mut() {
                let prev_diag = if !is_rev {
                    prev_pts[0].1 as i64 - prev_pts[0].0 as i64
                } else {
                    prev_pts[0].1 as i64 + prev_pts[0].0 as i64
                };
                let cur_diag = if !is_rev {
                    pts[0].1 as i64 - pts[0].0 as i64
                } else {
                    pts[0].1 as i64 + pts[0].0 as i64
                };
                if (cur_diag - prev_diag).abs() <= 384 {
                    prev_pts.append(&mut pts);
                    continue;
                }
            }
            merged.push(pts);
        }

        for mut pts in merged {
            pts.sort_unstable_by_key(|&(r, t)| (r, t));
            pts.dedup();
            if pts.len() < 3 {
                continue;
            }
            let mut used = vec![false; pts.len()];
            loop {
                let Some(start_i) = used.iter().position(|&u| !u) else {
                    break;
                };
                let mut chain = vec![pts[start_i]];
                used[start_i] = true;
                let (mut last_r, mut last_t) = pts[start_i];
                for j in (start_i + 1)..pts.len() {
                    if used[j] {
                        continue;
                    }
                    let (r, t) = pts[j];
                    if r <= last_r {
                        continue;
                    }
                    let dr = (r - last_r) as i64;
                    let dt = if !is_rev {
                        t as i64 - last_t as i64
                    } else {
                        last_t as i64 - t as i64
                    };
                    if dt > 0 && dr <= 2500 && dt <= 2500 && (dr - dt).abs() <= 200 {
                        chain.push((r, t));
                        used[j] = true;
                        last_r = r;
                        last_t = t;
                    }
                }
                if chain.len() >= 3 {
                    let r0 = chain.first().unwrap().0;
                    let r1 = chain.last().unwrap().0 + k;
                    let mut cov = 0u32;
                    let mut cur_end = 0u32;
                    for &(r, _) in &chain {
                        let s = r.max(cur_end);
                        let e = r + k;
                        if e > s {
                            cov += e - s;
                            cur_end = e;
                        }
                    }
                    let span = r1.saturating_sub(r0);
                    if span >= 150 && cov * 5 >= span {
                        chains.push(AnchorChain {
                            contig: cid,
                            is_rev,
                            anchors: chain,
                            read_start: r0,
                            read_end: r1,
                            covered_bases: cov,
                        });
                    }
                }
            }
        }
    }

    chains.sort_by(|a, b| {
        b.covered_bases
            .cmp(&a.covered_bases)
            .then_with(|| (b.read_end - b.read_start).cmp(&(a.read_end - a.read_start)))
    });
    chains
}

fn align_short_gap(r_slice: &[u8], t_slice: &[u8]) -> Vec<(usize, usize)> {
    let n = r_slice.len();
    let m = t_slice.len();
    if n == 0 || m == 0 {
        return Vec::new();
    }
    if n == m {
        return (0..n).map(|i| (i, i)).collect();
    }
    if n > 450 || m > 450 {
        let len = n.min(m);
        return (0..len).map(|i| (i, i)).collect();
    }
    let cols = m + 1;
    let neg_inf = i32::MIN / 4;
    let mut dp = vec![neg_inf; (n + 1) * cols];
    let mut bt = vec![0u8; (n + 1) * cols];
    let band = n.abs_diff(m) + 24;
    dp[0] = 0;
    for j in 1..=m.min(band) {
        dp[j] = -(j as i32) * 2;
        bt[j] = 3;
    }
    for i in 1..=n {
        let j_min = i.saturating_sub(band).max(1);
        let j_max = (i + band).min(m);
        if i <= band {
            dp[i * cols] = -(i as i32) * 2;
            bt[i * cols] = 2;
        }
        for j in j_min..=j_max {
            let score_m = if r_slice[i - 1].eq_ignore_ascii_case(&t_slice[j - 1]) {
                2
            } else {
                -2
            };
            let d = dp[(i - 1) * cols + (j - 1)] + score_m;
            let u = dp[(i - 1) * cols + j] - 2;
            let l = dp[i * cols + (j - 1)] - 2;
            let idx = i * cols + j;
            if d >= u && d >= l {
                dp[idx] = d;
                bt[idx] = 1;
            } else if u >= l {
                dp[idx] = u;
                bt[idx] = 2;
            } else {
                dp[idx] = l;
                bt[idx] = 3;
            }
        }
    }
    let (mut i, mut j) = (n, m);
    let mut pairs = Vec::with_capacity(n.min(m));
    while i > 0 && j > 0 {
        match bt[i * cols + j] {
            1 => {
                pairs.push((i - 1, j - 1));
                i -= 1;
                j -= 1;
            }
            2 => i -= 1,
            3 => j -= 1,
            _ => break,
        }
    }
    pairs.reverse();
    pairs
}

fn build_base_map(
    idx: &TargetIndex,
    orig_seq: &[u8],
    selected: &[AnchorChain],
) -> Vec<Option<(u16, u32, bool)>> {
    let mut map = vec![None; orig_seq.len()];
    let k = idx.k;
    for ch in selected {
        let contig = &idx.contigs[ch.contig as usize];
        let true_len = contig.seq.len();
        let ext = &contig.ext_seq;
        let cid = ch.contig;
        if !ch.is_rev {
            for w in 0..ch.anchors.len() {
                let (r, t) = ch.anchors[w];
                let (ru, tu) = (r as usize, t as usize);
                for off in 0..k {
                    if ru + off < map.len() && tu + off < ext.len() && map[ru + off].is_none() {
                        map[ru + off] = Some((cid, ((tu + off) % true_len) as u32, true));
                    }
                }
                if w + 1 < ch.anchors.len() {
                    let (rn, tn) = ch.anchors[w + 1];
                    let r_gap_s = ru + k;
                    let r_gap_e = rn as usize;
                    let t_gap_s = tu + k;
                    let t_gap_e = tn as usize;
                    if r_gap_e > r_gap_s && t_gap_e > t_gap_s && r_gap_e <= orig_seq.len() && t_gap_e <= ext.len() {
                        for (ro, to) in align_short_gap(&orig_seq[r_gap_s..r_gap_e], &ext[t_gap_s..t_gap_e]) {
                            if map[r_gap_s + ro].is_none() {
                                map[r_gap_s + ro] = Some((cid, ((t_gap_s + to) % true_len) as u32, true));
                            }
                        }
                    }
                }
            }
        } else {
            for w in 0..ch.anchors.len() {
                let (r, t) = ch.anchors[w];
                let (ru, tu) = (r as usize, t as usize);
                for off in 0..k {
                    if ru + off < map.len() && tu + k - 1 >= off && tu + k - 1 - off < ext.len() && map[ru + off].is_none() {
                        let tp = (tu + k - 1 - off) % true_len;
                        map[ru + off] = Some((cid, tp as u32, false));
                    }
                }
                if w + 1 < ch.anchors.len() {
                    let (rn, tn) = ch.anchors[w + 1];
                    let r_gap_s = ru + k;
                    let r_gap_e = rn as usize;
                    let t_high = tu;
                    let t_low = (tn as usize) + k;
                    if r_gap_e > r_gap_s && t_high > t_low && r_gap_e <= orig_seq.len() && t_high <= ext.len() {
                        let mut t_rc = ext[t_low..t_high].to_vec();
                        revcomp_in_place(&mut t_rc);
                        let gap_len = t_high - t_low;
                        for (ro, to) in align_short_gap(&orig_seq[r_gap_s..r_gap_e], &t_rc) {
                            if map[r_gap_s + ro].is_none() {
                                let ext_pos = t_low + (gap_len - 1 - to);
                                map[r_gap_s + ro] = Some((cid, (ext_pos % true_len) as u32, false));
                            }
                        }
                    }
                }
            }
        }
    }
    map
}

// =============================================================================================
// Placed Read Representation, NUMT/NUPT Separation, and Site Pileup
// =============================================================================================

#[derive(Clone, Copy, Debug)]
struct ProjectedSite {
    contig: u16,
    pos: u32,
    plus_strand: bool,
    context: ModContext,
    qual: u8,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ReadClass {
    Organelle,
    NumtHypermethylated,
    NumtChimeric,
}

impl ReadClass {
    fn as_str(self) -> &'static str {
        match self {
            Self::Organelle => "organelle",
            Self::NumtHypermethylated => "numt_hypermethylated",
            Self::NumtChimeric => "numt_chimeric",
        }
    }
}

struct PlacedRead {
    name: String,
    read_len: u32,
    primary_contig: u16,
    organelle_span: u32,
    nuclear_overhang: u32,
    summary: ReadModSummary,
    sites: Vec<ProjectedSite>,
}

fn process_candidate_read(
    idx: &TargetIndex,
    read: &BamRead,
    prob_threshold: f64,
) -> Option<PlacedRead> {
    let chains = chain_read(idx, &read.seq);
    if chains.is_empty() {
        return None;
    }
    let top_cov = chains[0].covered_bases;
    let top_span = chains[0].read_end - chains[0].read_start;
    if top_span < 250 {
        return None;
    }
    let tied: Vec<usize> = chains
        .iter()
        .enumerate()
        .filter(|(_, c)| {
            c.read_start.abs_diff(chains[0].read_start) <= 64
                && c.read_end.abs_diff(chains[0].read_end) <= 64
                && c.covered_bases * 100 >= top_cov * 98
        })
        .map(|(i, _)| i)
        .collect();
    let chosen_first = if tied.len() > 1 {
        tied[(fnv1a(&read.name) as usize) % tied.len()]
    } else {
        0
    };

    let mut selected = vec![chains[chosen_first].clone()];
    for (i, c) in chains.iter().enumerate() {
        if i == chosen_first || tied.contains(&i) {
            continue;
        }
        let overlap = selected
            .iter()
            .map(|s| s.read_end.min(c.read_end).saturating_sub(s.read_start.max(c.read_start)))
            .max()
            .unwrap_or(0);
        let c_span = c.read_end - c.read_start;
        if c_span >= 200 && overlap * 4 < c_span {
            selected.push(c.clone());
            if selected.len() >= 4 {
                break;
            }
        }
    }

    let rlen = read.seq.len() as u32;
    let min_start = selected.iter().map(|c| c.read_start).min().unwrap_or(0);
    let max_end = selected.iter().map(|c| c.read_end).max().unwrap_or(0);
    let mut intervals: Vec<(u32, u32)> = selected.iter().map(|c| (c.read_start, c.read_end)).collect();
    intervals.sort_unstable();
    let mut covered_span = 0u32;
    let mut cur_e = 0u32;
    for (s, e) in intervals {
        let s2 = s.max(cur_e);
        if e > s2 {
            covered_span += e - s2;
            cur_e = e;
        }
    }
    let nuclear_overhang = min_start.max(rlen.saturating_sub(max_end));

    let calls = match (read.mm.as_deref(), read.ml.as_deref()) {
        (Some(mm), Some(ml)) => decode_read_mods(&read.seq, mm, ml),
        _ => Vec::new(),
    };
    let summary = summarize_calls(&calls, prob_threshold);
    let base_map = build_base_map(idx, &read.seq, &selected);

    let mut sites = Vec::with_capacity(calls.len());
    for call in &calls {
        let Some((cid, t_pos, read_plus_on_target)) = base_map.get(call.read_pos as usize).and_then(|&x| x) else {
            continue;
        };
        let contig = &idx.contigs[cid as usize];
        let true_len = contig.seq.len();
        let plus_on_target = read_plus_on_target == call.on_orig_strand;
        match call.context {
            ModContext::Cg | ModContext::Chg | ModContext::Chh => {
                let Some(ref_ctx) =
                    contig_c_context(&contig.seq, contig.circular, t_pos as usize, plus_on_target)
                else {
                    continue;
                };
                if ref_ctx == ModContext::Cg {
                    let dyad_pos = if plus_on_target {
                        t_pos as usize
                    } else if t_pos > 0 {
                        (t_pos - 1) as usize
                    } else if contig.circular {
                        true_len - 1
                    } else {
                        continue;
                    };
                    sites.push(ProjectedSite {
                        contig: cid,
                        pos: dyad_pos as u32,
                        plus_strand: true,
                        context: ModContext::Cg,
                        qual: call.qual,
                    });
                } else {
                    sites.push(ProjectedSite {
                        contig: cid,
                        pos: t_pos,
                        plus_strand: plus_on_target,
                        context: ref_ctx,
                        qual: call.qual,
                    });
                }
            }
            ModContext::M6a => {
                let expect = if plus_on_target { b'A' } else { b'T' };
                if contig.seq[t_pos as usize] == expect {
                    sites.push(ProjectedSite {
                        contig: cid,
                        pos: t_pos,
                        plus_strand: plus_on_target,
                        context: ModContext::M6a,
                        qual: call.qual,
                    });
                }
            }
        }
    }

    Some(PlacedRead {
        name: String::from_utf8_lossy(&read.name).into_owned(),
        read_len: rlen,
        primary_contig: selected[0].contig,
        organelle_span: covered_span,
        nuclear_overhang,
        summary,
        sites,
    })
}

fn quantile_sorted(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let idx = ((sorted.len() - 1) as f64 * q.clamp(0.0, 1.0)).round() as usize;
    sorted[idx.min(sorted.len() - 1)]
}

fn roc_auc(neg_organelle: &[f64], pos_nuclear: &[f64]) -> Option<f64> {
    if neg_organelle.is_empty() || pos_nuclear.is_empty() {
        return None;
    }
    let mut all: Vec<(f64, bool)> = Vec::with_capacity(neg_organelle.len() + pos_nuclear.len());
    for &v in neg_organelle {
        all.push((v, false));
    }
    for &v in pos_nuclear {
        all.push((v, true));
    }
    all.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut rank_sum_pos = 0.0f64;
    let mut i = 0usize;
    while i < all.len() {
        let mut j = i + 1;
        while j < all.len() && all[j].0 == all[i].0 {
            j += 1;
        }
        let avg_rank = (i + 1 + j) as f64 * 0.5;
        let pos_cnt = all[i..j].iter().filter(|x| x.1).count() as f64;
        rank_sum_pos += pos_cnt * avg_rank;
        i = j;
    }
    let n_pos = pos_nuclear.len() as f64;
    let n_neg = neg_organelle.len() as f64;
    let u = rank_sum_pos - n_pos * (n_pos + 1.0) * 0.5;
    Some((u / (n_pos * n_neg)).clamp(0.0, 1.0))
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum EpigeneticRegime {
    PlantCgChgRddm,
    VertebrateGlobalCpg,
    InvertebrateMosaicCpg,
    InvertebrateUnmethylated,
    FungiVariable5mC,
}

impl EpigeneticRegime {
    fn as_str(self) -> &'static str {
        match self {
            Self::PlantCgChgRddm => "plant_cg_chg_rddm",
            Self::VertebrateGlobalCpg => "vertebrate_global_cpg",
            Self::InvertebrateMosaicCpg => "invertebrate_mosaic_cpg",
            Self::InvertebrateUnmethylated => "invertebrate_unmethylated",
            Self::FungiVariable5mC => "fungi_variable_5mc",
        }
    }

    fn note(self) -> &'static str {
        match self {
            Self::PlantCgChgRddm => {
                "Land-plant MET1/CMT3/RdDM regime: nuclear DNA carries heavy CG + CHG (+ CHH) 5mC while mtDNA/ptDNA stay unmethylated."
            }
            Self::VertebrateGlobalCpg => {
                "Vertebrate global DNMT1/DNMT3 CpG regime: nuclear NUMTs are heavily CpG-methylated (~65-85%) while mtDNA stays unmethylated (~0.5-2.5%); CHG rescue disabled."
            }
            Self::InvertebrateMosaicCpg => {
                "Invertebrate mosaic CpG regime: gene-body-integrated hypermethylated NUMTs are purged; unmethylated intergenic NUMTs require overhang/depth/graph evidence."
            }
            Self::InvertebrateUnmethylated => {
                "DNMT-deficient / unmethylated animal regime (e.g. Diptera, Nematoda): nuclear DNA lacks 5mC; conservative guard prevents false pruning of AT-rich mtDNA reads."
            }
            Self::FungiVariable5mC => {
                "Fungal regime: repeat-targeted 5mC varies across clades (present in Neurospora/Basidiomycota, absent in Saccharomycotina); CG-only filter applied."
            }
        }
    }
}

fn infer_epigenetic_regime(
    lineage: crate::seeddb::Lineage,
    nuclear_cgs: &[f64],
    placed_full_cgs: &[f64],
) -> EpigeneticRegime {
    match lineage {
        crate::seeddb::Lineage::Landplants => EpigeneticRegime::PlantCgChgRddm,
        crate::seeddb::Lineage::Fungi => EpigeneticRegime::FungiVariable5mC,
        crate::seeddb::Lineage::Metazoa => {
            let mut nuc = nuclear_cgs.to_vec();
            nuc.sort_by(|a, b| a.total_cmp(b));
            let mut plc = placed_full_cgs.to_vec();
            plc.sort_by(|a, b| a.total_cmp(b));
            let high_placed_frac = if plc.is_empty() {
                0.0
            } else {
                plc.iter().filter(|&&x| x >= 0.45).count() as f64 / plc.len() as f64
            };
            if nuc.len() >= 10 {
                let med = quantile_sorted(&nuc, 0.50);
                let p90 = quantile_sorted(&nuc, 0.90);
                if med >= 0.45 {
                    EpigeneticRegime::VertebrateGlobalCpg
                } else if p90 >= 0.25 || med >= 0.06 || high_placed_frac >= 0.02 {
                    EpigeneticRegime::InvertebrateMosaicCpg
                } else {
                    EpigeneticRegime::InvertebrateUnmethylated
                }
            } else if high_placed_frac >= 0.03 {
                EpigeneticRegime::VertebrateGlobalCpg
            } else if high_placed_frac >= 0.01 || quantile_sorted(&plc, 0.90) >= 0.15 {
                EpigeneticRegime::InvertebrateMosaicCpg
            } else {
                EpigeneticRegime::InvertebrateUnmethylated
            }
        }
    }
}

fn parse_cg_threshold_policy(
    spec: &str,
    organelle_full_cgs: &[f64],
    regime: EpigeneticRegime,
) -> Result<Option<(f64, &'static str)>> {
    let s = spec.trim();
    if s.eq_ignore_ascii_case("none") || s.eq_ignore_ascii_case("off") {
        return Ok(None);
    }
    if s.eq_ignore_ascii_case("auto") {
        let mut v = organelle_full_cgs.to_vec();
        v.sort_by(|a, b| a.total_cmp(b));
        let baseline = if v.is_empty() {
            match regime {
                EpigeneticRegime::VertebrateGlobalCpg
                | EpigeneticRegime::InvertebrateMosaicCpg
                | EpigeneticRegime::InvertebrateUnmethylated => 0.03,
                _ => 0.12,
            }
        } else {
            quantile_sorted(&v, 0.35)
        };
        let (thr, tag) = match regime {
            EpigeneticRegime::PlantCgChgRddm => ((baseline + 0.28).clamp(0.38, 0.55), "auto"),
            EpigeneticRegime::VertebrateGlobalCpg => {
                ((baseline + 0.22).clamp(0.25, 0.45), "auto_vertebrate_cpg")
            }
            EpigeneticRegime::InvertebrateMosaicCpg => {
                ((baseline + 0.26).clamp(0.32, 0.50), "auto_invertebrate_mosaic")
            }
            EpigeneticRegime::InvertebrateUnmethylated => {
                // Unmethylated nuclear genome (e.g. Drosophila, C. elegans): avoid discarding
                // AT-rich mtDNA reads with 5-8 CpGs that have 1-2 noisy basecaller calls.
                (0.55, "auto_unmethylated_clade_guard")
            }
            EpigeneticRegime::FungiVariable5mC => ((baseline + 0.28).clamp(0.35, 0.55), "auto_fungi"),
        };
        return Ok(Some((thr, tag)));
    }
    let val: f64 = s
        .parse()
        .with_context(|| format!("invalid --max-cg-methylation {spec:?}: expected 'auto', 'none', or float in [0, 1]"))?;
    ensure!((0.0..=1.0).contains(&val), "--max-cg-methylation must be between 0.0 and 1.0");
    Ok(Some((val, "explicit")))
}

#[derive(Clone, Copy, Default)]
struct SiteAcc {
    mod_cnt: u32,
    unmod_cnt: u32,
}

/// Run the full `ovasm check methylation` pipeline and return the JSON report.
pub fn run_methylation(args: &MethylationArgs) -> Result<Value> {
    let t0 = Instant::now();
    if let Some(n) = args.threads {
        let _ = rayon::ThreadPoolBuilder::new().num_threads(n).build_global();
    }
    ensure!(
        args.target.is_some() || args.gfa.is_some(),
        "give --target <molecules.fasta> and/or --gfa <organelle.gfa>"
    );
    ensure!((15..=31).contains(&args.kmer), "--kmer must be in 15..=31");
    for b in &args.bam {
        ensure!(b.is_file(), "BAM file not found: {}", b.display());
    }
    fs::create_dir_all(&args.out)
        .with_context(|| format!("cannot create {}", args.out.display()))?;

    let contigs = if let Some(target_fa) = &args.target {
        load_targets_from_fasta(target_fa, args.circular)?
    } else {
        load_targets_from_gfa(args.gfa.as_ref().unwrap())?
    };
    let idx = TargetIndex::build(contigs, args.kmer);

    let mut total_bam_reads = 0u64;
    let mut reads_with_mods = 0u64;
    let mut placed_all: Vec<PlacedRead> = Vec::new();
    let mut nuclear_baseline_cgs: Vec<f64> = Vec::with_capacity(NUCLEAR_SAMPLE_MAX);

    let batch_size = 512usize;
    let mut cand_batch: Vec<BamRead> = Vec::with_capacity(batch_size);

    let flush_batch = |batch: &mut Vec<BamRead>, out: &mut Vec<PlacedRead>| {
        if batch.is_empty() {
            return;
        }
        let placed: Vec<PlacedRead> = batch
            .par_iter()
            .filter_map(|r| process_candidate_read(&idx, r, args.prob_threshold))
            .collect();
        out.extend(placed);
        batch.clear();
    };

    for bam_path in &args.bam {
        for_each_bam_read(bam_path, true, |r| {
            total_bam_reads += 1;
            if r.mm.is_some() && r.ml.is_some() {
                reads_with_mods += 1;
            }
            if idx.quick_screen(&r.seq) {
                cand_batch.push(r);
                if cand_batch.len() >= batch_size {
                    flush_batch(&mut cand_batch, &mut placed_all);
                }
            } else if nuclear_baseline_cgs.len() < NUCLEAR_SAMPLE_MAX && (total_bam_reads & 31) == 0 {
                if let (Some(mm), Some(ml)) = (r.mm.as_deref(), r.ml.as_deref()) {
                    let s = summarize_read_mods(&r.seq, mm, ml, args.prob_threshold);
                    if s.cg_calls >= args.min_read_cpgs.max(10) {
                        nuclear_baseline_cgs.push(s.cg_frac());
                    }
                }
            }
            Ok(true)
        })?;
    }
    flush_batch(&mut cand_batch, &mut placed_all);

    let full_len_cgs: Vec<f64> = placed_all
        .iter()
        .filter(|p| p.nuclear_overhang < 100 && p.summary.cg_calls >= args.min_read_cpgs)
        .map(|p| p.summary.cg_frac())
        .collect();
    let regime = infer_epigenetic_regime(args.lineage, &nuclear_baseline_cgs, &full_len_cgs);
    let thr_policy = parse_cg_threshold_policy(&args.max_cg_methylation, &full_len_cgs, regime)?;
    let allow_plant_chg_rescue = args.lineage == crate::seeddb::Lineage::Landplants;

    let prob_u8 = (args.prob_threshold.clamp(0.0, 1.0) * 255.0).round() as u8;

    let mut unfilt_sites: BTreeMap<(u16, u32, bool, ModContext), SiteAcc> = BTreeMap::new();
    let mut clean_sites: BTreeMap<(u16, u32, bool, ModContext), SiteAcc> = BTreeMap::new();

    let tsv_path = args.out.join("reads.methylation.tsv");
    let mut tsv = BufWriter::with_capacity(256 * 1024, File::create(&tsv_path)?);
    writeln!(
        tsv,
        "read_name\tprimary_target\tread_len\torganelle_span\tnuclear_overhang\tclass\tcg_calls\tcg_frac\tchg_calls\tchg_frac\tchh_calls\tchh_frac\tm6a_calls\tm6a_frac"
    )?;

    let mut n_organelle_reads = 0u64;
    let mut n_numt_hyper = 0u64;
    let mut n_numt_chimeric = 0u64;
    let mut clean_read_cgs = Vec::new();
    let mut numt_read_cgs = Vec::new();

    for p in &placed_all {
        let is_chimeric = p.nuclear_overhang as usize >= args.max_overhang;
        let is_hyper = match thr_policy {
            Some((thr, _)) => {
                (p.summary.cg_calls >= args.min_read_cpgs && p.summary.cg_frac() > thr)
                    || (allow_plant_chg_rescue
                        && p.summary.cg_calls == 0
                        && p.summary.chg_calls >= 10
                        && p.summary.chg_frac() > 0.35)
            }
            None => false,
        };
        let class = if is_hyper {
            n_numt_hyper += 1;
            if p.summary.cg_calls >= args.min_read_cpgs {
                numt_read_cgs.push(p.summary.cg_frac());
            }
            ReadClass::NumtHypermethylated
        } else if is_chimeric {
            n_numt_chimeric += 1;
            if p.summary.cg_calls >= args.min_read_cpgs {
                numt_read_cgs.push(p.summary.cg_frac());
            }
            ReadClass::NumtChimeric
        } else {
            n_organelle_reads += 1;
            if p.summary.cg_calls >= args.min_read_cpgs {
                clean_read_cgs.push(p.summary.cg_frac());
            }
            ReadClass::Organelle
        };

        writeln!(
            tsv,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}\t{:.4}\t{}\t{:.4}\t{}\t{:.4}",
            p.name,
            idx.contigs[p.primary_contig as usize].name,
            p.read_len,
            p.organelle_span,
            p.nuclear_overhang,
            class.as_str(),
            p.summary.cg_calls,
            p.summary.cg_frac(),
            p.summary.chg_calls,
            p.summary.chg_frac(),
            p.summary.chh_calls,
            p.summary.chh_frac(),
            p.summary.m6a_calls,
            p.summary.m6a_frac(),
        )?;

        for s in &p.sites {
            let key = (s.contig, s.pos, s.plus_strand, s.context);
            let is_mod = s.qual >= prob_u8 && s.qual > 0;
            let u_acc = unfilt_sites.entry(key).or_default();
            if is_mod {
                u_acc.mod_cnt += 1;
            } else {
                u_acc.unmod_cnt += 1;
            }
            if class == ReadClass::Organelle {
                let c_acc = clean_sites.entry(key).or_default();
                if is_mod {
                    c_acc.mod_cnt += 1;
                } else {
                    c_acc.unmod_cnt += 1;
                }
            }
        }
    }
    tsv.flush()?;

    let bed_path = args.out.join("methylation.bed");
    let unfilt_bed_path = args.out.join("methylation.unfiltered.bed");
    let write_bed = |path: &Path, table: &BTreeMap<(u16, u32, bool, ModContext), SiteAcc>| -> Result<()> {
        let mut w = BufWriter::with_capacity(256 * 1024, File::create(path)?);
        writeln!(
            w,
            "#chrom\tstart\tend\tcontext\tstrand\tmod_rate\tcov\tmod_count\tunmod_count"
        )?;
        for (&(cid, pos, plus, ctx), acc) in table {
            let cov = acc.mod_cnt + acc.unmod_cnt;
            if cov == 0 {
                continue;
            }
            let rate = acc.mod_cnt as f64 / cov as f64;
            let strand = match ctx {
                ModContext::Cg => '.',
                _ if plus => '+',
                _ => '-',
            };
            let end = if ctx == ModContext::Cg { pos + 2 } else { pos + 1 };
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}",
                idx.contigs[cid as usize].name,
                pos,
                end,
                ctx.as_str(),
                strand,
                rate,
                cov,
                acc.mod_cnt,
                acc.unmod_cnt
            )?;
        }
        w.flush()?;
        Ok(())
    };
    write_bed(&bed_path, &clean_sites)?;
    write_bed(&unfilt_bed_path, &unfilt_sites)?;

    let mut contig_reports = Vec::with_capacity(idx.contigs.len());
    for (cid, c) in idx.contigs.iter().enumerate() {
        let cid_u16 = cid as u16;
        let summarize_ctx = |table: &BTreeMap<(u16, u32, bool, ModContext), SiteAcc>, ctx: ModContext, true_sites: usize| -> Value {
            let mut rates = Vec::new();
            let mut covs = Vec::new();
            let mut sites_ge_20 = 0usize;
            let mut sites_ge_50 = 0usize;
            for (&(tcid, _, _, tctx), acc) in table {
                if tcid == cid_u16 && tctx == ctx {
                    let cov = acc.mod_cnt + acc.unmod_cnt;
                    if cov > 0 {
                        let r = acc.mod_cnt as f64 / cov as f64;
                        rates.push(r);
                        covs.push(cov as f64);
                        if r >= 0.20 {
                            sites_ge_20 += 1;
                        }
                        if r >= 0.50 {
                            sites_ge_50 += 1;
                        }
                    }
                }
            }
            rates.sort_by(|a, b| a.total_cmp(b));
            covs.sort_by(|a, b| a.total_cmp(b));
            json!({
                "called_sites": rates.len(),
                "total_genomic_motifs": true_sites,
                "motif_coverage_fraction": if true_sites > 0 { (rates.len() as f64 / true_sites as f64).min(1.0) } else { 0.0 },
                "median_depth": quantile_sorted(&covs, 0.5),
                "median_mod_rate": quantile_sorted(&rates, 0.5),
                "p95_mod_rate": quantile_sorted(&rates, 0.95),
                "max_mod_rate": rates.last().copied().unwrap_or(0.0),
                "sites_ge_0_20": sites_ge_20,
                "sites_ge_0_50": sites_ge_50,
            })
        };

        contig_reports.push(json!({
            "name": c.name,
            "length": c.seq.len(),
            "circular": c.circular,
            "clean": {
                "CG": summarize_ctx(&clean_sites, ModContext::Cg, c.true_cpg_dyads),
                "CHG": summarize_ctx(&clean_sites, ModContext::Chg, c.true_chg_sites),
                "CHH": summarize_ctx(&clean_sites, ModContext::Chh, c.true_chh_sites),
                "6mA": summarize_ctx(&clean_sites, ModContext::M6a, c.true_a_sites),
            },
            "unfiltered_cg": summarize_ctx(&unfilt_sites, ModContext::Cg, c.true_cpg_dyads),
        }));
    }

    let annotated_gfa = if let Some(gfa_path) = &args.gfa {
        let out_gfa = args.out.join("organelle.methyl.gfa");
        annotate_gfa_segments(gfa_path, &out_gfa, &idx, &clean_sites)?;
        Some(out_gfa)
    } else {
        None
    };

    clean_read_cgs.sort_by(|a, b| a.total_cmp(b));
    numt_read_cgs.sort_by(|a, b| a.total_cmp(b));
    nuclear_baseline_cgs.sort_by(|a, b| a.total_cmp(b));

    let auc_vs_nuclear = roc_auc(&clean_read_cgs, &nuclear_baseline_cgs);
    let auc_vs_numts = roc_auc(&clean_read_cgs, &numt_read_cgs);

    let primary_contexts = match args.lineage {
        crate::seeddb::Lineage::Landplants => json!(["CG", "CHG", "CHH", "6mA"]),
        crate::seeddb::Lineage::Metazoa | crate::seeddb::Lineage::Fungi => json!(["CG", "6mA"]),
    };

    let report = json!({
        "schema": SCHEMA,
        "lineage": args.lineage.name(),
        "epigenetic_regime": regime.as_str(),
        "regime_note": regime.note(),
        "primary_contexts": primary_contexts,
        "parameters": {
            "lineage": args.lineage.name(),
            "prob_threshold": args.prob_threshold,
            "cg_threshold_policy": thr_policy.map_or("none", |(_, m)| m),
            "effective_max_cg_methylation": thr_policy.map(|(v, _)| v),
            "chg_rescue_enabled": allow_plant_chg_rescue,
            "min_read_cpgs": args.min_read_cpgs,
            "max_overhang_bp": args.max_overhang,
            "kmer": args.kmer,
        },
        "reads": {
            "total_bam_primary_reads": total_bam_reads,
            "reads_with_mm_ml_tags": reads_with_mods,
            "organelle_placed_reads": placed_all.len(),
            "clean_organelle_reads": n_organelle_reads,
            "numt_nupt_hypermethylated_reads": n_numt_hyper,
            "numt_nupt_chimeric_reads": n_numt_chimeric,
            "organelle_cg_median": quantile_sorted(&clean_read_cgs, 0.5),
            "organelle_cg_p95": quantile_sorted(&clean_read_cgs, 0.95),
            "numt_nupt_cg_median": if numt_read_cgs.is_empty() { Value::Null } else { json!(quantile_sorted(&numt_read_cgs, 0.5)) },
            "nuclear_baseline_sampled_reads": nuclear_baseline_cgs.len(),
            "nuclear_baseline_cg_median": if nuclear_baseline_cgs.is_empty() { Value::Null } else { json!(quantile_sorted(&nuclear_baseline_cgs, 0.5)) },
            "roc_auc_organelle_vs_nuclear": auc_vs_nuclear,
            "roc_auc_organelle_vs_numt": auc_vs_numts,
        },
        "contigs": contig_reports,
        "files": {
            "bed": bed_path.to_string_lossy(),
            "unfiltered_bed": unfilt_bed_path.to_string_lossy(),
            "per_read_tsv": tsv_path.to_string_lossy(),
            "annotated_gfa": annotated_gfa.as_ref().map(|p| p.to_string_lossy()),
        },
        "elapsed_seconds": t0.elapsed().as_secs_f64(),
    });

    let json_path = args.out.join("methylation.json");
    fs::write(&json_path, serde_json::to_string_pretty(&report)?)?;
    eprintln!(
        "[ovasm] methylation: {} clean organelle reads, {} NUMT/NUPT reads filtered ({} hypermethylated, {} chimeric); organelle CG median {:.3}; {:.1}s",
        n_organelle_reads,
        n_numt_hyper + n_numt_chimeric,
        n_numt_hyper,
        n_numt_chimeric,
        quantile_sorted(&clean_read_cgs, 0.5),
        t0.elapsed().as_secs_f64()
    );
    Ok(report)
}

fn annotate_gfa_segments(
    gfa_in: &Path,
    gfa_out: &Path,
    idx: &TargetIndex,
    clean_sites: &BTreeMap<(u16, u32, bool, ModContext), SiteAcc>,
) -> Result<()> {
    let text = fs::read_to_string(gfa_in)
        .with_context(|| format!("cannot read GFA {}", gfa_in.display()))?;
    let mut out = BufWriter::new(File::create(gfa_out)?);
    let name_to_cid: FxHashMap<&str, u16> = idx
        .contigs
        .iter()
        .enumerate()
        .map(|(i, c)| (c.name.as_str(), i as u16))
        .collect();

    for line in text.lines() {
        if !line.starts_with("S\t") {
            writeln!(out, "{line}")?;
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 || cols[2] == "*" {
            writeln!(out, "{line}")?;
            continue;
        }
        let seg_name = cols[1];
        let seg_seq = cols[2].as_bytes();
        let mut cg_rates = Vec::new();
        let mut cg_covs = Vec::new();
        let mut m6a_rates = Vec::new();

        if let Some(&cid) = name_to_cid.get(seg_name) {
            for (&(tcid, _, _, ctx), acc) in clean_sites {
                if tcid == cid {
                    let cov = acc.mod_cnt + acc.unmod_cnt;
                    if cov > 0 {
                        let r = acc.mod_cnt as f64 / cov as f64;
                        if ctx == ModContext::Cg {
                            cg_rates.push(r);
                            cg_covs.push(cov as f64);
                        } else if ctx == ModContext::M6a {
                            m6a_rates.push(r);
                        }
                    }
                }
            }
        } else if let Some(ch) = chain_read(idx, seg_seq).into_iter().next() {
            let base_map = build_base_map(idx, seg_seq, &[ch]);
            for mapped in base_map.into_iter().flatten() {
                let (cid, tpos, _) = mapped;
                if let Some(acc) = clean_sites.get(&(cid, tpos, true, ModContext::Cg)) {
                    let cov = acc.mod_cnt + acc.unmod_cnt;
                    if cov > 0 {
                        cg_rates.push(acc.mod_cnt as f64 / cov as f64);
                        cg_covs.push(cov as f64);
                    }
                }
                for s in [true, false] {
                    if let Some(acc) = clean_sites.get(&(cid, tpos, s, ModContext::M6a)) {
                        let cov = acc.mod_cnt + acc.unmod_cnt;
                        if cov > 0 {
                            m6a_rates.push(acc.mod_cnt as f64 / cov as f64);
                        }
                    }
                }
            }
        }

        cg_rates.sort_by(|a, b| a.total_cmp(b));
        cg_covs.sort_by(|a, b| a.total_cmp(b));
        m6a_rates.sort_by(|a, b| a.total_cmp(b));

        write!(
            out,
            "{line}\tMC:f:{:.4}\tMP:i:{}\tMD:f:{:.1}",
            quantile_sorted(&cg_rates, 0.5),
            cg_rates.len(),
            quantile_sorted(&cg_covs, 0.5)
        )?;
        if !m6a_rates.is_empty() {
            write!(out, "\tMA:f:{:.4}", quantile_sorted(&m6a_rates, 0.5))?;
        }
        writeln!(out)?;
    }
    out.flush()?;
    Ok(())
}

// =============================================================================================
// Unit Tests (Synthetic BAM + SAMtags MM/ML + Circular & 2-IR Plastid / NUMT Separation)
// =============================================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;

    fn pack_bam_seq(seq: &[u8]) -> Vec<u8> {
        let code = |b: u8| -> u8 {
            match b.to_ascii_uppercase() {
                b'A' => 1,
                b'C' => 2,
                b'G' => 4,
                b'T' => 8,
                _ => 15,
            }
        };
        let mut out = vec![0u8; (seq.len() + 1) / 2];
        for (i, &b) in seq.iter().enumerate() {
            if i & 1 == 0 {
                out[i >> 1] |= code(b) << 4;
            } else {
                out[i >> 1] |= code(b);
            }
        }
        out
    }

    fn write_synthetic_bam(path: &Path, records: &[(&str, u16, &[u8], Option<&str>, Option<&[u8]>)]) {
        let mut raw = Vec::new();
        raw.extend_from_slice(b"BAM\x01");
        let hdr = b"@HD\tVN:1.6\n";
        raw.extend_from_slice(&(hdr.len() as i32).to_le_bytes());
        raw.extend_from_slice(hdr);
        raw.extend_from_slice(&0i32.to_le_bytes());

        for &(name, flag, orig_seq, mm_opt, ml_opt) in records {
            let mut stored_seq = orig_seq.to_vec();
            if flag & 0x10 != 0 {
                revcomp_in_place(&mut stored_seq);
            }
            let mut rec = Vec::new();
            rec.extend_from_slice(&(-1i32).to_le_bytes());
            rec.extend_from_slice(&(-1i32).to_le_bytes());
            let l_name = (name.len() + 1) as u8;
            rec.push(l_name);
            rec.push(60u8);
            rec.extend_from_slice(&0u16.to_le_bytes());
            rec.extend_from_slice(&0u16.to_le_bytes());
            rec.extend_from_slice(&flag.to_le_bytes());
            rec.extend_from_slice(&(stored_seq.len() as i32).to_le_bytes());
            rec.extend_from_slice(&(-1i32).to_le_bytes());
            rec.extend_from_slice(&(-1i32).to_le_bytes());
            rec.extend_from_slice(&0i32.to_le_bytes());
            rec.extend_from_slice(name.as_bytes());
            rec.push(0);
            rec.extend_from_slice(&pack_bam_seq(&stored_seq));
            rec.extend(std::iter::repeat(30u8).take(stored_seq.len()));
            if let Some(mm) = mm_opt {
                rec.extend_from_slice(b"MMZ");
                rec.extend_from_slice(mm.as_bytes());
                rec.push(0);
            }
            if let Some(ml) = ml_opt {
                rec.extend_from_slice(b"MLBC");
                rec.extend_from_slice(&(ml.len() as u32).to_le_bytes());
                rec.extend_from_slice(ml);
            }
            raw.extend_from_slice(&(rec.len() as i32).to_le_bytes());
            raw.extend_from_slice(&rec);
        }
        let mut gz = GzEncoder::new(File::create(path).unwrap(), Compression::fast());
        gz.write_all(&raw).unwrap();
        gz.finish().unwrap();
    }

    #[test]
    fn decodes_plant_cg_chg_chh_and_6ma_across_forward_and_reverse_bam_records() {
        let seq = b"ACGTCAGTCATAACG";
        let mm = "C+m.,0,0,0,0;A+a?,1,2;";
        let ml: &[u8] = &[10, 230, 200, 5, 240, 12];

        let dir = std::env::temp_dir().join(format!("ovasm_meth_unit_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("test.bam");
        write_synthetic_bam(
            &bam,
            &[
                ("r_fwd", 4, seq, Some(mm), Some(ml)),
                ("r_rev", 16, seq, Some(mm), Some(ml)),
                ("r_sec", 0x100, seq, Some(mm), Some(ml)),
            ],
        );

        let mut seen = Vec::new();
        for_each_bam_read(&bam, true, |r| {
            let calls = decode_read_mods(&r.seq, r.mm.as_deref().unwrap(), r.ml.as_deref().unwrap());
            seen.push((String::from_utf8(r.name).unwrap(), r.seq, calls));
            Ok(true)
        })
        .unwrap();

        assert_eq!(seen.len(), 2, "secondary alignment must be skipped");
        for (name, rseq, calls) in &seen {
            assert_eq!(rseq.as_slice(), seq, "failed orientation restoration for {name}");
            let s = summarize_calls(calls, 0.5);
            assert_eq!(s.cg_calls, 2);
            assert_eq!(s.cg_mod, 0);
            assert_eq!(s.chg_calls, 1);
            assert_eq!(s.chg_mod, 1);
            assert_eq!(s.chh_calls, 1);
            assert_eq!(s.chh_mod, 1);
            assert_eq!(s.m6a_calls, 2);
            assert_eq!(s.m6a_mod, 1);
        }
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn preserves_plastid_inverted_repeats_circular_origin_and_purges_numts() {
        let mut rng: u64 = 0x1234_5678_9abc_def0;
        let mut rand_seq = |len: usize| -> Vec<u8> {
            (0..len)
                .map(|_| {
                    rng = hash64(rng.wrapping_add(1));
                    b"ACGT"[(rng & 3) as usize]
                })
                .collect()
        };
        let lsc = rand_seq(14_000);
        let ira = rand_seq(8_000);
        let ssc = rand_seq(6_000);
        let mut irb = ira.clone();
        revcomp_in_place(&mut irb);

        let mut plastid = Vec::with_capacity(36_000);
        plastid.extend_from_slice(&lsc);
        plastid.extend_from_slice(&ira);
        plastid.extend_from_slice(&ssc);
        plastid.extend_from_slice(&irb);

        let nuclear_flank = rand_seq(2_000);

        let dir = std::env::temp_dir().join(format!("ovasm_meth_ir_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let target_fa = dir.join("plastid.fasta");
        fs::write(
            &target_fa,
            format!(">plastid circular=true\n{}\n", String::from_utf8(plastid.clone()).unwrap()),
        )
        .unwrap();

        let make_cpg_mm_ml = |seq: &[u8], ml_val: u8| -> (String, Vec<u8>) {
            let mut deltas = Vec::new();
            let mut skip = 0usize;
            let mut mls = Vec::new();
            for i in 0..seq.len() {
                if seq[i] == b'C' {
                    if i + 1 < seq.len() && seq[i + 1] == b'G' {
                        deltas.push(skip.to_string());
                        mls.push(ml_val);
                        skip = 0;
                    } else {
                        skip += 1;
                    }
                }
            }
            let mm = if deltas.is_empty() {
                String::new()
            } else {
                format!("C+m?,{};", deltas.join(","))
            };
            (mm, mls)
        };

        let mut raw_records: Vec<(String, u16, Vec<u8>, String, Vec<u8>)> = Vec::new();
        for i in 0..72 {
            let start = (i * 500) % plastid.len();
            let len = 5_000usize;
            let mut rseq: Vec<u8> = (0..len).map(|k| plastid[(start + k) % plastid.len()]).collect();
            let flag = if i % 2 == 0 { 4u16 } else { 16u16 };
            if i % 3 == 0 {
                revcomp_in_place(&mut rseq);
            }
            let (mm, ml) = make_cpg_mm_ml(&rseq, 15);
            raw_records.push((format!("org_{i}"), flag, rseq, mm, ml));
        }
        for i in 0..15 {
            let rseq = plastid[2_000..7_000].to_vec();
            let (mm, ml) = make_cpg_mm_ml(&rseq, 235);
            raw_records.push((format!("nupt_internal_{i}"), 4, rseq, mm, ml));
        }
        for i in 0..5 {
            let mut rseq = nuclear_flank.clone();
            rseq.extend_from_slice(&plastid[0..3_000]);
            let (mm, ml) = make_cpg_mm_ml(&rseq, 220);
            raw_records.push((format!("nupt_chimera_{i}"), 4, rseq, mm, ml));
        }

        let bam_path = dir.join("reads.bam");
        let refs: Vec<(&str, u16, &[u8], Option<&str>, Option<&[u8]>)> = raw_records
            .iter()
            .map(|(n, fl, sq, mm, ml)| (n.as_str(), *fl, sq.as_slice(), Some(mm.as_str()), Some(ml.as_slice())))
            .collect();
        write_synthetic_bam(&bam_path, &refs);

        let out_dir = dir.join("out");
        let report = run_methylation(&MethylationArgs {
            target: Some(target_fa),
            gfa: None,
            bam: vec![bam_path],
            out: out_dir,
            prob_threshold: 0.5,
            max_cg_methylation: "auto".into(),
            min_read_cpgs: 5,
            max_overhang: 300,
            kmer: 31,
            circular: CircularMode::Auto,
            lineage: crate::seeddb::Lineage::Landplants,
            threads: Some(2),
        })
        .unwrap();

        assert_eq!(report["reads"]["clean_organelle_reads"], json!(72));
        assert_eq!(report["reads"]["numt_nupt_hypermethylated_reads"], json!(20));
        let cg = &report["contigs"][0]["clean"]["CG"];
        let cov_frac = cg["motif_coverage_fraction"].as_f64().unwrap();
        assert!(
            cov_frac >= 0.995,
            "expected >=99.5% CpG coverage across LSC + 8kb IRa/IRb + SSC + circular origin, got {cov_frac}"
        );
        assert_eq!(cg["sites_ge_0_20"], json!(0), "all NUPT false methylation spikes must be purged");
        assert!(
            report["contigs"][0]["unfiltered_cg"]["sites_ge_0_20"].as_u64().unwrap() > 50,
            "unfiltered pileup should show the NUPT spike before filtering"
        );
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn animal_mitogenome_adapts_to_vertebrate_rolling_circle_mega_numts_and_unmethylated_insects() {
        let mut rng: u64 = 0xfeed_face_cafe_beef;
        let mut rand_seq = |len: usize| -> Vec<u8> {
            (0..len)
                .map(|_| {
                    rng = hash64(rng.wrapping_add(1));
                    b"ACGT"[(rng & 3) as usize]
                })
                .collect()
        };
        // 16.5 kb circular animal mitogenome (~human/vertebrate scale)
        let mito = rand_seq(16_569);
        let dir = std::env::temp_dir().join(format!("ovasm_meth_animal_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let target_fa = dir.join("animal_mt.fasta");
        fs::write(
            &target_fa,
            format!(">mt circular=true\n{}\n", String::from_utf8(mito.clone()).unwrap()),
        )
        .unwrap();

        let make_cpg_mm_ml = |seq: &[u8], ml_val: u8| -> (String, Vec<u8>) {
            let mut deltas = Vec::new();
            let mut skip = 0usize;
            let mut mls = Vec::new();
            for i in 0..seq.len() {
                if seq[i] == b'C' {
                    if i + 1 < seq.len() && seq[i + 1] == b'G' {
                        deltas.push(skip.to_string());
                        mls.push(ml_val);
                        skip = 0;
                    } else {
                        skip += 1;
                    }
                }
            }
            let mm = if deltas.is_empty() {
                String::new()
            } else {
                format!("C+m?,{};", deltas.join(","))
            };
            (mm, mls)
        };

        // Case 1: Vertebrate global-CpG regime with 38 kb (>2x full 16.5 kb circle!) rolling-circle
        // mtDNA reads + full-length 16.5 kb hypermethylated polymorphic mega-NUMT reads (~80% CG).
        let mut vert_records: Vec<(String, u16, Vec<u8>, String, Vec<u8>)> = Vec::new();
        for i in 0..24 {
            let start = (i * 1_100) % mito.len();
            let len = if i < 6 { 38_000 } else { 14_000 }; // 6 multi-lap rolling-circle ONT reads!
            let rseq: Vec<u8> = (0..len).map(|k| mito[(start + k) % mito.len()]).collect();
            let (mm, ml) = make_cpg_mm_ml(&rseq, 10);
            vert_records.push((format!("vert_mt_{i}"), 4, rseq, mm, ml));
        }
        for i in 0..8 {
            // Full-circle 16.5 kb concatenated mega-NUMT (zero nuclear overhang; only CpG catches it!)
            let rseq = mito.clone();
            let (mm, ml) = make_cpg_mm_ml(&rseq, 235);
            vert_records.push((format!("vert_meganumt_{i}"), 4, rseq, mm, ml));
        }
        let vert_bam = dir.join("vert.bam");
        let vert_refs: Vec<_> = vert_records
            .iter()
            .map(|(n, fl, sq, mm, ml)| (n.as_str(), *fl, sq.as_slice(), Some(mm.as_str()), Some(ml.as_slice())))
            .collect();
        write_synthetic_bam(&vert_bam, &vert_refs);

        let vert_rep = run_methylation(&MethylationArgs {
            target: Some(target_fa.clone()),
            gfa: None,
            bam: vec![vert_bam],
            out: dir.join("vert_out"),
            prob_threshold: 0.5,
            max_cg_methylation: "auto".into(),
            min_read_cpgs: 5,
            max_overhang: 300,
            kmer: 31,
            circular: CircularMode::Auto,
            lineage: crate::seeddb::Lineage::Metazoa,
            threads: Some(2),
        })
        .unwrap();
        assert_eq!(vert_rep["lineage"], json!("metazoa"));
        assert_eq!(vert_rep["epigenetic_regime"], json!("vertebrate_global_cpg"));
        assert_eq!(vert_rep["parameters"]["chg_rescue_enabled"], json!(false));
        assert_eq!(
            vert_rep["reads"]["clean_organelle_reads"],
            json!(24),
            "all 24 vertebrate mtDNA reads (including 38 kb >2-lap rolling-circle reads) must be kept"
        );
        assert_eq!(vert_rep["reads"]["numt_nupt_chimeric_reads"], json!(0));
        assert_eq!(vert_rep["reads"]["numt_nupt_hypermethylated_reads"], json!(8));

        // Case 2: DNMT-deficient unmethylated insect regime (e.g. Drosophila): an AT-rich read with
        // zero CG calls and noisy CHG calls must NOT be discarded by plant CHG rescue, and a sparse
        // 5-CpG read with 2 noisy calls (40% CG) must NOT be falsely pruned under unmethylated guard.
        let mut ins_records: Vec<(String, u16, Vec<u8>, String, Vec<u8>)> = Vec::new();
        for i in 0..20 {
            let start = (i * 800) % mito.len();
            let rseq: Vec<u8> = (0..6_000).map(|k| mito[(start + k) % mito.len()]).collect();
            let (mm, ml) = make_cpg_mm_ml(&rseq, 8);
            ins_records.push((format!("ins_mt_{i}"), 4, rseq, mm, ml));
        }
        // Read with 5 called CpGs where 2/5 (0.40) happen to have noisy basecaller prob > 0.5
        let sparse_seq = mito[0..1_200].to_vec();
        let mut deltas = Vec::new();
        let mut mls = Vec::new();
        let mut skip = 0usize;
        for i in 0..sparse_seq.len() {
            if sparse_seq[i] == b'C' {
                if i + 1 < sparse_seq.len() && sparse_seq[i + 1] == b'G' && mls.len() < 5 {
                    deltas.push(skip.to_string());
                    mls.push(if mls.len() < 2 { 200u8 } else { 10u8 }); // 2/5 = 0.40
                    skip = 0;
                } else {
                    skip += 1;
                }
            }
        }
        ins_records.push((
            "ins_noisy_sparse".into(),
            4,
            sparse_seq,
            format!("C+m?,{};", deltas.join(",")),
            mls,
        ));

        let ins_bam = dir.join("ins.bam");
        let ins_refs: Vec<_> = ins_records
            .iter()
            .map(|(n, fl, sq, mm, ml)| (n.as_str(), *fl, sq.as_slice(), Some(mm.as_str()), Some(ml.as_slice())))
            .collect();
        write_synthetic_bam(&ins_bam, &ins_refs);

        let ins_rep = run_methylation(&MethylationArgs {
            target: Some(target_fa),
            gfa: None,
            bam: vec![ins_bam],
            out: dir.join("ins_out"),
            prob_threshold: 0.5,
            max_cg_methylation: "auto".into(),
            min_read_cpgs: 5,
            max_overhang: 300,
            kmer: 31,
            circular: CircularMode::Auto,
            lineage: crate::seeddb::Lineage::Metazoa,
            threads: Some(2),
        })
        .unwrap();
        assert_eq!(ins_rep["epigenetic_regime"], json!("invertebrate_unmethylated"));
        assert_eq!(ins_rep["reads"]["clean_organelle_reads"], json!(21));
        assert_eq!(ins_rep["reads"]["numt_nupt_hypermethylated_reads"], json!(0));

        let _ = fs::remove_dir_all(&dir);
    }
}
