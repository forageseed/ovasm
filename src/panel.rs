//! Stage 1c: which organelle a set of reads comes from, by conserved proteins (`ovasm panel`).
//!
//! `ovasm discover` finds organelle reads by depth alone and does not say which depth cluster
//! is which organelle; seeds say so, but only for species close enough to the seeds to share
//! 21-mers. Protein sequence stays recognisable far beyond that. Reads are translated in six
//! frames; every stretch between stop codons of at least `min_orf` residues is looked up by
//! amino-acid k-mers in a panel of organelle proteins (`data/organelle_proteins.faa`:
//! mitochondrial protein references and the CDS translations of land-plant plastomes), and goes
//! to the panel gene it shares most k-mers with, if that is at least `min_hits`. Non-coding
//! sequence rarely holds a stop-free stretch that long (about 2% of the stretches of a random
//! sequence), and the panel's k-mers are about 1% of all 6-mers, so chance hits stay far
//! below `min_hits`. A set of reads is called by how many distinct genes of each class it
//! holds, so a mitochondrial cluster is told from a plastid one and from nuclear repeats
//! (which hold none) without any reference genome of the species or a close relative.

use std::borrow::Cow;
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::Serialize;

use crate::evidence::revcomp;
use crate::recruit::{peak_rss_mb, read_batches};

/// The panel shipped in the binary.
const PANEL: &str = include_str!("../data/organelle_proteins.faa");

/// Standard genetic code, codon index = 16 * b1 + 4 * b2 + b3 with A, C, G, T = 0..3.
const CODE: &[u8; 64] = b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
const STOP: u8 = 20;
const RESIDUES: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

#[derive(Debug, Clone, Serialize)]
pub struct PanelParams {
    /// Amino-acid k-mer length (at most 12).
    pub k: usize,
    /// Shortest stop-free stretch that is looked up.
    pub min_orf: usize,
    /// Shared k-mers a stretch needs with its best gene.
    pub min_hits: u32,
    /// Stretches a gene needs to count as present in a set of reads.
    pub min_gene_segments: u64,
    /// Distinct genes a class needs for a call.
    pub min_genes: usize,
    /// Stretches assigned to a class's genes per read it needs for a call. Organelle reads
    /// carry many (plastid 2-9 per read, mitochondrion 0.7-3 in the v3 clusters); the
    /// plastid-derived and mitochondrion-derived sequence scattered through the nuclear
    /// genome holds many different genes too, but a few reads each (0.03-0.16 per read).
    pub min_density: f64,
}

impl Default for PanelParams {
    fn default() -> Self {
        Self {
            k: 6,
            min_orf: 80,
            min_hits: 5,
            min_gene_segments: 2,
            min_genes: 5,
            min_density: 0.3,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Class {
    Mito,
    Plastid,
}

struct Gene {
    class: Class,
    name: String,
}

pub struct Panel {
    genes: Vec<Gene>,
    /// amino-acid k-mer -> genes (any species) holding it
    index: FxHashMap<u64, Vec<u16>>,
    k: usize,
    proteins: usize,
}

#[inline]
fn residue_code(b: u8) -> Option<u8> {
    RESIDUES
        .iter()
        .position(|&r| r == b.to_ascii_uppercase())
        .map(|i| i as u8)
}

#[inline]
fn base_code(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Residue code of the codon at the start of `s` (`STOP` for a stop codon, `None` when a base
/// is not A, C, G or T).
#[inline]
fn translate_codon(s: &[u8]) -> Option<u8> {
    let i = 16 * base_code(s[0])? + 4 * base_code(s[1])? + base_code(s[2])?;
    let aa = CODE[i];
    if aa == b'*' {
        Some(STOP)
    } else {
        residue_code(aa)
    }
}

impl Panel {
    /// Parse a panel FASTA with headers `>class|gene|species`; proteins of species whose name
    /// starts with an `exclude` prefix are left out (for testing against a genus, say).
    pub fn parse(text: &str, k: usize, exclude: &[String]) -> Result<Self> {
        if !(2..=12).contains(&k) {
            bail!("the amino-acid k must be between 2 and 12");
        }
        let mut genes: Vec<Gene> = Vec::new();
        let mut gene_ids: FxHashMap<(u8, String), u16> = FxHashMap::default();
        let mut index: FxHashMap<u64, Vec<u16>> = FxHashMap::default();
        let mut proteins = 0;
        let mut add = |header: &str, seq: &[u8]| -> Result<()> {
            let mut f = header.split('|');
            let (class, gene, species) = match (f.next(), f.next(), f.next()) {
                (Some(c), Some(g), Some(s)) => (c, g, s),
                _ => bail!("panel header {header:?} is not class|gene|species"),
            };
            if exclude.iter().any(|p| species.starts_with(p.as_str())) {
                return Ok(());
            }
            let class = match class {
                "mito" => Class::Mito,
                "plastid" => Class::Plastid,
                other => bail!("panel class {other:?} is neither mito nor plastid"),
            };
            let n = gene_ids.len();
            let id = *gene_ids
                .entry((class as u8, gene.to_string()))
                .or_insert_with(|| {
                    genes.push(Gene {
                        class,
                        name: gene.to_string(),
                    });
                    n as u16
                });
            proteins += 1;
            let codes: Vec<u8> = seq.iter().filter_map(|&b| residue_code(b)).collect();
            let mask = (1u64 << (5 * k)) - 1;
            let mut code = 0u64;
            for (i, &a) in codes.iter().enumerate() {
                code = ((code << 5) | a as u64) & mask;
                if i + 1 >= k {
                    let v = index.entry(code).or_default();
                    if !v.contains(&id) {
                        v.push(id);
                    }
                }
            }
            Ok(())
        };
        let (mut header, mut seq): (Option<String>, Vec<u8>) = (None, Vec::new());
        for line in text.lines() {
            if let Some(h) = line.strip_prefix('>') {
                if let Some(prev) = header.take() {
                    add(&prev, &seq)?;
                }
                header = Some(h.trim().to_string());
                seq.clear();
            } else {
                seq.extend_from_slice(line.trim().as_bytes());
            }
        }
        if let Some(prev) = header {
            add(&prev, &seq)?;
        }
        if genes.is_empty() {
            bail!("the protein panel is empty");
        }
        Ok(Self {
            genes,
            index,
            k,
            proteins,
        })
    }

    pub fn embedded(k: usize, exclude: &[String]) -> Result<Self> {
        Self::parse(PANEL, k, exclude)
    }

    /// Proteins and distinct genes in the panel.
    pub(crate) fn sizes(&self) -> (usize, usize) {
        (self.proteins, self.genes.len())
    }
}

/// What a set of reads holds: stretches assigned to each panel gene.
#[derive(Clone)]
struct Tally {
    gene_segments: Vec<u64>,
    reads: u64,
    bases: u64,
    orfs: u64,
    assigned: u64,
}

impl Tally {
    fn new(genes: usize) -> Self {
        Self {
            gene_segments: vec![0; genes],
            reads: 0,
            bases: 0,
            orfs: 0,
            assigned: 0,
        }
    }

    fn merge(mut self, other: Tally) -> Tally {
        for (a, b) in self.gene_segments.iter_mut().zip(other.gene_segments) {
            *a += b;
        }
        self.reads += other.reads;
        self.bases += other.bases;
        self.orfs += other.orfs;
        self.assigned += other.assigned;
        self
    }
}

/// Look up one stop-free stretch; assign it to its best gene.
fn assign(aa: &[u8], panel: &Panel, p: &PanelParams, hits: &mut [u32], tally: &mut Tally) {
    tally.orfs += 1;
    let mask = (1u64 << (5 * panel.k)) - 1;
    let mut code = 0u64;
    let mut touched: Vec<u16> = Vec::new();
    for (i, &a) in aa.iter().enumerate() {
        code = ((code << 5) | a as u64) & mask;
        if i + 1 < panel.k {
            continue;
        }
        if let Some(genes) = panel.index.get(&code) {
            for &g in genes {
                if hits[g as usize] == 0 {
                    touched.push(g);
                }
                hits[g as usize] += 1;
            }
        }
    }
    let mut best: Option<(u16, u32)> = None;
    for &g in &touched {
        let h = hits[g as usize];
        // the lower gene id wins a tie, so the call does not depend on hash order
        if best.map_or(true, |(bg, bh)| h > bh || (h == bh && g < bg)) {
            best = Some((g, h));
        }
        hits[g as usize] = 0;
    }
    if let Some((g, h)) = best {
        if h >= p.min_hits {
            tally.gene_segments[g as usize] += 1;
            tally.assigned += 1;
        }
    }
}

fn scan_read(seq: &[u8], panel: &Panel, p: &PanelParams, hits: &mut [u32], tally: &mut Tally) {
    tally.reads += 1;
    tally.bases += seq.len() as u64;
    for strand in 0..2 {
        let s: Cow<[u8]> = if strand == 0 {
            Cow::Borrowed(seq)
        } else {
            Cow::Owned(revcomp(seq))
        };
        for frame in 0..3 {
            let mut aa: Vec<u8> = Vec::with_capacity(s.len() / 3 + 1);
            let mut i = frame;
            while i + 3 <= s.len() {
                match translate_codon(&s[i..i + 3]) {
                    Some(a) if a < STOP => aa.push(a),
                    _ => {
                        if aa.len() >= p.min_orf {
                            assign(&aa, panel, p, hits, tally);
                        }
                        aa.clear();
                    }
                }
                i += 3;
            }
            if aa.len() >= p.min_orf {
                assign(&aa, panel, p, hits, tally);
            }
        }
    }
}

#[derive(Serialize)]
pub struct ClassReport {
    /// Distinct genes with at least `min_gene_segments` stretches.
    pub genes: usize,
    /// Stretches assigned to those genes.
    pub segments: u64,
    /// Gene and stretches, most first.
    pub gene_segments: Vec<(String, u64)>,
}

#[derive(Serialize)]
pub struct SetReport {
    pub file: String,
    pub reads: u64,
    pub bases: u64,
    /// Stop-free stretches of at least `min_orf` residues, six frames.
    pub orfs: u64,
    pub assigned_orfs: u64,
    pub mitochondrion: ClassReport,
    pub plastid: ClassReport,
    /// Gene stretches per read, by class (see `PanelParams::min_density`).
    pub mitochondrion_density: f64,
    pub plastid_density: f64,
    /// "mitochondrion", "plastid", "mixed" (both classes present) or "none".
    pub call: &'static str,
}

#[derive(Serialize)]
pub struct PanelReport {
    pub tool: &'static str,
    pub version: &'static str,
    pub params: PanelParams,
    pub panel_proteins: usize,
    pub panel_genes: usize,
    pub excluded_species: Vec<String>,
    pub sets: Vec<SetReport>,
    pub elapsed_seconds: f64,
    pub peak_rss_mb: Option<f64>,
}

fn class_report(panel: &Panel, tally: &Tally, class: Class, p: &PanelParams) -> ClassReport {
    let mut gene_segments: Vec<(String, u64)> = panel
        .genes
        .iter()
        .zip(&tally.gene_segments)
        .filter(|(g, &n)| g.class == class && n >= p.min_gene_segments)
        .map(|(g, &n)| (g.name.clone(), n))
        .collect();
    gene_segments.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    ClassReport {
        genes: gene_segments.len(),
        segments: gene_segments.iter().map(|x| x.1).sum(),
        gene_segments,
    }
}

/// Panel genes in assembled sequences (the segments of one graph component, say), scanned
/// like reads: (mitochondrial, plastid), a gene counting from `min_gene_segments` stretches.
pub(crate) fn sequence_genes<'a>(
    seqs: impl IntoIterator<Item = &'a [u8]>,
    panel: &Panel,
    p: &PanelParams,
) -> (ClassReport, ClassReport) {
    let mut tally = Tally::new(panel.genes.len());
    let mut hits = vec![0u32; panel.genes.len()];
    for seq in seqs {
        scan_read(seq, panel, p, &mut hits, &mut tally);
    }
    (
        class_report(panel, &tally, Class::Mito, p),
        class_report(panel, &tally, Class::Plastid, p),
    )
}

/// The call from the two classes' (genes, stretches per read): a class counts when it has at
/// least `min_genes` genes at `min_density`; a counting class with at least twice the other's
/// genes (or alone) is the call, both counting without that margin is "mixed".
pub fn call_from(mito: (usize, f64), plastid: (usize, f64), p: &PanelParams) -> &'static str {
    let counts = |c: (usize, f64)| c.0 >= p.min_genes && c.1 >= p.min_density;
    let (m_ok, p_ok) = (counts(mito), counts(plastid));
    if p_ok && (!m_ok || plastid.0 >= 2 * mito.0) {
        "plastid"
    } else if m_ok && (!p_ok || mito.0 >= 2 * plastid.0) {
        "mitochondrion"
    } else if m_ok && p_ok {
        "mixed"
    } else {
        "none"
    }
}

pub fn panel_set(path: &Path, panel: &Panel, p: &PanelParams) -> Result<SetReport> {
    let (rx, reader) = read_batches(&[path.to_path_buf()]);
    let mut total = Tally::new(panel.genes.len());
    for batch in rx.iter() {
        let t = batch
            .par_iter()
            .fold(
                || (Tally::new(panel.genes.len()), vec![0u32; panel.genes.len()]),
                |(mut tally, mut hits), r| {
                    scan_read(&r.seq, panel, p, &mut hits, &mut tally);
                    (tally, hits)
                },
            )
            .map(|(t, _)| t)
            .reduce(|| Tally::new(panel.genes.len()), Tally::merge);
        total = total.merge(t);
    }
    reader.join().expect("reader thread panicked")?;
    let mitochondrion = class_report(panel, &total, Class::Mito, p);
    let plastid = class_report(panel, &total, Class::Plastid, p);
    let density = |c: &ClassReport| c.segments as f64 / total.reads.max(1) as f64;
    let (mitochondrion_density, plastid_density) = (density(&mitochondrion), density(&plastid));
    let call = call_from(
        (mitochondrion.genes, mitochondrion_density),
        (plastid.genes, plastid_density),
        p,
    );
    Ok(SetReport {
        file: path.display().to_string(),
        reads: total.reads,
        bases: total.bases,
        orfs: total.orfs,
        assigned_orfs: total.assigned,
        mitochondrion,
        plastid,
        mitochondrion_density,
        plastid_density,
        call,
    })
}

pub fn run(
    sets: &[PathBuf],
    panel_file: Option<&Path>,
    p: &PanelParams,
    exclude: &[String],
    out_json: &Path,
) -> Result<PanelReport> {
    let t0 = Instant::now();
    let panel = match panel_file {
        Some(f) => Panel::parse(
            &fs::read_to_string(f).with_context(|| format!("cannot read {}", f.display()))?,
            p.k,
            exclude,
        )?,
        None => Panel::embedded(p.k, exclude)?,
    };
    let mut reports = Vec::new();
    for s in sets {
        reports.push(panel_set(s, &panel, p)?);
    }
    let report = PanelReport {
        tool: "ovasm panel",
        version: env!("CARGO_PKG_VERSION"),
        params: p.clone(),
        panel_proteins: panel.proteins,
        panel_genes: panel.genes.len(),
        excluded_species: exclude.to_vec(),
        sets: reports,
        elapsed_seconds: t0.elapsed().as_secs_f64(),
        peak_rss_mb: peak_rss_mb(),
    };
    fs::write(out_json, serde_json::to_string_pretty(&report)?)?;
    Ok(report)
}

/// One line per set for the terminal.
pub fn summary(r: &PanelReport) -> String {
    let mut s = String::new();
    for x in &r.sets {
        let _ = writeln!(
            s,
            "[ovasm] {}: {} reads, {} stretches, {} assigned; mitochondrion {} genes ({:.2}/read), plastid {} genes ({:.2}/read) -> {}",
            x.file,
            x.reads,
            x.orfs,
            x.assigned_orfs,
            x.mitochondrion.genes,
            x.mitochondrion_density,
            x.plastid.genes,
            x.plastid_density,
            x.call
        );
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A deterministic protein of `n` residues.
    fn protein(n: usize, seed: u64) -> Vec<u8> {
        let mut s = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        (0..n)
            .map(|_| {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                RESIDUES[(s % 20) as usize]
            })
            .collect()
    }

    /// DNA that translates to `prot` (the first codon for each residue), then a stop.
    fn back_translate(prot: &[u8]) -> Vec<u8> {
        let bases = b"ACGT";
        let mut dna = Vec::new();
        for &a in prot {
            let i = CODE.iter().position(|&c| c == a).unwrap();
            dna.extend([bases[i / 16], bases[(i / 4) % 4], bases[i % 4]]);
        }
        dna.extend(b"TAA");
        dna
    }

    fn panel_text() -> String {
        let mut t = String::new();
        for (class, gene, seed) in [
            ("mito", "cox1", 1),
            ("mito", "atp1", 2),
            ("plastid", "rbcl", 3),
        ] {
            // two "species" per gene: the second differs at every tenth residue
            let a = protein(200, seed);
            let mut b = a.clone();
            for i in (0..b.len()).step_by(10) {
                b[i] = if b[i] == b'A' { b'C' } else { b'A' };
            }
            let _ = writeln!(t, ">{class}|{gene}|Sp_one\n{}", String::from_utf8_lossy(&a));
            let _ = writeln!(
                t,
                ">{class}|{gene}|Other_two\n{}",
                String::from_utf8_lossy(&b)
            );
        }
        t
    }

    #[test]
    fn translation_reads_both_strands_and_stops() {
        let prot = protein(100, 9);
        let dna = back_translate(&prot);
        let mut aa = Vec::new();
        for c in dna.chunks(3) {
            match translate_codon(c) {
                Some(a) if a < STOP => aa.push(RESIDUES[a as usize]),
                Some(_) => break,
                None => panic!("invalid codon"),
            }
        }
        assert_eq!(aa, prot);
        assert_eq!(translate_codon(b"ANA"), None);
        assert_eq!(translate_codon(b"TGA"), Some(STOP));
        assert_eq!(translate_codon(b"atg"), residue_code(b'M'));
    }

    #[test]
    fn a_read_goes_to_the_gene_it_shares_most_kmers_with() {
        let panel = Panel::parse(&panel_text(), 6, &[]).unwrap();
        let p = PanelParams::default();
        let mut hits = vec![0u32; panel.genes.len()];
        let gene_of = |name: &str, class: Class| {
            panel
                .genes
                .iter()
                .position(|g| g.name == name && g.class == class)
                .unwrap()
        };
        // a gene on the forward strand, with non-coding flanks
        let mut read = b"GCGCGCATATATGGCCATTAGGCGAAT".to_vec();
        read.extend(back_translate(&protein(200, 1)));
        read.extend(b"TTGACCATGGCATTAGC");
        let mut t = Tally::new(panel.genes.len());
        scan_read(&read, &panel, &p, &mut hits, &mut t);
        assert_eq!(t.gene_segments[gene_of("cox1", Class::Mito)], 1);
        assert_eq!(t.assigned, 1);
        // the same gene on the reverse strand
        let mut t = Tally::new(panel.genes.len());
        scan_read(&revcomp(&read), &panel, &p, &mut hits, &mut t);
        assert_eq!(t.gene_segments[gene_of("cox1", Class::Mito)], 1);
        // a protein the panel does not hold, and short ORFs, are not assigned
        let mut t = Tally::new(panel.genes.len());
        scan_read(
            &back_translate(&protein(200, 77)),
            &panel,
            &p,
            &mut hits,
            &mut t,
        );
        scan_read(
            &back_translate(&protein(60, 1)),
            &panel,
            &p,
            &mut hits,
            &mut t,
        );
        assert_eq!(t.assigned, 0);
        assert!(hits.iter().all(|&h| h == 0), "scratch counts are reset");
    }

    #[test]
    fn species_can_be_left_out_of_the_panel() {
        let all = Panel::parse(&panel_text(), 6, &[]).unwrap();
        let some = Panel::parse(&panel_text(), 6, &["Other".to_string()]).unwrap();
        assert_eq!((all.proteins, some.proteins), (6, 3));
        assert!(some.index.len() < all.index.len());
    }

    #[test]
    fn calls_need_enough_genes_density_and_a_margin() {
        let p = PanelParams::default();
        assert_eq!(call_from((0, 0.0), (30, 4.0), &p), "plastid");
        assert_eq!(call_from((12, 1.2), (3, 0.1), &p), "mitochondrion");
        // plastid DNA inside a mitochondrial set: a few plastid genes do not make it mixed
        assert_eq!(call_from((15, 1.0), (6, 0.5), &p), "mitochondrion");
        assert_eq!(call_from((8, 1.0), (6, 0.5), &p), "mixed");
        assert_eq!(call_from((2, 0.9), (1, 0.2), &p), "none");
        assert_eq!(call_from((4, 0.9), (0, 0.0), &p), "none");
        // many genes at a few reads each (plastid-derived nuclear DNA) is not a plastid
        assert_eq!(call_from((7, 0.1), (36, 0.12), &p), "none");
        assert_eq!(call_from((0, 0.0), (7, 0.15), &p), "none");
        // a plastid-like class too thin to count does not make a mitochondrial set mixed
        assert_eq!(call_from((37, 0.7), (13, 0.2), &p), "mitochondrion");
    }

    #[test]
    fn the_embedded_panel_parses_with_both_classes() {
        let panel = Panel::embedded(6, &[]).unwrap();
        assert!(panel.proteins > 2000);
        assert!(panel
            .genes
            .iter()
            .any(|g| g.class == Class::Mito && g.name == "cox1"));
        assert!(panel
            .genes
            .iter()
            .any(|g| g.class == Class::Plastid && g.name == "rbcl"));
    }

    #[test]
    fn a_set_of_reads_is_called_from_its_genes() {
        let dir = std::env::temp_dir().join(format!("ovasm-identify-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let mut fa = String::new();
        for (i, seed) in [(0, 1u64), (1, 2), (2, 1), (3, 2), (4, 1), (5, 2)] {
            let mut read = b"ATATGGCCGGTTAACCTTGGAACCTTAAGG".to_vec();
            read.extend(back_translate(&protein(200, seed)));
            let _ = writeln!(fa, ">r{i}\n{}", String::from_utf8_lossy(&read));
        }
        fs::write(dir.join("mito.fa"), &fa).unwrap();
        let panel = Panel::parse(&panel_text(), 6, &[]).unwrap();
        let p = PanelParams {
            min_genes: 2,
            ..PanelParams::default()
        };
        let r = panel_set(&dir.join("mito.fa"), &panel, &p).unwrap();
        assert_eq!(r.reads, 6);
        assert_eq!(
            r.mitochondrion.genes, 2,
            "{:?}",
            r.mitochondrion.gene_segments
        );
        assert_eq!(r.plastid.genes, 0);
        assert_eq!(r.call, "mitochondrion");
        fs::remove_dir_all(&dir).ok();
    }
}
