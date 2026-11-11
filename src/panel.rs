//! Stage 1c: which organelle a set of reads comes from, by conserved proteins (`ovasm find classify`).
//!
//! `ovasm find discover` finds organelle reads by depth alone and does not say which depth cluster
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
/// Mitochondrial core proteins of metazoa (551) and fungi (696), for protein seeds: from the RefSeq
/// mitochondrion release, chosen gene by gene as the fewest proteins that find held-out genera
/// (leave-genus-out greedy cover, scripts/polish_panel.py; no genus of the PMAT2 benchmark, see
/// data/lineage_panels.excluded_genera.txt). Larger covers find barely more and hit random
/// sequence far more often (see `Panel::lineage_min_hits`).
const METAZOA_PANEL: &str = include_str!("../data/metazoa_mito_proteins.faa");
const FUNGI_PANEL: &str = include_str!("../data/fungi_mito_proteins.faa");

/// Standard genetic code, codon index = 16 * b1 + 4 * b2 + b3 with A, C, G, T = 0..3.
const CODE: &[u8; 64] = b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";

/// An NCBI translation table as a codon table in `CODE`'s order: 1 (standard; land-plant
/// organelles), 4 (mold mitochondrial: TGA is Trp; most fungi), 5 (invertebrate mitochondrial: TGA
/// Trp, ATA Met, AGA/AGG Ser) or 2 (vertebrate mitochondrial: TGA Trp, ATA Met, AGA/AGG stop).
/// Mitochondrial genes of animals and fungi read TGA as Trp, so the standard table cuts them into
/// pieces shorter than `min_orf`.
pub fn genetic_code(table: u8) -> Result<[u8; 64]> {
    let mut c = *CODE;
    let at = |s: &[u8; 3]| {
        let b = |x: u8| b"ACGT".iter().position(|&y| y == x).unwrap();
        16 * b(s[0]) + 4 * b(s[1]) + b(s[2])
    };
    match table {
        1 => {}
        2 | 4 | 5 => {
            c[at(b"TGA")] = b'W';
            if table != 4 {
                c[at(b"ATA")] = b'M';
                let aga = if table == 5 { b'S' } else { b'*' };
                c[at(b"AGA")] = aga;
                c[at(b"AGG")] = aga;
            }
        }
        other => bail!("genetic code {other} is not supported (1, 2, 4 or 5)"),
    }
    Ok(c)
}
const STOP: u8 = 20;
const RESIDUES: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

/// BLOSUM62 in `RESIDUES` order, for the score of a gene stretch against its panel proteins.
#[rustfmt::skip]
const BLOSUM62: [[i8; 20]; 20] = [
    [ 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2], // A
    [ 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2], // C
    [-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3], // D
    [-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2], // E
    [-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3], // F
    [ 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3], // G
    [-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2], // H
    [-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1], // I
    [-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2], // K
    [-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1], // L
    [-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1], // M
    [-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2], // N
    [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3], // P
    [-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1], // Q
    [-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2], // R
    [ 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2], // S
    [ 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2], // T
    [ 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1], // V
    [-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2], // W
    [-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7], // Y
];

/// Karlin-Altschul parameters of ungapped BLOSUM62 alignments (BLAST's lambda and K for the
/// matrix with standard residue frequencies): bits = (LAMBDA * score - ln K) / ln 2.
const UNGAPPED_LAMBDA: f64 = 0.3176;
const UNGAPPED_K: f64 = 0.134;
/// X-drop of the ungapped extension, in raw BLOSUM62 score (about 7 bits, as BLAST's).
const XDROP: i32 = 16;
/// Chance alignments a gene stretch may expect at the score it reached against the whole panel
/// (`Panel::max_evalue`), animal and fungal panels. Shared amino-acid 6-mers alone do not say a
/// stretch is homologous: the hydrophobic, low-complexity peptides that AT-rich nuclear and
/// bacterial sequence translates to share many 6-mers with the short membrane proteins (atp8,
/// nad4l, nad6, nad3) by chance. PMAT2 benchmark, read sets without the mitochondrion: none of the
/// 31 whole-read seeds of Calliactis parasitica and 3 of the 46 of Myriapora truncata had any
/// homology to the panel (DIAMOND blastx, e 1e-5); they became a 28 kb and a 9.7 kb
/// "mitochondrion". At 1e-6 a stretch of 300 residues needs about 45 bits after what its
/// composition alone scores is taken off (`significant_bits`): of those seeds' 403 assigned
/// stretches the 12 of the 3 real homologs keep their gene (up to 220 bits), 348 fall short
/// forward and 43 more reach it forward (up to 120 bits) but score as much reversed; of the 84
/// assigned stretches of seven reference mitochondria (two sponges, a sea anemone, two bryozoans,
/// a moth, a sawfly; codes 4 and 5) 83 keep their gene (60-885 bits forward, at most 35 reversed),
/// the one lost a sawfly nad6 at 20 bits.
pub const LINEAGE_MAX_EVALUE: f64 = 1e-6;

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
    /// Codon table reads are translated with (`genetic_code`).
    code: [u8; 64],
    /// Each protein's gene and residue codes, and amino-acid k-mer -> (protein, position): the
    /// seeds of a stretch's alignment to its gene (`stretch_bits`; built by `with_max_evalue`,
    /// k-mers alone do not need it).
    protein_gene: Vec<u16>,
    protein_seqs: Vec<Vec<u8>>,
    occurrences: FxHashMap<u64, Vec<(u32, u32)>>,
    /// Residues in the panel: the database size of a stretch's E-value.
    residues: u64,
    /// A stretch keeps its gene only when its best ungapped alignment to that gene's proteins has
    /// an E-value against the whole panel of at most this (`LINEAGE_MAX_EVALUE`); None: shared
    /// k-mers alone decide, as for the land-plant panel.
    max_evalue: Option<f64>,
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
#[cfg(test)]
fn translate_codon(s: &[u8]) -> Option<u8> {
    translate_with(CODE, s)
}

/// `translate_codon` with the codon table `code`.
#[inline]
fn translate_with(code: &[u8; 64], s: &[u8]) -> Option<u8> {
    let i = 16 * base_code(s[0])? + 4 * base_code(s[1])? + base_code(s[2])?;
    let aa = code[i];
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
        let mut protein_gene: Vec<u16> = Vec::new();
        let mut protein_seqs: Vec<Vec<u8>> = Vec::new();
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
            protein_gene.push(id);
            protein_seqs.push(codes);
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
        let residues = protein_seqs.iter().map(|s| s.len() as u64).sum();
        Ok(Self {
            genes,
            index,
            k,
            proteins,
            code: *CODE,
            protein_gene,
            protein_seqs,
            occurrences: FxHashMap::default(),
            residues,
            max_evalue: None,
        })
    }

    pub fn embedded(k: usize, exclude: &[String]) -> Result<Self> {
        Self::parse(PANEL, k, exclude)
    }

    /// The built-in panel of a seed lineage with its genetic code: "landplants" (mitochondrion
    /// and plastid, standard code), "metazoa" (mitochondrion, code 5) or "fungi" (mitochondrion,
    /// code 4).
    pub fn for_lineage(lineage: &str, k: usize, exclude: &[String]) -> Result<Self> {
        let (text, table) = match lineage {
            "landplants" => (PANEL, 1),
            "metazoa" => (METAZOA_PANEL, 5),
            "fungi" => (FUNGI_PANEL, 4),
            other => bail!("no built-in protein panel for lineage {other:?}"),
        };
        Ok(Self::parse(text, k, exclude)?
            .with_code(table)?
            .with_max_evalue(Self::lineage_max_evalue(lineage)))
    }

    /// `Panel::max_evalue` of a lineage's panel: `LINEAGE_MAX_EVALUE` for the animal and fungal
    /// panels; None (shared k-mers decide) for the land-plant one, whose `min_hits` was set on
    /// random sequence and checked on five mosses (2026-11-05).
    pub fn lineage_max_evalue(lineage: &str) -> Option<f64> {
        match lineage {
            "metazoa" | "fungi" => Some(LINEAGE_MAX_EVALUE),
            _ => None,
        }
    }

    /// The same panel, a stretch keeping its gene only at an alignment E-value of at most `e`.
    pub fn with_max_evalue(mut self, e: Option<f64>) -> Self {
        self.max_evalue = e;
        if e.is_some() && self.occurrences.is_empty() {
            let mask = (1u64 << (5 * self.k)) - 1;
            for (pid, codes) in self.protein_seqs.iter().enumerate() {
                let mut code = 0u64;
                for (i, &a) in codes.iter().enumerate() {
                    code = ((code << 5) | a as u64) & mask;
                    if i + 1 >= self.k {
                        self.occurrences
                            .entry(code)
                            .or_default()
                            .push((pid as u32, (i + 1 - self.k) as u32));
                    }
                }
            }
        }
        self
    }

    /// `Panel::max_evalue` (for reports).
    pub fn max_evalue(&self) -> Option<f64> {
        self.max_evalue
    }

    /// `PanelParams::min_hits` for a lineage's built-in panel: the lowest threshold at which
    /// random DNA (six frames, 10 Mb, GC 0.30 and 0.45) gets at most 0.2 gene stretches per Mb (the
    /// land-plant panel: 0.4 at 5 hits), which keeps the most held-out genes
    /// (scripts/panel_tradeoff.py, 2026-11-05). Animal and fungal mitochondrial
    /// proteins are rich in the hydrophobic residues AT-rich random sequence translates to, so a
    /// larger panel needs a higher threshold and ends up finding fewer genes: metazoa 551 proteins
    /// at 5 hits found 98.96% of held-out genes, 13,868 proteins at 15 hits 98.08%; fungi 696 at 6
    /// hits 99.34%, 8,650 at 10 hits 99.10%.
    pub fn lineage_min_hits(lineage: &str) -> u32 {
        match lineage {
            "fungi" => 6,
            _ => 5,
        }
    }

    /// The same panel, reads translated with NCBI table `table` (`genetic_code`).
    pub fn with_code(mut self, table: u8) -> Result<Self> {
        self.code = genetic_code(table)?;
        Ok(self)
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

/// Look up one stop-free stretch: its best gene, if it shares at least `min_hits` k-mers with it.
fn best_gene(aa: &[u8], panel: &Panel, p: &PanelParams, hits: &mut [u32]) -> Option<u16> {
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
    let (g, h) = best?;
    if h < p.min_hits {
        return None;
    }
    match panel.max_evalue {
        None => Some(g),
        Some(e) => (significant_bits(aa, g, panel) >= bits_needed(aa.len(), panel.residues, e))
            .then_some(g),
    }
}

/// Bits a stretch of `m` residues needs against a panel of `n` residues for an E-value of at most
/// `e`: E = m n 2^-bits.
fn bits_needed(m: usize, n: u64, e: f64) -> f64 {
    ((m as f64) * (n as f64) / e).log2()
}

/// The stretch's best alignment to gene `g` in bits, less what its amino-acid composition alone
/// scores: the stretch reversed (residues and composition kept, homology gone) scores about
/// log2(m n) bits against unrelated sequence, and whatever it scores beyond that is the
/// composition's. Hydrophobic, low-complexity peptides score as high reversed as forward (43
/// stretches of the PMAT2 benchmark's false whole-read seeds reached the bits needed forward, up
/// to 120, and as much reversed); real genes at most 35 bits reversed (`LINEAGE_MAX_EVALUE`).
fn significant_bits(aa: &[u8], g: u16, panel: &Panel) -> f64 {
    let forward = stretch_bits(aa, g, panel);
    let reversed: Vec<u8> = aa.iter().rev().copied().collect();
    let chance = ((aa.len() as f64) * (panel.residues as f64)).log2();
    forward - (stretch_bits(&reversed, g, panel) - chance).max(0.0)
}

/// Score of the best ungapped alignment of `aa` to a protein of gene `g`, in bits: every distinct
/// diagonal of an amino-acid k-mer `aa` shares with a protein of `g` (the panel's k, the words
/// that picked the gene; stretch position minus protein position) is extended both ways from that
/// k-mer under BLOSUM62 with an X-drop (`XDROP`), as in BLAST's ungapped extension from word hits.
fn stretch_bits(aa: &[u8], g: u16, panel: &Panel) -> f64 {
    let k = panel.k;
    if aa.len() < k {
        return 0.0;
    }
    let mask = (1u64 << (5 * k)) - 1;
    let mut code = 0u64;
    let mut done: Vec<(u32, i64)> = Vec::new();
    let mut best = 0i32;
    for (i, &a) in aa.iter().enumerate() {
        code = ((code << 5) | a as u64) & mask;
        if i + 1 < k {
            continue;
        }
        let Some(occ) = panel.occurrences.get(&code) else {
            continue;
        };
        let qi = i + 1 - k;
        for &(pid, pj) in occ {
            if panel.protein_gene[pid as usize] != g {
                continue;
            }
            let d = qi as i64 - pj as i64;
            if done.contains(&(pid, d)) {
                continue;
            }
            done.push((pid, d));
            let s = ungapped_score(aa, &panel.protein_seqs[pid as usize], qi, pj as usize);
            best = best.max(s);
        }
    }
    if best <= 0 {
        return 0.0;
    }
    (UNGAPPED_LAMBDA * best as f64 - UNGAPPED_K.ln()) / std::f64::consts::LN_2
}

/// Best ungapped BLOSUM62 score of `q` against `s` through residue `qi` of `q` aligned to residue
/// `sj` of `s`: extended right from there and left from just before it, each side stopping once it
/// falls `XDROP` below its best.
fn ungapped_score(q: &[u8], s: &[u8], qi: usize, sj: usize) -> i32 {
    let score = |a: u8, b: u8| BLOSUM62[a as usize][b as usize] as i32;
    let (mut run, mut right) = (0i32, 0i32);
    let (mut i, mut j) = (qi, sj);
    while i < q.len() && j < s.len() {
        run += score(q[i], s[j]);
        right = right.max(run);
        if right - run > XDROP {
            break;
        }
        i += 1;
        j += 1;
    }
    let (mut run, mut left) = (0i32, 0i32);
    let (mut i, mut j) = (qi, sj);
    while i > 0 && j > 0 {
        i -= 1;
        j -= 1;
        run += score(q[i], s[j]);
        left = left.max(run);
        if left - run > XDROP {
            break;
        }
    }
    right + left
}

/// Every stop-free stretch of at least `min_orf` residues in the six frames of `seq`, with the
/// span it covers on `seq`'s forward strand: `f(residues, start, end)`.
fn for_each_stretch(
    seq: &[u8],
    code: &[u8; 64],
    min_orf: usize,
    mut f: impl FnMut(&[u8], usize, usize),
) {
    let n = seq.len();
    for strand in 0..2 {
        let s: Cow<[u8]> = if strand == 0 {
            Cow::Borrowed(seq)
        } else {
            Cow::Owned(revcomp(seq))
        };
        let mut emit = |aa: &[u8], start: usize, end: usize| {
            if strand == 0 {
                f(aa, start, end)
            } else {
                f(aa, n - end, n - start)
            }
        };
        for frame in 0..3 {
            let mut aa: Vec<u8> = Vec::with_capacity(s.len() / 3 + 1);
            let (mut i, mut start) = (frame, frame);
            while i + 3 <= s.len() {
                match translate_with(code, &s[i..i + 3]) {
                    Some(a) if a < STOP => {
                        if aa.is_empty() {
                            start = i;
                        }
                        aa.push(a)
                    }
                    _ => {
                        if aa.len() >= min_orf {
                            emit(&aa, start, i);
                        }
                        aa.clear();
                    }
                }
                i += 3;
            }
            if aa.len() >= min_orf {
                emit(&aa, start, i);
            }
        }
    }
}

fn scan_read(seq: &[u8], panel: &Panel, p: &PanelParams, hits: &mut [u32], tally: &mut Tally) {
    tally.reads += 1;
    tally.bases += seq.len() as u64;
    for_each_stretch(seq, &panel.code, p.min_orf, |aa, _, _| {
        tally.orfs += 1;
        if let Some(g) = best_gene(aa, panel, p, hits) {
            tally.gene_segments[g as usize] += 1;
            tally.assigned += 1;
        }
    });
}

/// The stretches of `seq` assigned to a panel gene: (gene, start, end) on the forward strand,
/// for cropping a read to its genes (`protein_seed`).
pub(crate) fn gene_stretches(
    seq: &[u8],
    panel: &Panel,
    p: &PanelParams,
    hits: &mut [u32],
) -> Vec<(u16, usize, usize)> {
    let mut out = Vec::new();
    for_each_stretch(seq, &panel.code, p.min_orf, |aa, start, end| {
        if let Some(g) = best_gene(aa, panel, p, hits) {
            out.push((g, start, end));
        }
    });
    out
}

impl Panel {
    /// Distinct panel genes (the length of a per-gene count vector).
    pub(crate) fn gene_count(&self) -> usize {
        self.genes.len()
    }

    /// Panel genes of an organelle ("mitochondrion" or "plastid"): 42 and 70 in the land-plant
    /// panel, 14 and 0 in the metazoa one.
    pub(crate) fn organelle_gene_count(&self, organelle: &str) -> usize {
        (0..self.genes.len())
            .filter(|&g| self.gene_organelle(g as u16) == organelle)
            .count()
    }

    /// A gene's organelle: "mitochondrion" or "plastid".
    pub(crate) fn gene_organelle(&self, g: u16) -> &'static str {
        match self.genes[g as usize].class {
            Class::Mito => "mitochondrion",
            Class::Plastid => "plastid",
        }
    }
}

/// The call of a set of `reads` reads holding `gene_segments` stretches per panel gene, as
/// `panel_set` makes it: (mitochondrion, plastid, their densities, call).
pub(crate) fn call_segments(
    panel: &Panel,
    gene_segments: Vec<u64>,
    reads: u64,
    p: &PanelParams,
) -> (ClassReport, ClassReport, f64, f64, &'static str) {
    let tally = Tally {
        gene_segments,
        reads,
        ..Tally::new(0)
    };
    let mitochondrion = class_report(panel, &tally, Class::Mito, p);
    let plastid = class_report(panel, &tally, Class::Plastid, p);
    let density = |c: &ClassReport| c.segments as f64 / reads.max(1) as f64;
    let (md, pd) = (density(&mitochondrion), density(&plastid));
    let call = call_from((mitochondrion.genes, md), (plastid.genes, pd), p);
    (mitochondrion, plastid, md, pd, call)
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
    /// The alignment E-value a gene stretch needs (`Panel::max_evalue`); null: k-mers alone.
    pub panel_max_evalue: Option<f64>,
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

/// `ovasm find classify`: each set called with `panel` when given (a lineage's panel), else with
/// `panel_file` or the embedded panel.
pub fn run(
    sets: &[PathBuf],
    panel_file: Option<&Path>,
    panel: Option<Panel>,
    p: &PanelParams,
    exclude: &[String],
    out_json: &Path,
) -> Result<PanelReport> {
    let t0 = Instant::now();
    let panel = match (panel, panel_file) {
        (Some(panel), _) => panel,
        (None, Some(f)) => Panel::parse(
            &fs::read_to_string(f).with_context(|| format!("cannot read {}", f.display()))?,
            p.k,
            exclude,
        )?,
        (None, None) => Panel::embedded(p.k, exclude)?,
    };
    let mut reports = Vec::new();
    for s in sets {
        reports.push(panel_set(s, &panel, p)?);
    }
    let report = PanelReport {
        tool: "ovasm find classify",
        version: env!("CARGO_PKG_VERSION"),
        params: p.clone(),
        panel_proteins: panel.proteins,
        panel_genes: panel.genes.len(),
        panel_max_evalue: panel.max_evalue,
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
    fn gene_stretches_are_spans_on_the_forward_strand() {
        let panel = Panel::parse(&panel_text(), 6, &[]).unwrap();
        let p = PanelParams::default();
        let mut hits = vec![0u32; panel.genes.len()];
        let cox1 = panel
            .genes
            .iter()
            .position(|g| g.name == "cox1")
            .unwrap() as u16;
        // frame 0 of the flank stops at TAG (bases 18-20): the stretch starts at 21 and runs to
        // the gene's stop codon (27 + 600)
        let mut read = b"GCGCGCATATATGGCCATTAGGCGAAT".to_vec();
        read.extend(back_translate(&protein(200, 1)));
        read.extend(b"TTGACCATGGCATTAGC");
        assert_eq!(
            gene_stretches(&read, &panel, &p, &mut hits),
            vec![(cox1, 21, 627)]
        );
        assert_eq!(
            gene_stretches(&revcomp(&read), &panel, &p, &mut hits),
            vec![(cox1, read.len() - 627, read.len() - 21)]
        );
    }

    #[test]
    fn mitochondrial_codes_read_tga_as_tryptophan() {
        let c = |t: u8, codon: &[u8; 3]| translate_with(&genetic_code(t).unwrap(), codon);
        assert_eq!(c(1, b"TGA"), Some(STOP));
        for t in [2, 4, 5] {
            assert_eq!(c(t, b"TGA"), residue_code(b'W'), "table {t}");
        }
        assert_eq!(c(4, b"ATA"), residue_code(b'I'));
        assert_eq!(c(5, b"ATA"), residue_code(b'M'));
        assert_eq!(c(5, b"AGG"), residue_code(b'S'));
        assert_eq!(c(2, b"AGA"), Some(STOP));
        assert_eq!(c(1, b"AGA"), residue_code(b'R'));
        assert!(genetic_code(3).is_err());
    }

    #[test]
    fn an_animal_gene_is_one_stretch_only_in_its_own_code() {
        let code5 = genetic_code(5).unwrap();
        // DNA that code 5 translates to the protein (its first codon per residue: TGA for Trp)
        let prot = protein(200, 1);
        let mut dna = Vec::new();
        for &a in &prot {
            let i = code5.iter().position(|&c| c == a).unwrap();
            dna.extend([b"ACGT"[i / 16], b"ACGT"[(i / 4) % 4], b"ACGT"[i % 4]]);
        }
        dna.extend(b"TAA");
        assert!(dna.windows(3).step_by(3).any(|w| w == b"TGA"));
        let p = PanelParams::default();
        let standard = Panel::parse(&panel_text(), 6, &[]).unwrap();
        let animal = Panel::parse(&panel_text(), 6, &[])
            .unwrap()
            .with_code(5)
            .unwrap();
        let mut hits = vec![0u32; standard.genes.len()];
        let found = gene_stretches(&dna, &animal, &p, &mut hits);
        assert_eq!(found.len(), 1);
        assert_eq!((found[0].1, found[0].2), (0, 600));
        assert!(gene_stretches(&dna, &standard, &p, &mut hits)
            .iter()
            .all(|&(_, a, b)| b - a < 600));
    }

    /// A panel of one mitochondrial gene whose protein is a random protein with a hydrophobic,
    /// low-complexity stretch (a transmembrane-like run) in its middle.
    fn membrane_panel() -> (String, Vec<u8>) {
        let tm: &[u8] = b"LILLFLLIVLFLILLFILLVLFLILLFLIVLLFILLFLIL";
        let mut prot = protein(120, 41);
        prot.extend_from_slice(tm);
        prot.extend(protein(120, 42));
        let text = format!(">mito|nad6|Sp_one\n{}\n", String::from_utf8_lossy(&prot));
        (text, prot)
    }

    fn codes(prot: &[u8]) -> Vec<u8> {
        prot.iter().map(|&b| residue_code(b).unwrap()).collect()
    }

    #[test]
    fn ungapped_extension_stops_at_the_x_drop() {
        let a = codes(b"WWWWWAAAAAPPPPPPPPPPPPPPPPWWWWWWWWWW");
        let b = codes(b"WWWWWAAAAAGGGGGGGGGGGGGGGGWWWWWWWWWW");
        // W-W scores 11, A-A 4, P-G -2: the run from the start reaches 5 * 11 + 5 * 4 = 75 and
        // falls more than XDROP below it inside the P/G block, so the Ws after it are not reached
        assert_eq!(ungapped_score(&a, &b, 0, 0), 75);
        // seeded inside the last W block: it extends right to the end and left into the P/G block
        // only until the X-drop
        assert_eq!(ungapped_score(&a, &b, 26, 26), 10 * 11);
        assert!((bits_needed(300, 100_000, 1e-6) - (3e7f64 / 1e-6).log2()).abs() < 1e-9);
    }

    #[test]
    fn a_diverged_homolog_keeps_its_gene_and_a_composition_match_does_not() {
        let (text, prot) = membrane_panel();
        let p = PanelParams::default();
        let kmers_only = Panel::parse(&text, 6, &[]).unwrap();
        let checked = Panel::parse(&text, 6, &[])
            .unwrap()
            .with_max_evalue(Some(LINEAGE_MAX_EVALUE));
        let mut hits = vec![0u32; 1];
        // the same gene in a distant species: every third residue changed, the conserved parts
        // still sharing 6-mers here and there
        let mut homolog = prot.clone();
        for (i, r) in homolog.iter_mut().enumerate() {
            if i % 3 == 1 && !(120..160).contains(&i) {
                *r = if *r == b'K' { b'R' } else { b'K' };
            }
            if i % 37 == 0 {
                *r = b'W';
            }
        }
        let aa = codes(&homolog);
        assert_eq!(best_gene(&aa, &kmers_only, &p, &mut hits), Some(0));
        assert_eq!(best_gene(&aa, &checked, &p, &mut hits), Some(0));
        assert!(significant_bits(&aa, 0, &checked) > 100.0);
        // an unrelated protein with a transmembrane-like run of its own, of the same residues: it
        // shares six 6-mers with the gene by composition alone, and scores about as much against
        // the gene reversed as forward
        let mut chance = protein(150, 77);
        chance.extend_from_slice(b"IILFLLLIIFFILILLFIFILLFLLFLILLFFILLFFIIILILIVIIL");
        chance.extend(protein(100, 78));
        let aa = codes(&chance);
        assert_eq!(best_gene(&aa, &kmers_only, &p, &mut hits), Some(0));
        assert_eq!(
            best_gene(&aa, &checked, &p, &mut hits),
            None,
            "forward {:.1} bits, significant {:.1}, needed {:.1}",
            stretch_bits(&aa, 0, &checked),
            significant_bits(&aa, 0, &checked),
            bits_needed(aa.len(), checked.residues, LINEAGE_MAX_EVALUE)
        );
        assert!(hits.iter().all(|&h| h == 0), "scratch counts are reset");
    }

    #[test]
    fn the_animal_and_fungal_panels_check_alignments_the_plant_panel_does_not() {
        for l in ["metazoa", "fungi"] {
            let panel = Panel::for_lineage(l, 6, &[]).unwrap();
            assert_eq!(panel.max_evalue(), Some(LINEAGE_MAX_EVALUE), "{l}");
            assert!(panel.residues > 50_000);
        }
        assert_eq!(Panel::lineage_max_evalue("landplants"), None);
        assert_eq!(Panel::embedded(6, &[]).unwrap().max_evalue(), None);
    }

    #[test]
    fn the_lineage_panels_parse_with_their_codes() {
        let metazoa = Panel::for_lineage("metazoa", 6, &[]).unwrap();
        let fungi = Panel::for_lineage("fungi", 6, &[]).unwrap();
        for (panel, table) in [(&metazoa, 5), (&fungi, 4)] {
            assert!(panel.proteins > 500);
            assert!(panel.genes.iter().all(|g| g.class == Class::Mito));
            assert!(panel.genes.iter().any(|g| g.name == "cox1"));
            assert_eq!(panel.code, genetic_code(table).unwrap());
        }
        assert!(Panel::for_lineage("protists", 6, &[]).is_err());
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
        let dir = std::env::temp_dir().join(format!("ovasm-panel-call-{}", std::process::id()));
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
