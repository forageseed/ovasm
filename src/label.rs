//! Graph-first target assignment (`ovasm label`): which assembled sequence of a sample is its
//! mitochondrion, which its plastid.
//!
//! Seeds from other species decide which reads a target gets; when they are distant (a plastid
//! seed against a heavily rearranged plastome, mitochondrial seeds carrying plastid-derived
//! stretches) the wrong reads or none are recruited. This module asks the sample instead.
//! `discover` splits the organelle-like reads (deep along their whole length) into depth
//! clusters; each cluster holds genomes of one depth, so it is assembled on its own and the
//! assembler's depth-relative cleaning applies (in one graph, a 14x mitochondrion beside a
//! 200x plastid is cleaned away). Each unitig is then labelled by the organelle genes it
//! encodes: every reference protein's amino-acid k-mers are looked up in the unitig's six
//! frames, and a gene counts when enough of its k-mers occur. A unitig with at least
//! `min_genes` genes of one organelle and `majority` times as many as of the other is that
//! organelle's; plastid-derived stretches in a mitochondrion (MTPTs) leave a few plastid genes
//! beside many mitochondrial ones. Unitigs without genes (intergenic mitochondrial sequence)
//! take the label of their connected component when its labelled unitigs agree. The labelled
//! unitigs are the sample's own seeds for `recruit`.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;

use crate::assemble::{self, AssembleParams};
use crate::evidence::{parse_gfa, Oriented};

/// Amino acids encoded in a protein k-mer (5 bits each).
const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

pub const ORGANELLES: [&str; 2] = ["mitochondrion", "plastid"];

#[derive(Debug, Clone, Copy, Serialize)]
pub struct LabelParams {
    /// Amino-acid k-mer length (at most 12).
    pub k: usize,
    /// Distinct k-mers of one reference gene a unitig must contain for the gene to count.
    pub min_gene_kmers: usize,
    /// Genes of one organelle a unitig needs to be labelled by its own genes.
    pub min_genes: usize,
    /// ... and at least this many times as many as of the other organelle.
    pub majority: f64,
    /// Node k-mer length of the per-cluster assemblies.
    pub assemble_k: usize,
}

impl Default for LabelParams {
    fn default() -> Self {
        Self {
            k: 8,
            min_gene_kmers: 10,
            min_genes: 2,
            majority: 3.0,
            assemble_k: 1001,
        }
    }
}

fn aa_code(a: u8) -> Option<u64> {
    AA.iter().position(|&x| x == a).map(|i| i as u64)
}

/// Codon -> amino acid (standard code; plastid table 11 differs only in start codons).
fn translate_codon(c: &[u8]) -> u8 {
    const T: &[u8; 64] = b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
    let mut i = 0usize;
    for &b in c {
        i = i << 2
            | match b {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => return b'X',
            };
    }
    T[i]
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Packed k-mers of a protein sequence; k-mers across a stop or an unknown residue are skipped.
fn aa_kmers(prot: &[u8], k: usize) -> impl Iterator<Item = u64> + '_ {
    let mask = if 5 * k >= 64 {
        u64::MAX
    } else {
        (1u64 << (5 * k)) - 1
    };
    let mut acc = 0u64;
    let mut valid = 0usize;
    prot.iter().filter_map(move |&a| match aa_code(a) {
        Some(c) => {
            acc = (acc << 5 | c) & mask;
            valid += 1;
            (valid >= k).then_some(acc)
        }
        None => {
            valid = 0;
            None
        }
    })
}

/// Reference proteins indexed by amino-acid k-mer.
pub struct GeneIndex {
    k: usize,
    /// (organelle index, gene name) of each gene.
    pub genes: Vec<(usize, String)>,
    index: FxHashMap<u64, Vec<u32>>,
}

impl GeneIndex {
    /// Protein FASTA with headers `organelle|gene|source` (organelle: mitochondrion or plastid).
    /// Records of one organelle and gene name from several species form one gene.
    pub fn load(path: &Path, k: usize) -> Result<Self> {
        if !(1..=12).contains(&k) {
            bail!("amino-acid k must be between 1 and 12");
        }
        let mut reader = needletail::parse_fastx_file(path)
            .with_context(|| format!("cannot read gene database {}", path.display()))?;
        let mut ids: FxHashMap<(usize, String), u32> = FxHashMap::default();
        let mut genes = Vec::new();
        let mut index: FxHashMap<u64, Vec<u32>> = FxHashMap::default();
        while let Some(record) = reader.next() {
            let record =
                record.with_context(|| format!("malformed record in {}", path.display()))?;
            let header = String::from_utf8_lossy(record.id()).to_string();
            let mut parts = header.split('|');
            let (Some(org), Some(gene)) = (parts.next(), parts.next()) else {
                bail!("gene database header {header:?} is not organelle|gene|source");
            };
            let Some(o) = ORGANELLES.iter().position(|&x| x == org) else {
                bail!("gene database header {header:?}: unknown organelle {org:?}");
            };
            let key = (o, gene.to_ascii_lowercase());
            let id = *ids.entry(key.clone()).or_insert_with(|| {
                genes.push(key);
                genes.len() as u32 - 1
            });
            let seq = record.seq().to_ascii_uppercase();
            for km in aa_kmers(&seq, k) {
                let e = index.entry(km).or_default();
                if !e.contains(&id) {
                    e.push(id);
                }
            }
        }
        Ok(Self { k, genes, index })
    }

    /// Genes with at least `min_kmers` distinct k-mers in one of the six frames of `seq`.
    pub fn genes_in(&self, seq: &[u8], min_kmers: usize) -> Vec<u32> {
        let mut hits: FxHashMap<u32, FxHashSet<u64>> = FxHashMap::default();
        let rc = revcomp(seq);
        for strand in [seq, rc.as_slice()] {
            for frame in 0..3 {
                let prot: Vec<u8> = strand[frame.min(strand.len())..]
                    .chunks_exact(3)
                    .map(translate_codon)
                    .collect();
                for km in aa_kmers(&prot, self.k) {
                    if let Some(ids) = self.index.get(&km) {
                        for &g in ids {
                            hits.entry(g).or_default().insert(km);
                        }
                    }
                }
            }
        }
        let mut out: Vec<u32> = hits
            .into_iter()
            .filter(|(_, kms)| kms.len() >= min_kmers)
            .map(|(g, _)| g)
            .collect();
        out.sort_unstable();
        out
    }
}

/// Label from gene counts per organelle: `Some(Ok(o))` organelle o, `Some(Err(()))` genes of
/// both without a clear majority, `None` no genes.
fn label_of(counts: [usize; 2], p: &LabelParams) -> Option<std::result::Result<usize, ()>> {
    if counts == [0, 0] {
        return None;
    }
    let best = if counts[0] >= counts[1] { 0 } else { 1 };
    if counts[best] >= p.min_genes && counts[best] as f64 >= p.majority * counts[1 - best] as f64 {
        Some(Ok(best))
    } else {
        Some(Err(()))
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct UnitigLabel {
    pub name: String,
    pub length: usize,
    /// Genes found per organelle (mitochondrion, plastid).
    pub genes: [Vec<String>; 2],
    /// "mitochondrion", "plastid", "ambiguous" or "none" by its own genes.
    pub by_genes: String,
    /// Final label after component inheritance ("none" when not an organelle seed).
    pub label: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct ClusterLabels {
    pub reads: String,
    pub graph: String,
    pub unitigs: Vec<UnitigLabel>,
    /// Seed bases per organelle (mitochondrion, plastid).
    pub seed_bp: [usize; 2],
}

#[derive(Debug, Clone, Serialize)]
pub struct LabelReport {
    pub tool: &'static str,
    pub version: &'static str,
    pub params: LabelParams,
    pub gene_database: String,
    pub genes_indexed: [usize; 2],
    pub clusters: Vec<ClusterLabels>,
    /// Seed FASTA per organelle (empty string when no unitig carries that label).
    pub seeds: [String; 2],
    pub seed_bp: [usize; 2],
}

fn components(n: usize, links: &[(Oriented, Oriented)]) -> Vec<usize> {
    let mut parent: Vec<usize> = (0..n).collect();
    fn root(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    for &((a, _), (b, _)) in links {
        let (ra, rb) = (root(&mut parent, a), root(&mut parent, b));
        parent[ra] = rb;
    }
    (0..n).map(|x| root(&mut parent, x)).collect()
}

/// Label the unitigs of one assembly graph.
pub fn label_graph(gfa: &Path, genes: &GeneIndex, p: &LabelParams) -> Result<Vec<UnitigLabel>> {
    let g = parse_gfa(gfa)?;
    let found: Vec<Vec<u32>> = g
        .seqs
        .par_iter()
        .map(|s| genes.genes_in(s, p.min_gene_kmers))
        .collect();
    let counts: Vec<[usize; 2]> = found
        .iter()
        .map(|f| {
            let mut c = [0usize; 2];
            for &id in f {
                c[genes.genes[id as usize].0] += 1;
            }
            c
        })
        .collect();
    let own: Vec<Option<std::result::Result<usize, ()>>> =
        counts.iter().map(|&c| label_of(c, p)).collect();
    let comp = components(g.names.len(), &g.links);
    // labelled bases per organelle in each component, and in the whole graph (one depth cluster)
    let mut votes: FxHashMap<usize, [usize; 2]> = FxHashMap::default();
    let mut cluster = [0usize; 2];
    for (u, l) in own.iter().enumerate() {
        if let Some(Ok(o)) = l {
            votes.entry(comp[u]).or_default()[*o] += g.seqs[u].len();
            cluster[*o] += g.seqs[u].len();
        }
    }
    // The organelle every clearly labelled unitig of this depth cluster belongs to, if one.
    let cluster_organelle = match cluster {
        [m, 0] if m > 0 => Some(0),
        [0, c] if c > 0 => Some(1),
        _ => None,
    };
    let name = |l: &Option<std::result::Result<usize, ()>>| match l {
        Some(Ok(o)) => ORGANELLES[*o].to_string(),
        Some(Err(())) => "ambiguous".to_string(),
        None => "none".to_string(),
    };
    Ok((0..g.names.len())
        .map(|u| {
            let by_component = || {
                votes.get(&comp[u]).and_then(|v| match v {
                    [m, 0] if *m > 0 => Some(0),
                    [0, c] if *c > 0 => Some(1),
                    _ => None,
                })
            };
            let final_label = match own[u] {
                Some(Ok(o)) => Some(o),
                // Genes of both organelles without a 3:1 majority: a mitochondrial chromosome
                // with large plastid-derived stretches (Narthecium's 100.8 kb second chromosome,
                // 15 mitochondrial and 7 plastid genes, its own component) is still the
                // organelle its depth cluster holds, where its own genes agree.
                Some(Err(())) => by_component()
                    .or(cluster_organelle.filter(|&o| counts[u][o] > counts[u][1 - o])),
                None => by_component(),
            };
            let mut per = [Vec::new(), Vec::new()];
            for &id in &found[u] {
                let (o, ref gname) = genes.genes[id as usize];
                per[o].push(gname.clone());
            }
            UnitigLabel {
                name: g.names[u].clone(),
                length: g.seqs[u].len(),
                genes: per,
                by_genes: name(&own[u]),
                label: final_label.map_or("none".to_string(), |o| ORGANELLES[o].to_string()),
            }
        })
        .collect())
}

pub fn run(
    clusters: &[PathBuf],
    gene_db: &Path,
    p: LabelParams,
    out_dir: &Path,
) -> Result<LabelReport> {
    std::fs::create_dir_all(out_dir)?;
    let genes = GeneIndex::load(gene_db, p.k)?;
    let mut genes_indexed = [0usize; 2];
    for (o, _) in &genes.genes {
        genes_indexed[*o] += 1;
    }
    let mut seed_text = [String::new(), String::new()];
    let mut seed_bp = [0usize; 2];
    let mut out = Vec::new();
    for (i, reads) in clusters.iter().enumerate() {
        let gfa = out_dir.join(format!("cluster{}.gfa", i + 1));
        assemble::run(
            vec![reads.clone()],
            AssembleParams {
                k: p.assemble_k,
                s: 31,
                min_count: None,
            },
            &gfa,
            &out_dir.join(format!("cluster{}.assembly.json", i + 1)),
        )?;
        let unitigs = label_graph(&gfa, &genes, &p)?;
        let g = parse_gfa(&gfa)?;
        let mut cl_bp = [0usize; 2];
        for (u, l) in unitigs.iter().enumerate() {
            if let Some(o) = ORGANELLES.iter().position(|&x| x == l.label) {
                writeln!(
                    seed_text[o],
                    ">cluster{}_{}\n{}",
                    i + 1,
                    l.name,
                    std::str::from_utf8(&g.seqs[u])?
                )?;
                cl_bp[o] += l.length;
            }
        }
        for o in 0..2 {
            seed_bp[o] += cl_bp[o];
        }
        out.push(ClusterLabels {
            reads: reads.display().to_string(),
            graph: gfa.display().to_string(),
            unitigs,
            seed_bp: cl_bp,
        });
    }
    let mut seeds = [String::new(), String::new()];
    for o in 0..2 {
        if seed_bp[o] > 0 {
            let path = out_dir.join(format!("seed_{}.fasta", ORGANELLES[o]));
            std::fs::write(&path, &seed_text[o])?;
            seeds[o] = path.display().to_string();
        }
    }
    let report = LabelReport {
        tool: "ovasm label",
        version: env!("CARGO_PKG_VERSION"),
        params: p,
        gene_database: gene_db.display().to_string(),
        genes_indexed,
        clusters: out,
        seeds,
        seed_bp,
    };
    std::fs::write(
        out_dir.join("labels.json"),
        serde_json::to_string_pretty(&report)?,
    )?;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_codon_table_is_the_standard_code() {
        assert_eq!(translate_codon(b"ATG"), b'M');
        assert_eq!(translate_codon(b"TGG"), b'W');
        for stop in [b"TAA", b"TAG", b"TGA"] {
            assert_eq!(translate_codon(stop), b'*');
        }
        assert_eq!(translate_codon(b"GCT"), b'A');
        assert_eq!(translate_codon(b"TTT"), b'F');
        assert_eq!(translate_codon(b"AGA"), b'R');
        assert_eq!(translate_codon(b"ANA"), b'X');
    }

    /// DNA encoding `prot` with fixed codons.
    fn back_translate(prot: &[u8]) -> Vec<u8> {
        let table: FxHashMap<u8, &[u8]> = [
            (b'A', &b"GCT"[..]),
            (b'C', b"TGT"),
            (b'D', b"GAT"),
            (b'E', b"GAA"),
            (b'F', b"TTT"),
            (b'G', b"GGT"),
            (b'H', b"CAT"),
            (b'I', b"ATT"),
            (b'K', b"AAA"),
            (b'L', b"CTT"),
            (b'M', b"ATG"),
            (b'N', b"AAT"),
            (b'P', b"CCT"),
            (b'Q', b"CAA"),
            (b'R', b"CGT"),
            (b'S', b"TCT"),
            (b'T', b"ACT"),
            (b'V', b"GTT"),
            (b'W', b"TGG"),
            (b'Y', b"TAT"),
        ]
        .into_iter()
        .collect();
        prot.iter().flat_map(|a| table[a].iter().copied()).collect()
    }

    fn protein(n: usize, seed: u64) -> Vec<u8> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                AA[((s >> 33) % 20) as usize]
            })
            .collect()
    }

    #[test]
    fn genes_are_found_on_either_strand_and_any_frame() {
        let dir = std::env::temp_dir().join(format!("ovasm-label-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (cox1, psba) = (protein(300, 1), protein(300, 2));
        let db = dir.join("genes.faa");
        std::fs::write(
            &db,
            format!(
                ">mitochondrion|cox1|A\n{}\n>plastid|psbA|B\n{}\n",
                std::str::from_utf8(&cox1).unwrap(),
                std::str::from_utf8(&psba).unwrap()
            ),
        )
        .unwrap();
        let idx = GeneIndex::load(&db, 8).unwrap();
        let mut dna = b"AC".to_vec(); // frame 2
        dna.extend(back_translate(&cox1));
        let rc_psba = revcomp(&back_translate(&psba));
        dna.extend_from_slice(b"GGGGG");
        dna.extend(rc_psba);
        let found: Vec<&str> = idx
            .genes_in(&dna, 10)
            .iter()
            .map(|&g| idx.genes[g as usize].1.as_str())
            .collect();
        assert_eq!(found, vec!["cox1", "psba"]);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn labels_need_genes_and_a_clear_majority() {
        let p = LabelParams::default();
        assert_eq!(label_of([0, 0], &p), None);
        assert_eq!(
            label_of([5, 1], &p),
            Some(Ok(0)),
            "an MTPT leaves a plastid gene"
        );
        assert_eq!(
            label_of([1, 0], &p),
            Some(Err(())),
            "one gene is not enough"
        );
        assert_eq!(label_of([2, 2], &p), Some(Err(())));
        assert_eq!(label_of([0, 7], &p), Some(Ok(1)));
    }

    #[test]
    fn a_mixed_unitig_takes_its_clusters_organelle_where_its_genes_agree() {
        let dir = std::env::temp_dir().join(format!("ovasm-label-c-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mt: Vec<Vec<u8>> = (0..6).map(|i| protein(150, 30 + i)).collect();
        let pt: Vec<Vec<u8>> = (0..5).map(|i| protein(150, 50 + i)).collect();
        let mut db = String::new();
        for (i, g) in mt.iter().enumerate() {
            db.push_str(&format!(
                ">mitochondrion|m{i}|A\n{}\n",
                std::str::from_utf8(g).unwrap()
            ));
        }
        for (i, g) in pt.iter().enumerate() {
            db.push_str(&format!(
                ">plastid|p{i}|B\n{}\n",
                std::str::from_utf8(g).unwrap()
            ));
        }
        std::fs::write(dir.join("genes.faa"), db).unwrap();
        let dna = |genes: &[&Vec<u8>]| -> String {
            String::from_utf8(genes.iter().flat_map(|g| back_translate(g)).collect()).unwrap()
        };
        // u0: clearly mitochondrial; u1 (own component): 3 mitochondrial + 2 plastid genes;
        // u2 (own component): 2 mitochondrial + 3 plastid genes
        let gfa = dir.join("g.gfa");
        std::fs::write(
            &gfa,
            format!(
                "H\tVN:Z:1.0\nS\tu0\t{}\nS\tu1\t{}\nS\tu2\t{}\n",
                dna(&[&mt[0], &mt[1], &mt[2]]),
                dna(&[&mt[3], &mt[4], &mt[5], &pt[0], &pt[1]]),
                dna(&[&mt[0], &mt[1], &pt[2], &pt[3], &pt[4]]),
            ),
        )
        .unwrap();
        let idx = GeneIndex::load(&dir.join("genes.faa"), 8).unwrap();
        let l = label_graph(&gfa, &idx, &LabelParams::default()).unwrap();
        assert_eq!(l[0].label, "mitochondrion");
        assert_eq!(
            (l[1].by_genes.as_str(), l[1].label.as_str()),
            ("ambiguous", "mitochondrion")
        );
        assert_eq!(
            (l[2].by_genes.as_str(), l[2].label.as_str()),
            ("ambiguous", "none")
        );
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn unitigs_without_genes_take_their_components_label() {
        let dir = std::env::temp_dir().join(format!("ovasm-label-g-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let genes: Vec<Vec<u8>> = (0..3).map(|i| protein(200, 10 + i)).collect();
        let db = dir.join("genes.faa");
        let text: String = genes
            .iter()
            .enumerate()
            .map(|(i, g)| {
                format!(
                    ">mitochondrion|nad{i}|A\n{}\n",
                    std::str::from_utf8(g).unwrap()
                )
            })
            .collect();
        std::fs::write(&db, text).unwrap();
        let coding: Vec<u8> = genes.iter().flat_map(|g| back_translate(g)).collect();
        let s = |v: &[u8]| String::from_utf8(v.to_vec()).unwrap();
        let gfa = dir.join("g.gfa");
        std::fs::write(
            &gfa,
            format!(
                "H\tVN:Z:1.0\nS\tu0\t{}\nS\tu1\t{}\nS\tu2\t{}\nL\tu0\t+\tu1\t+\t0M\n",
                s(&coding),
                "ACGTTGCA".repeat(50),
                "TTGACCAG".repeat(50)
            ),
        )
        .unwrap();
        let idx = GeneIndex::load(&db, 8).unwrap();
        let l = label_graph(&gfa, &idx, &LabelParams::default()).unwrap();
        assert_eq!(l[0].by_genes, "mitochondrion");
        assert_eq!(l[0].label, "mitochondrion");
        assert_eq!(
            (l[1].by_genes.as_str(), l[1].label.as_str()),
            ("none", "mitochondrion")
        );
        assert_eq!(
            l[2].label, "none",
            "an unconnected piece without genes is not a seed"
        );
        std::fs::remove_dir_all(&dir).ok();
    }
}
