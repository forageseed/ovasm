//! Stage 2b: which connected components of an assembly graph are the target organelle
//! (`ovasm components`).
//!
//! Reads pooled from depth clusters (`ovasm discover`) can hold nuclear high-copy sequence
//! that shares the organelle's depth cluster, and the other organelle's reads; both assemble
//! into components of their own beside the organelle's. Each component gets its total length,
//! its length-weighted mean depth and the panel genes its segments encode (the `panel`
//! machinery on assembled sequence: one stop-free stretch is enough for a gene here, a read
//! set needs two). The component with the most target genes is the reference. A component is
//! kept when it is the reference, when it holds target genes at a depth that is not far below
//! the reference's, or when it holds no genes and its depth is close to the reference's (the
//! small circles of a multi-chromosome mitochondrion may hold no panel gene at all); it is
//! removed when its genes are the other organelle's or when its depth is elsewhere. Without
//! a reference nothing is removed. Every component, kept or removed, is in the JSON report
//! with its length, depth, genes and the reason.
//!
//! The thresholds were set on cases whose published genomes had been seen (see
//! `ComponentParams`); they are not validated on unseen data.

use std::fs;
use std::path::Path;

use anyhow::{bail, Context, Result};
use rustc_hash::FxHashSet;
use serde::Serialize;

use crate::evidence::parse_gfa;
use crate::panel::{sequence_genes, Panel, PanelParams};

/// The factors come from graphs whose published genomes had been seen (v3 holdout and two
/// later Oatk species, reads pooled from seeds and depth clusters), not from unseen data.
#[derive(Debug, Clone, Serialize)]
pub struct ComponentParams {
    /// A component without genes is kept when its depth is within this factor of the
    /// reference's, either way. Components without genes that are in no published genome sat
    /// at 0.09-0.45 of the reference depth (0.45: a 14.5 kb circle beside the Climacium
    /// mitochondrion) and at 2.5-3.3 times it (the nuclear rDNA unit, a 12.4 kb circle, beside
    /// the Ulota plastome). No published component was without genes in these graphs; those
    /// with genes sat at 0.51-1.36 (0.51: a 19 kb circle of the ten Galeopsis has).
    pub depth_factor: f64,
    /// A component with target genes is kept unless the reference is more than this many times
    /// deeper: organelle-derived sequence in the nuclear genome holds organelle genes too, one
    /// gene at 0.10-0.28 of the reference depth in the Climacium and Luzula mitochondrial
    /// graphs, against 0.51 and more for the published components.
    pub gene_depth_factor: f64,
    /// Genes of the other organelle that make a component the other organelle's, when they
    /// are also at least twice its target genes (the margin of `panel`'s calls). Fewer are
    /// what organelle-derived inserts leave (plastid genes in a mitochondrial chromosome).
    pub min_other_genes: usize,
    /// Amino-acid k-mer length, shortest stop-free stretch and shared k-mers per stretch of
    /// the gene search (as in `ovasm panel`).
    pub k: usize,
    pub min_orf: usize,
    pub min_hits: u32,
}

impl Default for ComponentParams {
    fn default() -> Self {
        Self {
            depth_factor: 2.0,
            gene_depth_factor: 2.5,
            min_other_genes: 5,
            k: 6,
            min_orf: 80,
            min_hits: 5,
        }
    }
}

/// What the rule looks at in one component.
#[derive(Debug, Clone, Copy)]
pub struct Stat {
    pub length: usize,
    pub depth: Option<f64>,
    /// Distinct panel genes of the target organelle and of the other one.
    pub target_genes: usize,
    pub other_genes: usize,
}

/// The reference component and, for every component, whether it stays and why.
pub fn decide(stats: &[Stat], p: &ComponentParams) -> (Option<usize>, Vec<(bool, &'static str)>) {
    // the reference: most target genes, more of them than of the other organelle's, with a
    // depth; the deeper one wins a tie (nuclear repeats are long, not deep), then the earlier
    let reference = stats
        .iter()
        .enumerate()
        .filter(|(_, s)| {
            s.target_genes > 0 && s.target_genes > s.other_genes && s.depth.is_some_and(|d| d > 0.0)
        })
        .max_by(|a, b| {
            (a.1.target_genes.cmp(&b.1.target_genes))
                .then(
                    a.1.depth
                        .unwrap_or(0.0)
                        .total_cmp(&b.1.depth.unwrap_or(0.0)),
                )
                .then(b.0.cmp(&a.0))
        })
        .map(|(i, _)| i);
    let Some(r) = reference else {
        return (None, vec![(true, "no_reference"); stats.len()]);
    };
    let ref_depth = stats[r].depth.unwrap_or(0.0);
    let calls = stats
        .iter()
        .enumerate()
        .map(|(i, s)| {
            if i == r {
                return (true, "reference");
            }
            let Some(d) = s.depth else {
                return (true, "no_depth");
            };
            if s.other_genes >= p.min_other_genes && s.other_genes >= 2 * s.target_genes {
                (false, "other_organelle")
            } else if s.target_genes > 0 {
                if d * p.gene_depth_factor >= ref_depth {
                    (true, "target_genes")
                } else {
                    (false, "target_genes_but_depth_far_below_reference")
                }
            } else if d * p.depth_factor < ref_depth {
                (false, "depth_below_reference")
            } else if d > ref_depth * p.depth_factor {
                (false, "depth_above_reference")
            } else {
                (true, "depth_like_reference")
            }
        })
        .collect();
    (reference, calls)
}

#[derive(Serialize)]
pub struct ComponentReport {
    pub id: usize,
    pub segments: Vec<String>,
    /// Sum of the segments' lengths.
    pub length: usize,
    /// Length-weighted mean of the segments' depths (None when a segment has no depth tag).
    pub depth: Option<f64>,
    /// `depth` over the reference component's.
    pub depth_ratio: Option<f64>,
    /// Panel genes of the target organelle on the segments: gene and stop-free stretches.
    pub target_genes: Vec<(String, u64)>,
    /// The same for the other organelle.
    pub other_genes: Vec<(String, u64)>,
    pub kept: bool,
    /// "reference", "target_genes", "depth_like_reference", "no_reference", "no_depth" (kept);
    /// "other_organelle", "target_genes_but_depth_far_below_reference", "depth_below_reference",
    /// "depth_above_reference" (removed).
    pub reason: &'static str,
}

#[derive(Serialize)]
pub struct ComponentsReport {
    pub tool: &'static str,
    pub version: &'static str,
    pub graph: String,
    pub organelle: String,
    pub params: ComponentParams,
    pub panel_proteins: usize,
    pub panel_genes: usize,
    /// `id` of the reference component (None: no component holds target genes, nothing removed).
    pub reference_component: Option<usize>,
    pub reference_depth: Option<f64>,
    pub kept_components: usize,
    pub kept_segments: usize,
    pub kept_length: usize,
    pub removed_components: usize,
    pub removed_segments: usize,
    pub removed_length: usize,
    /// GFA lines other than S left out because they name a removed segment.
    pub removed_lines: usize,
    pub components: Vec<ComponentReport>,
}

/// Connected components of the graph: for each, its segments in graph order; components in
/// the order of their first segment.
fn components_of(segments: usize, links: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let mut parent: Vec<usize> = (0..segments).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    for &(a, b) in links {
        let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
        // the smaller index is the root, so a component is found at its first segment
        parent[ra.max(rb)] = ra.min(rb);
    }
    let mut slot = vec![usize::MAX; segments];
    let mut out: Vec<Vec<usize>> = Vec::new();
    for s in 0..segments {
        let r = find(&mut parent, s);
        if slot[r] == usize::MAX {
            slot[r] = out.len();
            out.push(Vec::new());
        }
        out[slot[r]].push(s);
    }
    out
}

/// The GFA without the removed segments and every line that names one (links, containments,
/// jumps, paths); the other lines as they were. Returns the text and the non-S lines left out.
fn filter_gfa(text: &str, removed: &FxHashSet<&str>) -> (String, usize) {
    if removed.is_empty() {
        return (text.to_string(), 0);
    }
    let mut out = String::with_capacity(text.len());
    let mut dropped = 0;
    for line in text.split_inclusive('\n') {
        let f: Vec<&str> = line.trim_end_matches(['\n', '\r']).split('\t').collect();
        let gone = |name: &str| removed.contains(name);
        let drop = match f.first().copied() {
            Some("S") => {
                if f.get(1).is_some_and(|n| gone(n)) {
                    continue;
                }
                false
            }
            Some("L" | "C" | "J") => {
                f.get(1).is_some_and(|n| gone(n)) || f.get(3).is_some_and(|n| gone(n))
            }
            Some("P") => f.get(2).is_some_and(|steps| {
                steps
                    .split([',', ';'])
                    .any(|s| gone(s.trim_end_matches(['+', '-'])))
            }),
            _ => false,
        };
        if drop {
            dropped += 1;
        } else {
            out.push_str(line);
        }
    }
    (out, dropped)
}

pub fn run(
    gfa: &Path,
    organelle: &str,
    panel: &Panel,
    p: &ComponentParams,
    out_gfa: &Path,
    out_json: &Path,
) -> Result<ComponentsReport> {
    let target_is_mito = match organelle {
        "mitochondrion" => true,
        "plastid" => false,
        other => bail!("organelle {other:?} is neither mitochondrion nor plastid"),
    };
    if p.depth_factor.is_nan()
        || p.gene_depth_factor.is_nan()
        || p.depth_factor < 1.0
        || p.gene_depth_factor < 1.0
    {
        bail!("the depth factors must be at least 1");
    }
    let g = parse_gfa(gfa)?;
    let text = fs::read_to_string(gfa).with_context(|| format!("cannot read {}", gfa.display()))?;
    let tags = crate::unify::segment_tags(&text);
    let depth_of = |s: usize| {
        tags.get(&g.names[s])
            .and_then(|t| crate::unify::depth(t, g.seqs[s].len()))
            .map(|d| d.0)
    };
    let links: Vec<(usize, usize)> = g.links.iter().map(|&(a, b)| (a.0, b.0)).collect();
    let parts = components_of(g.names.len(), &links);
    // one stretch is a gene on assembled sequence, which holds each gene once
    let ip = PanelParams {
        k: p.k,
        min_orf: p.min_orf,
        min_hits: p.min_hits,
        min_gene_segments: 1,
        ..PanelParams::default()
    };
    let mut stats = Vec::new();
    let mut genes = Vec::new();
    for part in &parts {
        let length: usize = part.iter().map(|&s| g.seqs[s].len()).sum();
        let depth = part
            .iter()
            .map(|&s| depth_of(s).map(|d| d * g.seqs[s].len() as f64))
            .sum::<Option<f64>>()
            .map(|sum| sum / length.max(1) as f64);
        let (mito, plastid) =
            sequence_genes(part.iter().map(|&s| g.seqs[s].as_slice()), panel, &ip);
        let (target, other) = if target_is_mito {
            (mito, plastid)
        } else {
            (plastid, mito)
        };
        stats.push(Stat {
            length,
            depth,
            target_genes: target.genes,
            other_genes: other.genes,
        });
        genes.push((target.gene_segments, other.gene_segments));
    }
    let (reference, calls) = decide(&stats, p);
    let reference_depth = reference.and_then(|r| stats[r].depth);
    let mut removed: FxHashSet<&str> = FxHashSet::default();
    let mut components = Vec::new();
    for (id, (((part, stat), (target_genes, other_genes)), &(kept, reason))) in
        parts.iter().zip(&stats).zip(genes).zip(&calls).enumerate()
    {
        if !kept {
            removed.extend(part.iter().map(|&s| g.names[s].as_str()));
        }
        components.push(ComponentReport {
            id,
            segments: part.iter().map(|&s| g.names[s].clone()).collect(),
            length: stat.length,
            depth: stat.depth,
            depth_ratio: match (stat.depth, reference_depth) {
                (Some(d), Some(r)) => Some(d / r),
                _ => None,
            },
            target_genes,
            other_genes,
            kept,
            reason,
        });
    }
    let (filtered, removed_lines) = filter_gfa(&text, &removed);
    fs::create_dir_all(out_gfa.parent().unwrap_or(Path::new(".")))?;
    fs::write(out_gfa, filtered)?;
    let sum = |kept: bool, f: &dyn Fn(&ComponentReport) -> usize| -> usize {
        components.iter().filter(|c| c.kept == kept).map(f).sum()
    };
    let (panel_proteins, panel_genes) = panel.sizes();
    let report = ComponentsReport {
        tool: "ovasm components",
        version: env!("CARGO_PKG_VERSION"),
        graph: gfa.display().to_string(),
        organelle: organelle.to_string(),
        params: p.clone(),
        panel_proteins,
        panel_genes,
        reference_component: reference,
        reference_depth,
        kept_components: sum(true, &|_| 1),
        kept_segments: sum(true, &|c| c.segments.len()),
        kept_length: sum(true, &|c| c.length),
        removed_components: sum(false, &|_| 1),
        removed_segments: sum(false, &|c| c.segments.len()),
        removed_length: sum(false, &|c| c.length),
        removed_lines,
        components,
    };
    fs::write(out_json, serde_json::to_string_pretty(&report)?)?;
    Ok(report)
}

/// The outcome and one line per removed component for the terminal.
pub fn summary(r: &ComponentsReport) -> String {
    use std::fmt::Write as _;
    let mut s = String::new();
    let _ = writeln!(
        s,
        "[ovasm] {} components: kept {} ({} bp), removed {} ({} bp); reference depth {}",
        r.components.len(),
        r.kept_components,
        r.kept_length,
        r.removed_components,
        r.removed_length,
        r.reference_depth
            .map_or("none (nothing removed)".to_string(), |d| format!("{d:.1}")),
    );
    for c in r.components.iter().filter(|c| !c.kept) {
        let _ = writeln!(
            s,
            "[ovasm] removed component {}: {} segments, {} bp, depth {}, {} target / {} other genes: {}",
            c.id,
            c.segments.len(),
            c.length,
            c.depth.map_or("?".to_string(), |d| format!("{d:.1}")),
            c.target_genes.len(),
            c.other_genes.len(),
            c.reason
        );
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    const RESIDUES: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
    /// One codon per residue, in the order of `RESIDUES`.
    const CODONS: [&[u8; 3]; 20] = [
        b"GCT", b"TGT", b"GAT", b"GAA", b"TTT", b"GGT", b"CAT", b"ATT", b"AAA", b"CTT", b"ATG",
        b"AAT", b"CCT", b"CAA", b"CGT", b"TCT", b"ACT", b"GTT", b"TGG", b"TAT",
    ];

    fn xorshift(s: &mut u64) -> u64 {
        *s ^= *s << 13;
        *s ^= *s >> 7;
        *s ^= *s << 17;
        *s
    }

    /// A deterministic protein of `n` residues.
    fn protein(n: usize, seed: u64) -> Vec<u8> {
        let mut s = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        (0..n)
            .map(|_| RESIDUES[(xorshift(&mut s) % 20) as usize])
            .collect()
    }

    /// DNA that translates to protein `seed` of the panel, then a stop.
    fn gene(seed: u64) -> String {
        let mut dna = Vec::new();
        for a in protein(200, seed) {
            dna.extend(CODONS[RESIDUES.iter().position(|&r| r == a).unwrap()]);
        }
        dna.extend(b"TAA");
        String::from_utf8(dna).unwrap()
    }

    /// Sequence without genes: a 22 bp period holding TAA on both strands, so every frame
    /// meets a stop codon within 66 bp.
    fn noncoding(n: usize, seed: u64) -> String {
        let mut s = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        let mut out = String::new();
        while out.len() < n {
            out.push_str("TAAGTGAGTTAACTCACTTA");
            for _ in 0..2 {
                out.push(b"ACGT"[(xorshift(&mut s) % 4) as usize] as char);
            }
        }
        out.truncate(n);
        out
    }

    fn panel() -> Panel {
        let mut t = String::new();
        for (class, name, seed) in [
            ("mito", "cox1", 1),
            ("mito", "atp1", 2),
            ("mito", "nad1", 3),
            ("plastid", "rbcl", 4),
            ("plastid", "psba", 5),
        ] {
            t += &format!(
                ">{class}|{name}|Sp_one\n{}\n",
                String::from_utf8(protein(200, seed)).unwrap()
            );
        }
        Panel::parse(&t, 6, &[]).unwrap()
    }

    fn stat(length: usize, depth: f64, target_genes: usize, other_genes: usize) -> Stat {
        Stat {
            length,
            depth: Some(depth),
            target_genes,
            other_genes,
        }
    }

    fn reasons(stats: &[Stat]) -> (Option<usize>, Vec<&'static str>) {
        let (r, calls) = decide(stats, &ComponentParams::default());
        for (kept, reason) in &calls {
            let keeps = matches!(
                *reason,
                "reference" | "target_genes" | "depth_like_reference" | "no_reference" | "no_depth"
            );
            assert_eq!(*kept, keeps, "{reason}");
        }
        (r, calls.into_iter().map(|c| c.1).collect())
    }

    #[test]
    fn components_follow_the_links_in_segment_order() {
        let parts = components_of(6, &[(4, 1), (1, 2), (5, 5)]);
        assert_eq!(parts, vec![vec![0], vec![1, 2, 4], vec![3], vec![5]]);
        assert!(components_of(0, &[]).is_empty());
    }

    #[test]
    fn nuclear_repeats_at_another_depth_go_and_the_gene_rich_circle_stays() {
        // Climacium dendroides mitochondrion (v3 holdout, union arm): the genome and four of
        // its 27 other components
        let (r, why) = reasons(&[
            stat(119_515, 10.6, 1, 0),
            stat(719_779, 15.7, 1, 4),
            stat(104_863, 80.5, 38, 0),
            stat(14_521, 35.9, 0, 0),
            stat(14_696, 9.6, 0, 1),
        ]);
        assert_eq!(r, Some(2));
        assert_eq!(
            why,
            [
                "target_genes_but_depth_far_below_reference",
                "target_genes_but_depth_far_below_reference",
                "reference",
                "depth_below_reference",
                "depth_below_reference",
            ]
        );
        // Luzula sylvatica mitochondrion: 12 kb with one mitochondrial gene, not published
        let (r, why) = reasons(&[stat(642_734, 26.0, 37, 7), stat(12_309, 7.4, 1, 0)]);
        assert_eq!(r, Some(0));
        assert_eq!(why[1], "target_genes_but_depth_far_below_reference");
    }

    #[test]
    fn the_circles_of_a_multi_chromosome_mitochondrion_stay() {
        // Galeopsis tetrahit: gene-bearing circles down to half the reference's depth; the
        // last line is made up, a circle without panel genes at a comparable depth
        let (r, why) = reasons(&[
            stat(38_710, 27.5, 6, 0),
            stat(68_384, 28.8, 9, 0),
            stat(19_192, 14.7, 4, 1),
            stat(15_993, 18.8, 2, 1),
            stat(8_286, 11.1, 0, 0),
            stat(8_904, 20.0, 0, 0),
        ]);
        assert_eq!(r, Some(1));
        assert_eq!(
            why,
            [
                "target_genes",
                "reference",
                "target_genes",
                "target_genes",
                "depth_below_reference",
                "depth_like_reference",
            ]
        );
    }

    #[test]
    fn a_deeper_component_without_genes_and_the_other_organelle_go() {
        // Ulota crispa plastid: a 12 kb circle at 3.3 times the plastome's depth; the last
        // three lines are made up
        let (r, why) = reasons(&[
            stat(114_670, 93.4, 63, 0),
            stat(12_405, 304.9, 0, 0),
            stat(8_413, 36.2, 0, 0),
            // a mitochondrion beside the plastome, whatever its depth
            stat(100_000, 90.0, 2, 30),
            // an organelle-derived insert's few genes do not make a component the other's
            stat(30_000, 80.0, 1, 2),
            stat(30_000, 80.0, 0, 3),
        ]);
        assert_eq!(r, Some(0));
        assert_eq!(
            why,
            [
                "reference",
                "depth_above_reference",
                "depth_below_reference",
                "other_organelle",
                "target_genes",
                "depth_like_reference",
            ]
        );
    }

    #[test]
    fn the_factors_are_inclusive_and_ties_go_to_the_deeper_component() {
        let (r, why) = reasons(&[
            stat(500_000, 16.0, 3, 0),
            stat(50_000, 40.0, 3, 0),
            stat(1_000, 20.0, 0, 0),
            stat(1_000, 80.0, 0, 0),
            stat(1_000, 19.9, 0, 0),
            stat(1_000, 80.1, 0, 0),
            stat(1_000, 15.9, 1, 0),
            stat(1_000, 4_000.0, 1, 0),
        ]);
        assert_eq!(r, Some(1));
        assert_eq!(
            why,
            [
                "target_genes",
                "reference",
                "depth_like_reference",
                "depth_like_reference",
                "depth_below_reference",
                "depth_above_reference",
                "target_genes_but_depth_far_below_reference",
                "target_genes",
            ]
        );
    }

    #[test]
    fn nothing_is_removed_without_a_reference_or_without_depth() {
        // no target genes anywhere, or only beside more of the other organelle's
        let (r, why) = reasons(&[stat(9_000, 5.0, 0, 0), stat(9_000, 500.0, 1, 2)]);
        assert_eq!(r, None);
        assert_eq!(why, ["no_reference", "no_reference"]);
        let mut no_depth = stat(9_000, 0.0, 0, 0);
        no_depth.depth = None;
        let (r, why) = reasons(&[stat(9_000, 50.0, 4, 0), no_depth]);
        assert_eq!(r, Some(0));
        assert_eq!(why, ["reference", "no_depth"]);
    }

    #[test]
    fn filtering_keeps_the_other_lines_as_they_are() {
        let text = "H\tVN:Z:1.0\nS\ta\tACGT\tdp:f:9.0\nS\tb\tACGT\nS\tc\tGGCC\r\n\
                    L\ta\t+\tb\t-\t0M\nL\tc\t+\tc\t+\t0M\nP\tp1\ta+,b-\t*\nP\tp2\tc+\t*\n# note";
        let none = FxHashSet::default();
        assert_eq!(filter_gfa(text, &none), (text.to_string(), 0));
        let removed: FxHashSet<&str> = ["c"].into_iter().collect();
        let (out, dropped) = filter_gfa(text, &removed);
        assert_eq!(
            out,
            "H\tVN:Z:1.0\nS\ta\tACGT\tdp:f:9.0\nS\tb\tACGT\nL\ta\t+\tb\t-\t0M\nP\tp1\ta+,b-\t*\n# note"
        );
        assert_eq!(dropped, 2);
    }

    #[test]
    fn a_graph_is_filtered_by_genes_and_depth_and_everything_is_reported() {
        let dir = std::env::temp_dir().join(format!("ovasm-components-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        // "organelle": two linked segments holding two mitochondrial genes (one on the reverse
        // strand); "repeat": no genes, shallow; "circle": no genes, the organelle's depth;
        // "numt": a mitochondrial gene far below it; "plastid": plastid genes; "deep": no genes
        let organelle_a = noncoding(300, 1) + &gene(1) + &noncoding(200, 2);
        let organelle_b = String::from_utf8(crate::evidence::revcomp(gene(2).as_bytes())).unwrap()
            + &noncoding(900, 3);
        let segments = [
            ("organelle_a", organelle_a, 80.0),
            ("repeat", noncoding(2_000, 4), 12.0),
            ("organelle_b", organelle_b, 60.0),
            ("circle", noncoding(1_500, 5), 45.0),
            ("numt", noncoding(100, 6) + &gene(3), 9.0),
            ("plastid", gene(4) + &noncoding(50, 7) + &gene(5), 70.0),
            ("deep", noncoding(700, 8), 400.0),
        ];
        let mut gfa = String::from("H\tVN:Z:1.0\n");
        for (name, seq, depth) in &segments {
            gfa += &format!("S\t{name}\t{seq}\tLN:i:{}\tdp:f:{depth:.1}\n", seq.len());
        }
        gfa += "L\torganelle_a\t+\torganelle_b\t-\t0M\nL\tcircle\t+\tcircle\t+\t0M\n";
        gfa += "L\trepeat\t+\trepeat\t-\t0M\n";
        fs::write(dir.join("g.gfa"), &gfa).unwrap();
        let p = ComponentParams {
            min_other_genes: 2,
            ..ComponentParams::default()
        };
        let r = run(
            &dir.join("g.gfa"),
            "mitochondrion",
            &panel(),
            &p,
            &dir.join("out.gfa"),
            &dir.join("out.json"),
        )
        .unwrap();
        let by_first = |name: &str| r.components.iter().find(|c| c.segments[0] == name).unwrap();
        let organelle = by_first("organelle_a");
        assert_eq!(organelle.segments, ["organelle_a", "organelle_b"]);
        assert_eq!(r.reference_component, Some(organelle.id));
        let (la, lb) = (segments[0].1.len(), segments[2].1.len());
        let depth = (80.0 * la as f64 + 60.0 * lb as f64) / (la + lb) as f64;
        assert!((r.reference_depth.unwrap() - depth).abs() < 1e-9);
        assert_eq!(organelle.length, la + lb);
        let names = |v: &[(String, u64)]| v.iter().map(|x| x.0.clone()).collect::<Vec<_>>();
        assert_eq!(names(&organelle.target_genes), ["atp1", "cox1"]);
        assert!(organelle.other_genes.is_empty());
        let why: Vec<(&str, bool, &str)> = r
            .components
            .iter()
            .map(|c| (c.segments[0].as_str(), c.kept, c.reason))
            .collect();
        assert_eq!(
            why,
            [
                ("organelle_a", true, "reference"),
                ("repeat", false, "depth_below_reference"),
                ("circle", true, "depth_like_reference"),
                ("numt", false, "target_genes_but_depth_far_below_reference"),
                ("plastid", false, "other_organelle"),
                ("deep", false, "depth_above_reference"),
            ]
        );
        assert_eq!(names(&by_first("numt").target_genes), ["nad1"]);
        assert_eq!(names(&by_first("plastid").other_genes), ["psba", "rbcl"]);
        assert_eq!((r.kept_components, r.removed_components), (2, 4));
        assert_eq!((r.kept_segments, r.removed_segments), (3, 4));
        assert_eq!(r.kept_length, la + lb + 1_500);
        assert_eq!(
            r.kept_length + r.removed_length,
            segments.iter().map(|s| s.1.len()).sum::<usize>()
        );
        assert_eq!(r.removed_lines, 1);
        // the kept lines, unchanged and in order
        let expect: String = gfa
            .split_inclusive('\n')
            .filter(|l| l.starts_with('H') || l.contains("organelle_") || l.contains("\tcircle\t"))
            .collect();
        assert_eq!(fs::read_to_string(dir.join("out.gfa")).unwrap(), expect);
        let json: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(dir.join("out.json")).unwrap()).unwrap();
        assert_eq!(json["components"].as_array().unwrap().len(), 6);
        assert_eq!(json["components"][1]["reason"], "depth_below_reference");
        assert_eq!(json["components"][1]["length"], 2_000);
        assert_eq!(json["components"][1]["kept"], false);

        // the same graph for the plastid: the plastid component is the reference
        let r = run(
            &dir.join("g.gfa"),
            "plastid",
            &panel(),
            &p,
            &dir.join("out.gfa"),
            &dir.join("out.json"),
        )
        .unwrap();
        let reference = &r.components[r.reference_component.unwrap()];
        assert_eq!(reference.segments, ["plastid"]);
        let organelle = &r.components[0];
        assert_eq!(
            (organelle.kept, organelle.reason),
            (false, "other_organelle")
        );

        // a graph without target genes comes back byte for byte
        let plain = format!(
            "H\tVN:Z:1.0\nS\tx\t{}\tdp:f:5.0\nS\ty\t{}\tdp:f:900.0\n",
            noncoding(500, 9),
            noncoding(500, 10)
        );
        fs::write(dir.join("plain.gfa"), &plain).unwrap();
        let r = run(
            &dir.join("plain.gfa"),
            "mitochondrion",
            &panel(),
            &p,
            &dir.join("plain.out.gfa"),
            &dir.join("plain.json"),
        )
        .unwrap();
        assert_eq!((r.reference_component, r.removed_components), (None, 0));
        assert_eq!(
            fs::read_to_string(dir.join("plain.out.gfa")).unwrap(),
            plain
        );
        assert!(run(
            &dir.join("plain.gfa"),
            "nucleus",
            &panel(),
            &p,
            &dir.join("plain.out.gfa"),
            &dir.join("plain.json"),
        )
        .is_err());
        fs::remove_dir_all(&dir).ok();
    }
}
