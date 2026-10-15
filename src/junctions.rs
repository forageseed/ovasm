//! Restore observed junctions hidden inside a unitig when a repeat is shorter than graph k.
//!
//! A long exact k-mer is an error-free observation, not a count of all original reads
//! spanning a structural junction. A minority connection can therefore lose all long-word
//! nodes while many complete reads still agree on its two unique flanks. Enumerate maximal
//! sub-k inverted repeats from the product's own backbone, and compare normal/alternative
//! contexts with the existing evidence algorithm. The original automatic min_count is the
//! unchanged minimum number of distinct eligible original records. No external assembly,
//! reference sequence, fitted cutoff or change to graph cleaning is used.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use rustc_hash::FxHashMap;
use serde::Serialize;

use crate::evidence::{self, EvidenceParams, Graph, JunctionWitness, Oriented};
use crate::kmer::CanonicalKmers;

/// These are the existing HiFi `evidence` preset, not a separate rescue calibration.
pub(crate) fn hifi_evidence() -> EvidenceParams {
    EvidenceParams {
        short_read_adjacency: false,
        k: 21,
        flank: 2000,
        min_anchor: 500,
        diag_tol: 30,
        min_side_hits: 10,
        spacing_tol: 0.05,
        bootstrap: 0,
        seed: 0,
        max_repeat_len: 20_000,
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize)]
pub(crate) struct InvertedRepeat {
    pub segment: usize,
    pub left: usize,
    pub right: usize,
    pub length: usize,
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => 0,
    }
}

/// Maximal nonoverlapping inverted repeats with room for the existing unique-flank
/// anchors. A repeat already as long as k is not a hidden sub-k junction.
pub(crate) fn inverted_repeats(
    g: &Graph,
    k: usize,
    s: usize,
    anchor: usize,
) -> Vec<InvertedRepeat> {
    let mut found = BTreeSet::new();
    for (segment, seq) in g.seqs.iter().enumerate() {
        let mut words: FxHashMap<u64, Vec<(usize, bool)>> = FxHashMap::default();
        for (pos, word) in CanonicalKmers::new(seq, s) {
            let forward = seq[pos..pos + s] <= evidence::revcomp(&seq[pos..pos + s])[..];
            words.entry(word).or_default().push((pos, forward));
        }
        let mut covered: FxHashMap<usize, Vec<(usize, usize)>> = FxHashMap::default();
        for occurrences in words.values() {
            for (i, &(left, forward)) in occurrences.iter().enumerate() {
                for &(right, other_forward) in &occurrences[i + 1..] {
                    if forward == other_forward || right - left < s {
                        continue;
                    }
                    let diagonal = left + right;
                    if covered
                        .get(&diagonal)
                        .is_some_and(|v| v.iter().any(|&(a, b)| a <= left && left <= b))
                    {
                        continue;
                    }
                    let (mut a, mut b, mut length) = (left, right, s);
                    while a > 0
                        && b + length < seq.len()
                        && complement(seq[a - 1]) != 0
                        && seq[a - 1] == complement(seq[b + length])
                    {
                        a -= 1;
                        length += 1;
                    }
                    while a + length + 1 < b
                        && complement(seq[a + length]) != 0
                        && seq[a + length] == complement(seq[b - 1])
                    {
                        b -= 1;
                        length += 1;
                    }
                    covered
                        .entry(diagonal)
                        .or_default()
                        .push((a, a + length - s));
                    if length >= k
                        || a <= anchor
                        || b + length + anchor >= seq.len()
                        || b - (a + length) < 2 * anchor
                    {
                        continue;
                    }
                    found.insert(InvertedRepeat {
                        segment,
                        left: a,
                        right: b,
                        length,
                    });
                }
            }
        }
    }
    found.into_iter().collect()
}

fn hypothesis_graph(g: &Graph, repeats: &[InvertedRepeat], flank: usize) -> Graph {
    let mut out = Graph {
        names: Vec::new(),
        seqs: Vec::new(),
        links: Vec::new(),
        overlaps: Vec::new(),
    };
    for (index, r) in repeats.iter().enumerate() {
        let seq = &g.seqs[r.segment];
        let (a, b, n) = (r.left, r.right, r.length);
        let offset = out.seqs.len();
        out.names.extend([
            format!("j{index}_left"),
            format!("j{index}_middle"),
            format!("j{index}_right"),
        ]);
        out.seqs.extend([
            seq[a.saturating_sub(flank)..a + n].to_vec(),
            seq[a..b + n].to_vec(),
            seq[b..seq.len().min(b + n + flank)].to_vec(),
        ]);
        for (from, to) in [
            ((offset, true), (offset + 1, true)),
            ((offset + 1, true), (offset + 2, true)),
            ((offset, true), (offset + 1, false)),
            ((offset + 1, false), (offset + 2, true)),
        ] {
            out.links.push((from, to));
            out.overlaps.push(n);
        }
    }
    out
}

#[derive(Clone, Debug, Serialize)]
pub struct JunctionCall {
    pub segment: String,
    pub left_repeat: (usize, usize),
    pub right_repeat: (usize, usize),
    pub normal_reads: [u64; 2],
    pub alternative_reads: [u64; 2],
    pub nuclear_flanked_reads: [u64; 2],
    pub ambiguous_reads: [u64; 2],
    pub retained: [bool; 2],
    pub reason: [String; 2],
    pub witnesses: [Vec<JunctionWitness>; 2],
}

#[derive(Default, Serialize)]
pub struct JunctionReport {
    pub method: &'static str,
    pub source: &'static str,
    pub support_threshold: u32,
    pub support_threshold_changed: bool,
    pub records_scanned: u64,
    pub unique_records: u64,
    pub candidate_repeat_pairs: usize,
    pub retained_junctions: usize,
    pub backbone_segments: usize,
    pub backbone_links: usize,
    pub backbone_sequence_preserved: bool,
    pub physical_topology_identified: bool,
    pub calls: Vec<JunctionCall>,
    pub read_context: Option<crate::read_context::ContextSummary>,
    pub elapsed_seconds: f64,
}

#[derive(Clone)]
struct Bridge {
    repeat: usize,
    side: usize,
    segment: usize,
    from: (usize, bool),
    to: (usize, bool),
    seq: Vec<u8>,
    reads: u64,
    normal: u64,
}

/// Minimal complete unique-flank context, copied exclusively from the backbone.
fn bridge(
    g: &Graph,
    r: &InvertedRepeat,
    repeat: usize,
    side: usize,
    anchor: usize,
    reads: u64,
    normal: u64,
) -> Bridge {
    let seq = &g.seqs[r.segment];
    let (a, b, n) = (r.left, r.right, r.length);
    let (from, to, word) = if side == 0 {
        let mut word = seq[a - anchor..a + n].to_vec();
        word.extend(evidence::revcomp(&seq[b - anchor..b]));
        ((a - anchor, true), (b - anchor, false), word)
    } else {
        let mut word = evidence::revcomp(&seq[a + n..a + n + anchor]);
        word.extend_from_slice(&seq[b..b + n + anchor]);
        ((a + n + anchor, false), (b + n + anchor, true), word)
    };
    Bridge {
        repeat,
        side,
        segment: r.segment,
        from,
        to,
        seq: word,
        reads,
        normal,
    }
}

fn parent_tags(text: &str) -> BTreeMap<&str, Vec<&str>> {
    text.lines()
        .filter_map(|line| {
            let f: Vec<_> = line.split('\t').collect();
            (f.first() == Some(&"S")).then(|| {
                (
                    f[1],
                    f[3..]
                        .iter()
                        .copied()
                        .filter(|t| !t.starts_with("LN:") && !t.starts_with("NC:"))
                        .collect(),
                )
            })
        })
        .collect()
}

/// Split original nodes at bridge attachment coordinates, preserve every original base
/// and adjacency, and add only promoted context nodes. Never replace the backbone by a
/// reference haplotype or infer joint phasing from two separately supported junctions.
fn rewrite(g: &Graph, text: &str, bridges: &[Bridge]) -> Result<String> {
    let tags = parent_tags(text);
    let mut cuts: Vec<BTreeSet<usize>> = g
        .seqs
        .iter()
        .map(|s| BTreeSet::from([0, s.len()]))
        .collect();
    for b in bridges {
        cuts[b.segment].extend([b.from.0, b.to.0]);
    }
    let mut out = String::from("H\tVN:Z:1.0\n");
    let mut parts: Vec<Vec<(usize, usize, String)>> = Vec::new();
    for (i, cut) in cuts.iter().enumerate() {
        let positions: Vec<_> = cut.iter().copied().collect();
        let mut own = Vec::new();
        let mut restored = Vec::new();
        for (j, w) in positions.windows(2).enumerate() {
            let (a, b) = (w[0], w[1]);
            let seq = &g.seqs[i][a..b];
            if seq.is_empty() {
                bail!("empty backbone piece");
            }
            restored.extend_from_slice(seq);
            let name = if positions.len() == 2 {
                g.names[i].clone()
            } else {
                format!("{}__p{j}", g.names[i])
            };
            write!(
                out,
                "S\t{name}\t{}\tLN:i:{}",
                std::str::from_utf8(seq)?,
                seq.len()
            )?;
            for tag in &tags[g.names[i].as_str()] {
                write!(out, "\t{tag}")?;
            }
            if positions.len() > 2 {
                write!(out, "\tPS:Z:{}\tPI:Z:{a}-{b}", g.names[i])?;
            }
            out.push('\n');
            own.push((a, b, name));
        }
        if restored != g.seqs[i] {
            bail!("junction rewrite changed original backbone sequence");
        }
        parts.push(own);
    }
    let mut edges: BTreeMap<(String, bool, String, bool), usize> = BTreeMap::new();
    let mut insert = |a: String, af: bool, b: String, bf: bool, ov: usize| {
        let f = (a.clone(), af, b.clone(), bf);
        let r = (b, !bf, a, !af);
        edges.insert(f.min(r), ov);
    };
    let terminal = |o: Oriented, exit: bool| {
        let p = &parts[o.0];
        if o.1 == exit {
            p.last().unwrap().2.clone()
        } else {
            p.first().unwrap().2.clone()
        }
    };
    for (i, &(a, b)) in g.links.iter().enumerate() {
        insert(
            terminal(a, true),
            a.1,
            terminal(b, false),
            b.1,
            g.overlaps[i],
        );
    }
    for own in &parts {
        for w in own.windows(2) {
            insert(w[0].2.clone(), true, w[1].2.clone(), true, 0);
        }
    }
    for b in bridges {
        let name = format!("jr{}_{}", b.repeat, b.side);
        if g.names.contains(&name) {
            bail!("structural junction name collides with an input segment");
        }
        let depth = tags[g.names[b.segment].as_str()]
            .iter()
            .find_map(|t| t.strip_prefix("dp:f:")?.parse::<f64>().ok());
        write!(
            out,
            "S\t{name}\t{}\tLN:i:{}\tJS:i:{}\tNS:i:{}\tJR:Z:read_backed_alternative",
            std::str::from_utf8(&b.seq)?,
            b.seq.len(),
            b.reads,
            b.normal
        )?;
        if let Some(d) = depth {
            // Project an observed allele fraction into the existing node-depth units.
            // Counting raw reads as error-free long-k-mer depth would inflate minor arms.
            let allele = d * b.reads as f64 / (b.normal + b.reads).max(1) as f64;
            write!(out, "\tdp:f:{allele:.6}\tdm:f:{allele:.6}")?;
        }
        out.push('\n');
        let source = parts[b.segment]
            .iter()
            .find(|&&(a, z, _)| {
                if b.from.1 {
                    z == b.from.0
                } else {
                    a == b.from.0
                }
            })
            .context("missing source attachment")?;
        let dest = parts[b.segment]
            .iter()
            .find(|&&(a, z, _)| if b.to.1 { a == b.to.0 } else { z == b.to.0 })
            .context("missing destination attachment")?;
        insert(source.2.clone(), b.from.1, name.clone(), true, 0);
        insert(name, true, dest.2.clone(), b.to.1, 0);
    }
    for ((a, af, b, bf), ov) in edges {
        writeln!(
            out,
            "L\t{a}\t{}\t{b}\t{}\t{ov}M",
            if af { '+' } else { '-' },
            if bf { '+' } else { '-' }
        )?;
    }
    Ok(out)
}

pub fn restore<F>(
    gfa: &Path,
    inputs: &[PathBuf],
    k: usize,
    s: usize,
    threshold: u32,
    nuclear_flanked: F,
) -> Result<JunctionReport>
where
    F: Fn(&crate::recruit::Rec, i64, i64) -> Result<bool> + Sync,
{
    let start = Instant::now();
    let g = evidence::parse_gfa(gfa)?;
    let text = fs::read_to_string(gfa)?;
    let p = hifi_evidence();
    let repeats = inverted_repeats(&g, k, s, p.min_anchor);
    let mut report = JunctionReport {
        method: "sub-k inverted-repeat complete-read junction support",
        source: "current backbone and original input records; no external assembly",
        support_threshold: threshold,
        candidate_repeat_pairs: repeats.len(),
        backbone_segments: g.seqs.len(),
        backbone_links: g.links.len(),
        backbone_sequence_preserved: true,
        ..Default::default()
    };
    if repeats.is_empty() {
        report.elapsed_seconds = start.elapsed().as_secs_f64();
        return Ok(report);
    }
    let hypotheses = hypothesis_graph(&g, &repeats, p.flank);
    let counts = evidence::junction_support(&hypotheses, inputs, &p, nuclear_flanked)?;
    report.records_scanned = counts.records_scanned;
    report.unique_records = counts.unique_records;
    let eligible = |i: usize| {
        counts.witnesses[i]
            .iter()
            .filter(|r| !r.nuclear_flanked)
            .count() as u64
    };
    let mut bridges = Vec::new();
    for (index, r) in repeats.iter().enumerate() {
        let normal = [eligible(index * 4), eligible(index * 4 + 1)];
        let alternative = [eligible(index * 4 + 2), eligible(index * 4 + 3)];
        let nuclear = [
            counts.witnesses[index * 4 + 2].len() as u64 - alternative[0],
            counts.witnesses[index * 4 + 3].len() as u64 - alternative[1],
        ];
        let mut call = JunctionCall {
            segment: g.names[r.segment].clone(),
            left_repeat: (r.left, r.left + r.length),
            right_repeat: (r.right, r.right + r.length),
            normal_reads: normal,
            alternative_reads: alternative,
            nuclear_flanked_reads: nuclear,
            ambiguous_reads: [
                counts.ambiguous_reads[index * 4 + 2],
                counts.ambiguous_reads[index * 4 + 3],
            ],
            retained: [false; 2],
            reason: [String::new(), String::new()],
            witnesses: [
                counts.witnesses[index * 4 + 2].clone(),
                counts.witnesses[index * 4 + 3].clone(),
            ],
        };
        for side in 0..2 {
            if alternative[side] < threshold as u64 {
                call.reason[side] = "fewer_eligible_original_reads_than_unchanged_min_count".into();
                continue;
            }
            let b = bridge(
                &g,
                r,
                index,
                side,
                p.min_anchor,
                alternative[side],
                normal[side],
            );
            // An end-overlap cannot be longer than a newly split end piece. Do not
            // silently corrupt such a graph or manufacture upstream anchor sequence.
            let mut left = 0;
            let mut right = 0;
            for (li, &(a, z)) in g.links.iter().enumerate() {
                for (o, exit) in [(a, true), (z, false)] {
                    if o.0 == r.segment {
                        if o.1 == exit {
                            right = right.max(g.overlaps[li]);
                        } else {
                            left = left.max(g.overlaps[li]);
                        }
                    }
                }
            }
            let lo = b.from.0.min(b.to.0);
            let hi = b.from.0.max(b.to.0);
            if lo <= left || g.seqs[r.segment].len() - hi <= right {
                call.reason[side] = "unresolved_existing_graph_end_overlap".into();
                continue;
            }
            call.retained[side] = true;
            call.reason[side] = "independently_supported_complete_read_context".into();
            bridges.push(b);
            report.retained_junctions += 1;
        }
        report.calls.push(call);
    }
    if !bridges.is_empty() {
        let restored = rewrite(&g, &text, &bridges)?;
        let backup = gfa.with_extension("backbone.gfa");
        fs::write(&backup, &text)?;
        let pending = gfa.with_extension("junction-pending.gfa");
        fs::write(&pending, restored)?;
        // Parse and validate every original/added GFA overlap before publishing a report.
        let rewritten = evidence::parse_gfa(&pending)?;
        for (i, &(a, b)) in rewritten.links.iter().enumerate() {
            let ov = rewritten.overlaps[i];
            let x = rewritten.oriented_seq(a);
            let y = rewritten.oriented_seq(b);
            if x[x.len() - ov..] != y[..ov] {
                bail!("junction restoration made an inconsistent GFA overlap");
            }
        }
        fs::rename(&pending, gfa)?;
    }
    report.elapsed_seconds = start.elapsed().as_secs_f64();
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn seq(n: usize, mut x: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x & 3) as usize]
            })
            .collect()
    }
    fn graph(s: Vec<u8>) -> Graph {
        Graph {
            names: vec!["B".into()],
            seqs: vec![s],
            links: vec![],
            overlaps: vec![],
        }
    }
    fn toy() -> (Graph, InvertedRepeat) {
        let mut s = seq(2600, 19);
        let repeat = seq(83, 71);
        let a = s.len();
        s.extend(&repeat);
        s.extend(seq(5000, 29));
        let b = s.len();
        s.extend(evidence::revcomp(&repeat));
        s.extend(seq(2600, 31));
        let g = graph(s);
        let r = inverted_repeats(&g, 1001, 31, 500)
            .into_iter()
            .find(|r| r.left <= a && r.right <= b && r.length >= 83)
            .unwrap();
        (g, r)
    }
    #[test]
    fn discovers_maximal_subk_repeat_from_backbone_only() {
        let (g, r) = toy();
        assert_eq!(
            &g.seqs[0][r.left..r.left + r.length],
            &evidence::revcomp(&g.seqs[0][r.right..r.right + r.length])
        );
        assert!(inverted_repeats(&g, r.length, 31, 500)
            .iter()
            .all(|x| x.length < r.length));
    }
    #[test]
    fn normal_and_alternative_context_overlaps_are_exact() {
        let (g, r) = toy();
        let h = hypothesis_graph(&g, &[r], 2000);
        for (i, &(a, b)) in h.links.iter().enumerate() {
            let x = h.oriented_seq(a);
            let y = h.oriented_seq(b);
            let n = h.overlaps[i];
            assert_eq!(x[x.len() - n..], y[..n]);
        }
    }
    #[test]
    fn rewrites_multiple_connections_without_changing_original_bases() {
        let (g, r) = toy();
        let bridges = vec![
            bridge(&g, &r, 0, 0, 500, 20, 200),
            bridge(&g, &r, 0, 1, 500, 30, 300),
        ];
        let text = format!(
            "H\tVN:Z:1.0\nS\tB\t{}\tdp:f:100.0\n",
            std::str::from_utf8(&g.seqs[0]).unwrap()
        );
        let out = rewrite(&g, &text, &bridges).unwrap();
        let backbone: Vec<u8> = out
            .lines()
            .filter_map(|l| {
                let f: Vec<_> = l.split('\t').collect();
                (f[0] == "S" && f[1].starts_with("B__p")).then(|| f[2].as_bytes().to_vec())
            })
            .flatten()
            .collect();
        assert_eq!(backbone, g.seqs[0]);
        assert!(out.contains("dp:f:9.090909"));
        assert_eq!(out.lines().filter(|l| l.starts_with("S\tjr")).count(), 2);
    }

    fn workspace(tag: &str) -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "ovasm-junction-{tag}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir_all(&path).unwrap();
        path
    }
    fn fixture(dir: &Path) -> (PathBuf, PathBuf) {
        let (g, r) = toy();
        let path = dir.join("assembly.gfa");
        fs::write(
            &path,
            format!(
                "H\tVN:Z:1.0\nS\tB\t{}\tdp:f:100.0\n",
                std::str::from_utf8(&g.seqs[0]).unwrap()
            ),
        )
        .unwrap();
        let h = hypothesis_graph(&g, &[r], 2000);
        let reads = dir.join("reads.fa");
        let mut text = String::new();
        for side in 0..2 {
            let i = 2 + side;
            let (a, b) = h.links[i];
            let mut seq = h.oriented_seq(a).to_vec();
            seq.extend_from_slice(&h.oriented_seq(b)[h.overlaps[i]..]);
            for n in 0..3 {
                writeln!(
                    text,
                    ">side{side}_{n}\n{}",
                    std::str::from_utf8(&seq).unwrap()
                )
                .unwrap();
            }
        }
        fs::write(&reads, text).unwrap();
        (path, reads)
    }
    #[test]
    fn distinct_complete_reads_promote_both_links_without_lowering_quorum() {
        let dir = workspace("support");
        let (gfa, reads) = fixture(&dir);
        let report = restore(&gfa, &[reads], 1001, 31, 3, |_, _, _| Ok(false)).unwrap();
        assert_eq!(report.retained_junctions, 2);
        assert!(report
            .calls
            .iter()
            .any(|r| r.alternative_reads == [3, 3] && r.retained == [true, true]));
        assert!(gfa.with_extension("backbone.gfa").is_file());
        fs::remove_dir_all(dir).unwrap();
    }
    #[test]
    fn duplicate_records_cannot_inflate_support() {
        let dir = workspace("duplicate");
        let (gfa, reads) = fixture(&dir);
        let report = restore(&gfa, &[reads.clone(), reads], 1001, 31, 4, |_, _, _| {
            Ok(false)
        })
        .unwrap();
        assert_eq!(report.records_scanned, 12);
        assert_eq!(report.unique_records, 6);
        assert_eq!(report.retained_junctions, 0);
        fs::remove_dir_all(dir).unwrap();
    }
    #[test]
    fn whole_read_nuclear_context_does_not_vote_as_organelle_junction_support() {
        let dir = workspace("nuclear");
        let (gfa, reads) = fixture(&dir);
        let before = fs::read(&gfa).unwrap();
        let report = restore(&gfa, &[reads], 1001, 31, 3, |_, _, _| Ok(true)).unwrap();
        assert_eq!(report.retained_junctions, 0);
        assert!(report
            .calls
            .iter()
            .any(|r| r.nuclear_flanked_reads == [3, 3]));
        assert_eq!(before, fs::read(&gfa).unwrap());
        fs::remove_dir_all(dir).unwrap();
    }
}
