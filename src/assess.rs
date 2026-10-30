//! `ovasm check assess`: evidence about an assembly that its own run already holds, no other
//! tool. From the assembly graph (segment depths), the read evidence (link support), the
//! linearization (representative molecules), the chromosome calls and the built-in protein panel:
//!
//! - depth: the single-copy depth (median over the bases of segments the molecules walk once)
//!   and each segment's depth in those units;
//! - copy number: how often the molecules walk each segment against its depth, and the share
//!   of bases where the two agree (a molecule that walks a repeat too few or too many times is
//!   the wrong structure, however well each junction is supported);
//! - shallow sequence: segments far below the single-copy depth (NUMT, contamination, errors);
//! - junctions: the reads behind every link the molecules use, the weakest and any with none;
//! - complexity, as assembly graphs are classed elsewhere (simple: one segment; complex: more;
//!   multipartite: more than one connected component; nodes, links, dead ends, double forks),
//!   and backed by reads: repeats the reads show in both pairings (recombining), in one, or
//!   that no read crosses, and the chromosome count;
//! - genes: panel genes of each organelle on each molecule and the call they make;
//! - plastid: the quadripartite structure (LSC, SSC and the inverted repeat) of a circular
//!   plastid molecule.
//!
//! Nothing here re-reads the reads or calls another program; it derives from what the run wrote.

use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use rustc_hash::{FxHashMap, FxHashSet};
use serde_json::{json, Value};

use crate::evidence::{parse_gfa, Graph, Oriented};
use crate::panel::{sequence_genes, Panel, PanelParams};

/// Depth below this share of the single-copy depth: shallow (not the organelle's own copy).
const SHALLOW: f64 = 0.25;
/// Walk count and depth agree when the depth in single-copy units is within this of the count.
const COPY_TOLERANCE: f64 = 0.5;
/// Links listed as the weakest.
const WEAKEST: usize = 5;

/// One step of a molecule's path: segment and orientation (true: forward).
type Step = Oriented;

pub struct Molecule {
    pub name: String,
    pub path: Vec<Step>,
    pub circular: bool,
    pub length: u64,
}

/// "u1+ u0- u2-" -> steps, with the segment indices of `g`.
fn parse_path(path: &str, index: &FxHashMap<&str, usize>) -> Option<Vec<Step>> {
    path.split_whitespace()
        .map(|w| {
            let (name, o) = w.split_at(w.len().checked_sub(1)?);
            let fwd = match o {
                "+" => true,
                "-" => false,
                _ => return None,
            };
            Some((*index.get(name)?, fwd))
        })
        .collect()
}

/// The representative molecules of a linearization report, in its order.
pub fn molecules(linear: &Value, g: &Graph) -> Vec<Molecule> {
    let index: FxHashMap<&str, usize> =
        g.names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let mut out = Vec::new();
    for c in linear["components"].as_array().into_iter().flatten() {
        for m in c["best"]["molecules"].as_array().into_iter().flatten() {
            let Some(path) = m["path"].as_str().and_then(|p| parse_path(p, &index)) else {
                continue;
            };
            out.push(Molecule {
                name: format!("linear.{}", out.len() + 1),
                path,
                circular: m["circular"].as_bool().unwrap_or(false),
                length: m["length"].as_u64().unwrap_or(0),
            });
        }
    }
    out
}

/// Link key independent of the strand it is read on.
fn link_key(a: Step, b: Step) -> (Step, Step) {
    let rev = ((b.0, !b.1), (a.0, !a.1));
    std::cmp::min((a, b), rev)
}

fn parse_step(s: &str, index: &FxHashMap<&str, usize>) -> Option<Step> {
    parse_path(s, index).and_then(|v| v.first().copied())
}

/// Overlap of the link a -> b in the graph (0 when the graph has no such link).
fn overlap_map(g: &Graph) -> FxHashMap<(Step, Step), usize> {
    g.links
        .iter()
        .zip(&g.overlaps)
        .map(|(&(a, b), &o)| (link_key(a, b), o))
        .collect()
}

/// Bases each step adds to its molecule: the segment less its overlap with the step before.
fn step_bases(m: &Molecule, g: &Graph, overlaps: &FxHashMap<(Step, Step), usize>) -> Vec<u64> {
    let n = m.path.len();
    (0..n)
        .map(|i| {
            let len = g.seqs[m.path[i].0].len();
            let prev = if i > 0 {
                Some(m.path[i - 1])
            } else if m.circular && n > 0 {
                Some(m.path[n - 1])
            } else {
                None
            };
            let ov = prev.map_or(0, |p| overlaps.get(&link_key(p, m.path[i])).copied().unwrap_or(0));
            len.saturating_sub(ov) as u64
        })
        .collect()
}

/// Weighted median: the value at half the total weight.
fn weighted_median(mut v: Vec<(f64, u64)>) -> Option<f64> {
    v.retain(|x| x.0.is_finite() && x.1 > 0);
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.0.total_cmp(&b.0));
    let total: u64 = v.iter().map(|x| x.1).sum();
    let mut acc = 0;
    for (x, w) in &v {
        acc += w;
        if 2 * acc >= total {
            return Some(*x);
        }
    }
    v.last().map(|x| x.0)
}

/// Weakly connected components of the graph's segments.
fn components(g: &Graph) -> Vec<Vec<usize>> {
    let n = g.names.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    for &(a, b) in &g.links {
        let (ra, rb) = (find(&mut parent, a.0), find(&mut parent, b.0));
        if ra != rb {
            parent[ra] = rb;
        }
    }
    let mut groups: FxHashMap<usize, Vec<usize>> = FxHashMap::default();
    for s in 0..n {
        let r = find(&mut parent, s);
        groups.entry(r).or_default().push(s);
    }
    let mut out: Vec<Vec<usize>> = groups.into_values().collect();
    out.sort_by_key(|c| c[0]);
    out
}

/// Graph shape: dead ends (segment ends with no link) and double forks (segments with at
/// least two links at each end, the repeat nodes recombination acts on).
fn shape(g: &Graph) -> (usize, usize) {
    let n = g.names.len();
    // links at the right end (leaving the segment forward) and the left end
    let mut right: Vec<FxHashSet<Step>> = vec![FxHashSet::default(); n];
    let mut left: Vec<FxHashSet<Step>> = vec![FxHashSet::default(); n];
    for &(a, b) in &g.links {
        // a's exit end gets b; b's entry end gets a
        if a.1 { right[a.0].insert(b) } else { left[a.0].insert(b) };
        if b.1 { left[b.0].insert(a) } else { right[b.0].insert(a) };
    }
    let dead = (0..n).map(|s| usize::from(right[s].is_empty()) + usize::from(left[s].is_empty())).sum();
    let forks = (0..n).filter(|&s| right[s].len() >= 2 && left[s].len() >= 2).count();
    (dead, forks)
}

fn read_fasta(path: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    let text = fs::read_to_string(path).with_context(|| format!("cannot read {}", path.display()))?;
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            out.push((h.split_whitespace().next().unwrap_or("").to_string(), Vec::new()));
        } else if let Some(last) = out.last_mut() {
            last.1.extend(line.trim().bytes());
        }
    }
    Ok(out)
}

/// LSC, SSC and inverted repeat of a circular plastid molecule: the steps whose segment the
/// molecule walks twice in opposite orientations are the repeat; around the circle they must
/// form two runs of equal length, separating two single-copy runs.
fn quadripartite(m: &Molecule, bases: &[u64]) -> Value {
    if !m.circular {
        return json!({"status": "not_assessed", "reason": "the molecule is not circular"});
    }
    let mut orient: FxHashMap<usize, (u32, u32)> = FxHashMap::default();
    for &(s, f) in &m.path {
        let e = orient.entry(s).or_default();
        if f { e.0 += 1 } else { e.1 += 1 }
    }
    let ir: Vec<bool> = m.path.iter().map(|&(s, _)| orient[&s] == (1, 1)).collect();
    let n = ir.len();
    if !ir.iter().any(|&x| x) {
        return json!({"status": "no_inverted_repeat", "reason": "no segment is walked once in each orientation"});
    }
    if ir.iter().all(|&x| x) {
        return json!({"status": "no_single_copy", "reason": "every segment is in the repeat"});
    }
    // runs around the circle, starting after a change
    let start = (0..n).find(|&i| ir[i] != ir[(i + n - 1) % n]).unwrap_or(0);
    let mut runs: Vec<(bool, u64)> = Vec::new();
    for j in 0..n {
        let i = (start + j) % n;
        match runs.last_mut() {
            Some(r) if r.0 == ir[i] => r.1 += bases[i],
            _ => runs.push((ir[i], bases[i])),
        }
    }
    let repeats: Vec<u64> = runs.iter().filter(|r| r.0).map(|r| r.1).collect();
    let single: Vec<u64> = runs.iter().filter(|r| !r.0).map(|r| r.1).collect();
    if repeats.len() != 2 || single.len() != 2 {
        return json!({
            "status": "not_quadripartite",
            "reason": format!("{} repeat runs and {} single-copy runs around the circle", repeats.len(), single.len()),
        });
    }
    let (lsc, ssc) = (single[0].max(single[1]), single[0].min(single[1]));
    json!({
        "status": if repeats[0] == repeats[1] { "quadripartite" } else { "unequal_repeat_copies" },
        "lsc": lsc,
        "ssc": ssc,
        "ir": repeats[0].max(repeats[1]),
        "ir_copies": repeats,
        "length": lsc + ssc + repeats[0] + repeats[1],
    })
}

pub struct Inputs<'a> {
    pub gfa: &'a Path,
    pub evidence: &'a Value,
    pub linearization: &'a Value,
    pub chromosomes: Option<&'a Value>,
    /// The molecules' FASTA, for the panel genes (none: genes not assessed).
    pub molecules_fasta: Option<&'a Path>,
    pub organelle: &'a str,
}

pub fn assess(inp: &Inputs, panel: Option<&Panel>) -> Result<Value> {
    let g = parse_gfa(inp.gfa)?;
    let index: FxHashMap<&str, usize> =
        g.names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let n = g.names.len();
    let depth: Vec<Option<f64>> = g
        .names
        .iter()
        .map(|nm| inp.evidence["segment_depths"][nm].as_f64())
        .collect();
    let mols = molecules(inp.linearization, &g);
    let overlaps = overlap_map(&g);

    // walks of each segment over all molecules
    let mut walks = vec![0u32; n];
    for m in &mols {
        for &(s, _) in &m.path {
            walks[s] += 1;
        }
    }
    let graph_bases: u64 = g.seqs.iter().map(|s| s.len() as u64).sum();
    let len = |s: usize| g.seqs[s].len() as u64;

    // single-copy depth: over segments walked once; else the evidence report's median
    let unit = weighted_median(
        (0..n)
            .filter(|&s| walks[s] == 1)
            .filter_map(|s| depth[s].map(|d| (d, len(s))))
            .collect(),
    )
    .or_else(|| inp.evidence["median_depth"].as_f64())
    .filter(|&d| d > 0.0);

    // copy number against walks
    let mut agree_bases = 0u64;
    let mut walked_bases = 0u64;
    let mut disagreements = Vec::new();
    let mut shallow = Vec::new();
    let mut shallow_bases = 0u64;
    let mut unwalked_bases = 0u64;
    if let Some(u) = unit {
        for s in 0..n {
            let Some(d) = depth[s] else { continue };
            let copies = d / u;
            if copies < SHALLOW {
                shallow_bases += len(s);
                shallow.push(json!({"segment": g.names[s], "length": len(s), "depth": d, "copies": copies, "walks": walks[s]}));
            }
            if walks[s] == 0 {
                unwalked_bases += len(s);
                continue;
            }
            walked_bases += len(s);
            if (copies - f64::from(walks[s])).abs() <= COPY_TOLERANCE {
                agree_bases += len(s);
            } else {
                disagreements.push(json!({
                    "segment": g.names[s],
                    "length": len(s),
                    "depth": d,
                    "copies_by_depth": copies,
                    "walks": walks[s],
                }));
            }
        }
    }
    disagreements.sort_by(|a, b| b["length"].as_u64().cmp(&a["length"].as_u64()));
    shallow.sort_by(|a, b| b["length"].as_u64().cmp(&a["length"].as_u64()));

    // reads behind the links the molecules use
    let mut reads_of: FxHashMap<(Step, Step), (u64, u64)> = FxHashMap::default();
    for l in inp.evidence["links"].as_array().into_iter().flatten() {
        let (Some(a), Some(b)) = (
            l["from"].as_str().and_then(|x| parse_step(x, &index)),
            l["to"].as_str().and_then(|x| parse_step(x, &index)),
        ) else {
            continue;
        };
        let r = l["reads"].as_u64().unwrap_or(0) + l["adjacency_reads"].as_u64().unwrap_or(0);
        let amb = l["ambiguous_reads"].as_u64().unwrap_or(0);
        let e = reads_of.entry(link_key(a, b)).or_default();
        e.0 += r;
        e.1 += amb;
    }
    let mut used: Vec<(u64, String)> = Vec::new();
    let mut seen: FxHashSet<(Step, Step)> = FxHashSet::default();
    let step_name = |s: Step| format!("{}{}", g.names[s.0], if s.1 { '+' } else { '-' });
    for m in &mols {
        let k = m.path.len();
        let pairs = if m.circular { k } else { k.saturating_sub(1) };
        for i in 0..pairs {
            let (a, b) = (m.path[i], m.path[(i + 1) % k]);
            if k == 1 && !m.circular {
                break;
            }
            let key = link_key(a, b);
            if seen.insert(key) {
                let r = reads_of.get(&key).map_or(0, |x| x.0);
                used.push((r, format!("{} {}", step_name(a), step_name(b))));
            }
        }
    }
    used.sort();
    let zero = used.iter().filter(|x| x.0 == 0).count();
    let median_reads = if used.is_empty() { None } else { Some(used[used.len() / 2].0) };

    // complexity
    let comps = components(&g);
    let (dead_ends, double_forks) = shape(&g);
    let class = if comps.len() > 1 {
        "multipartite"
    } else if n == 1 {
        "simple"
    } else {
        "complex"
    };
    let mut repeat_calls: FxHashMap<String, usize> = FxHashMap::default();
    let mut recombining = Vec::new();
    let mut chromosome_range = Value::Null;
    if let Some(c) = inp.chromosomes {
        chromosome_range = c["range"].clone();
        for ch in c["chromosomes"].as_array().into_iter().flatten() {
            for r in ch["repeats"].as_array().into_iter().flatten() {
                let call = r["call"].as_str().unwrap_or("unknown").to_string();
                if call == "interconvertible" {
                    if let Some(p) = r["pairings"].as_array() {
                        let v: Vec<u64> = p.iter().filter_map(Value::as_u64).collect();
                        if v.len() == 2 && v[0] + v[1] > 0 {
                            recombining.push(json!({
                                "segment": r["segment"],
                                "length": r["length"],
                                "pairings": v,
                                "minor_share": v[0].min(v[1]) as f64 / (v[0] + v[1]) as f64,
                            }));
                        }
                    }
                }
                *repeat_calls.entry(call).or_default() += 1;
            }
        }
    }

    // each molecule's depth in single-copy units: its steps' depth shared among the walks of
    // their segment, weighted by the bases each step adds; a contaminant or a plastid caught in
    // a mitochondrial assembly stands far from 1
    let molecule_depth: Vec<Option<f64>> = mols
        .iter()
        .map(|m| {
            let u = unit?;
            let b = step_bases(m, &g, &overlaps);
            weighted_median(
                m.path
                    .iter()
                    .zip(&b)
                    .filter_map(|(&(s, _), &w)| depth[s].map(|d| (d / u / f64::from(walks[s].max(1)), w)))
                    .collect(),
            )
        })
        .collect();

    // genes on each molecule
    let genes = match (panel, inp.molecules_fasta) {
        (Some(panel), Some(f)) if f.exists() => {
            let p = PanelParams { min_gene_segments: 1, ..PanelParams::default() };
            let records = read_fasta(f)?;
            let same_order = records.len() == mols.len();
            let mut per = Vec::new();
            for (i, (name, seq)) in records.into_iter().enumerate() {
                let (mito, plastid) = sequence_genes([seq.as_slice()], panel, &p);
                let (target, other) = if inp.organelle == "plastid" {
                    (&plastid, &mito)
                } else {
                    (&mito, &plastid)
                };
                let other_name = if inp.organelle == "plastid" { "mitochondrion" } else { "plastid" };
                // A molecule is the target's when it holds enough of the target's genes, whatever
                // else it holds: plant mitochondria carry plastid genes (MTPTs) by the dozen. It is
                // the other organelle's only when that has enough genes and twice the target's.
                let call = if target.genes >= p.min_genes {
                    inp.organelle
                } else if other.genes >= p.min_genes && other.genes >= 2 * target.genes {
                    other_name
                } else {
                    "few_genes"
                };
                per.push(json!({
                    "molecule": name,
                    "length": seq.len(),
                    "depth_copies": if same_order { json!(molecule_depth[i]) } else { Value::Null },
                    "target_genes": target.gene_segments.iter().map(|x| x.0.clone()).collect::<Vec<_>>(),
                    "other_genes": other.gene_segments.iter().map(|x| x.0.clone()).collect::<Vec<_>>(),
                    "call": call,
                }));
            }
            json!(per)
        }
        _ => Value::Null,
    };

    // plastid structure, on the longest circular molecule
    let plastid = if inp.organelle == "plastid" {
        mols.iter()
            .filter(|m| m.circular)
            .max_by_key(|m| m.length)
            .map(|m| {
                let mut q = quadripartite(m, &step_bases(m, &g, &overlaps));
                q["molecule"] = json!(m.name);
                q
            })
            .unwrap_or_else(|| json!({"status": "not_assessed", "reason": "no circular molecule"}))
    } else {
        Value::Null
    };

    let molecule_bases: u64 = mols.iter().map(|m| step_bases(m, &g, &overlaps).iter().sum::<u64>()).sum();
    Ok(json!({
        "tool": "ovasm check assess",
        "organelle": inp.organelle,
        "scope": "derived from the run's graph, read evidence, molecules and built-in panel; no other tool",
        "depth": {
            "single_copy_depth": unit,
            "source": "length-weighted median over segments the molecules walk once",
        },
        "copy_number": {
            "walked_bases": walked_bases,
            "agreeing_bases": agree_bases,
            "agreement": (walked_bases > 0).then(|| agree_bases as f64 / walked_bases as f64),
            "tolerance": COPY_TOLERANCE,
            "disagreements": disagreements,
        },
        "shallow": {
            "below_share_of_single_copy": SHALLOW,
            "bases": shallow_bases,
            "segments": shallow,
        },
        "coverage_of_graph": {
            "graph_bases": graph_bases,
            "molecule_bases": molecule_bases,
            "unwalked_bases": unwalked_bases,
        },
        "junctions": {
            "used": used.len(),
            "without_reads": zero,
            "min_reads": used.first().map(|x| x.0),
            "median_reads": median_reads,
            "weakest": used.iter().take(WEAKEST).map(|(r, l)| json!({"link": l, "reads": r})).collect::<Vec<_>>(),
        },
        "complexity": {
            "class": class,
            "segments": n,
            "links": g.links.len(),
            "components": comps.len(),
            "dead_ends": dead_ends,
            "double_forks": double_forks,
            "repeat_calls": repeat_calls,
            "recombining_repeats": recombining,
            "chromosome_range": chromosome_range,
        },
        "genes": genes,
        "plastid_structure": plastid,
    }))
}

/// One line for logs.
pub fn summary(r: &Value) -> String {
    format!(
        "copy-number agreement {}; {} links used, {} without reads (min {}); {} {}",
        r["copy_number"]["agreement"].as_f64().map_or("n/a".into(), |a| format!("{:.3}", a)),
        r["junctions"]["used"],
        r["junctions"]["without_reads"],
        r["junctions"]["min_reads"],
        r["complexity"]["class"].as_str().unwrap_or("?"),
        match r["plastid_structure"]["status"].as_str() {
            Some(s) => format!("plastid {s}"),
            None => String::new(),
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write(dir: &Path, name: &str, text: &str) -> std::path::PathBuf {
        let p = dir.join(name);
        fs::write(&p, text).unwrap();
        p
    }

    fn seq(n: usize, seed: u64) -> String {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4] as char
            })
            .collect()
    }

    /// A plastid graph like the Carex one: LSC u1, IR u0 (twice the depth), SSC u2.
    fn plastid_case(dir: &Path, ir_depth: f64, ssc_reads: u64) -> (std::path::PathBuf, Value, Value) {
        let gfa = format!(
            "H\tVN:Z:1.0\nS\tu0\t{}\tdp:f:{ir_depth}\nS\tu1\t{}\tdp:f:100\nS\tu2\t{}\tdp:f:95\n\
             L\tu0\t-\tu2\t-\t0M\nL\tu0\t-\tu2\t+\t0M\nL\tu0\t+\tu1\t-\t0M\nL\tu0\t+\tu1\t+\t0M\n",
            seq(300, 1),
            seq(1000, 2),
            seq(100, 3)
        );
        let g = write(dir, "a.gfa", &gfa);
        let evidence = json!({
            "segment_depths": {"u0": ir_depth, "u1": 100.0, "u2": 95.0},
            "median_depth": 100.0,
            "links": [
                {"from": "u0-", "to": "u2-", "reads": ssc_reads},
                {"from": "u0-", "to": "u2+", "reads": 40},
                {"from": "u0+", "to": "u1-", "reads": 30},
                {"from": "u0+", "to": "u1+", "reads": 35},
            ],
        });
        let linear = json!({"components": [{"best": {"molecules": [
            {"path": "u1+ u0- u2- u0+", "circular": true, "length": 1700}
        ]}}]});
        (g, evidence, linear)
    }

    #[test]
    fn a_plastid_with_its_repeat_at_twice_the_depth_agrees_and_is_quadripartite() {
        let d = tempfile_dir("assess-pt");
        let (g, ev, lin) = plastid_case(&d, 200.0, 25);
        let r = assess(&Inputs { gfa: &g, evidence: &ev, linearization: &lin, chromosomes: None,
            molecules_fasta: None, organelle: "plastid" }, None).unwrap();
        assert_eq!(r["copy_number"]["agreement"].as_f64(), Some(1.0));
        assert_eq!(r["plastid_structure"]["status"], "quadripartite");
        assert_eq!(r["plastid_structure"]["lsc"], 1000);
        assert_eq!(r["plastid_structure"]["ssc"], 100);
        assert_eq!(r["plastid_structure"]["ir"], 300);
        assert_eq!(r["junctions"]["used"], 4);
        assert_eq!(r["junctions"]["without_reads"], 0);
        assert_eq!(r["junctions"]["min_reads"], 25);
        assert_eq!(r["complexity"]["class"], "complex");
        assert_eq!(r["complexity"]["double_forks"], 1);
        assert_eq!(r["complexity"]["dead_ends"], 0);
    }

    #[test]
    fn a_repeat_at_single_copy_depth_walked_twice_disagrees_and_an_unread_link_counts() {
        let d = tempfile_dir("assess-bad");
        let (g, ev, lin) = plastid_case(&d, 100.0, 0);
        let r = assess(&Inputs { gfa: &g, evidence: &ev, linearization: &lin, chromosomes: None,
            molecules_fasta: None, organelle: "plastid" }, None).unwrap();
        let a = r["copy_number"]["agreement"].as_f64().unwrap();
        assert!((a - 1100.0 / 1400.0).abs() < 1e-9, "{a}");
        assert_eq!(r["copy_number"]["disagreements"][0]["segment"], "u0");
        assert_eq!(r["junctions"]["without_reads"], 1);
        assert_eq!(r["junctions"]["weakest"][0]["reads"], 0);
    }

    #[test]
    fn separate_components_are_multipartite_and_shallow_segments_are_listed() {
        let d = tempfile_dir("assess-multi");
        let gfa = format!("H\tVN:Z:1.0\nS\ta\t{}\nS\tb\t{}\nS\tc\t{}\nL\ta\t+\ta\t+\t0M\nL\tb\t+\tb\t+\t0M\n",
            seq(500, 4), seq(400, 5), seq(200, 6));
        let g = write(&d, "m.gfa", &gfa);
        let ev = json!({"segment_depths": {"a": 50.0, "b": 52.0, "c": 3.0}, "median_depth": 50.0,
            "links": [{"from": "a+", "to": "a+", "reads": 12}, {"from": "b+", "to": "b+", "reads": 9}]});
        let lin = json!({"components": [
            {"best": {"molecules": [{"path": "a+", "circular": true, "length": 500}]}},
            {"best": {"molecules": [{"path": "b+", "circular": true, "length": 400}]}}]});
        let r = assess(&Inputs { gfa: &g, evidence: &ev, linearization: &lin, chromosomes: None,
            molecules_fasta: None, organelle: "mitochondrion" }, None).unwrap();
        assert_eq!(r["complexity"]["class"], "multipartite");
        assert_eq!(r["complexity"]["components"], 3);
        assert_eq!(r["shallow"]["segments"][0]["segment"], "c");
        assert_eq!(r["coverage_of_graph"]["unwalked_bases"], 200);
        assert_eq!(r["junctions"]["used"], 2);
        assert_eq!(r["junctions"]["min_reads"], 9);
        assert!(r["plastid_structure"].is_null());
    }

    fn tempfile_dir(name: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("ovasm-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }
}
