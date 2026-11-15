//! Stage 2b: unringing, representative molecules from read-supported paths (`ovasm solve linearize`).
//!
//! Unique segments (those the evidence report does not list as repeats) are anchors. Two anchors
//! are bridged when a read-supported phased path runs from one to the other through repeats, or
//! when a supported link joins them directly; a bridge back to its own anchor (a self-link, or a
//! path through repeats) can close a one-anchor circle. The assembly graph holds every
//! conformation at once; unringing picks molecules out of it. Anchors that no bridge joins cannot
//! share a molecule, so each connected component of the anchor graph is solved on its own
//! (several chromosomes come out as several components). Within a component, a cover is a set of
//! molecules that together visit every anchor exactly once through bridges (with one molecule, a
//! Hamiltonian path in the small anchor graph); repeats are laid down as often as the chosen
//! bridges cross them. A molecule is circular when a bridge leads from its last anchor back to
//! its first. Covers of every permitted size compete together; `cover_order` ranks them by
//! circularity, junction evidence, and `path_score`. `path_score` is the likelihood
//! of the paths when each anchor's exit is
//! chosen in proportion to
//! the reads on its bridges. A best path whose end leads back into its own interior is
//! reported as the circle and its tail (`cut_tail`). An anchor without any bridge is a
//! linear molecule of its own when it is at least `MIN_UNBRIDGED_BP` long.
//!
//! What no cover explains is not dropped: a component without a cover, and every required
//! anchor without a bridge that is too short to be a molecule, comes out as partial paths
//! (`partial_paths`), in the report and in a FASTA of their own. They are not molecules and
//! the component stays unsolved.

use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

use anyhow::{bail, Context, Result};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use serde_json::Value;

use crate::evidence::{flip, parse_gfa, Graph, Oriented};

#[cfg(test)]
#[path = "multichrom_diagnostics.rs"]
mod multichrom_diagnostics;

#[derive(Clone, Debug)]
struct Bridge {
    /// Walk from one anchor to another; everything between them is a repeat.
    walk: Vec<Oriented>,
    reads: u64,
    /// From the graph alone, through repeats no read crosses (see `bridges`).
    topological: bool,
    /// Taken at an anchor end no bridge of `min_reads` leaves: every pairing the graph allows
    /// there (see `bridges`); a molecule through it is not decisive.
    weak: bool,
}

/// Reads a topological bridge counts as: all pairings through an uncrossed repeat alike.
const TOPO_READS: u64 = 1;
/// Longest run of uncrossed repeats a topological bridge passes.
const TOPO_MAX_REPEATS: usize = 8;
/// Reads crossing repeats from an anchor end below which a pairing no read took is not evidence
/// that it is absent: the end also gets the graph's other pairings (topological).
const ABSENCE_READS: u64 = 10;
/// Penalty (ln) per repeat copy a cover lays down more or fewer times than the repeat's depth
/// says (see `copy_deviation`): a factor e^-2 per copy.
const COPY_PENALTY_LN: f64 = 2.0;
/// Shortest anchor without any bridge reported as a molecule of its own (see `run`).
const MIN_UNBRIDGED_BP: usize = 1000;

#[derive(Clone, Copy)]
pub struct LinearizeParams {
    /// Bridges supported by fewer reads are ignored.
    pub min_reads: u64,
    /// Upper bound on search expansions per component and cover size (the anchor graph is small,
    /// but stay bounded).
    pub max_steps: u64,
    /// Automatic mode may retry a small, near-complete search once with a bounded budget.
    pub adaptive_steps: bool,
    /// Explicit cover-size bound; zero uses the existing four-cover search envelope
    /// raised to the component's necessary path-end lower bound.
    pub max_molecules: usize,
    pub max_alternatives: usize,
    /// An anchor shallower than this fraction of the median depth is optional: a low-frequency
    /// configuration the representative path need not visit. In a component shallower than
    /// the graph's anchors as a whole the threshold comes down in proportion, unless the
    /// component is itself below the fraction (see `run`). Without depth in the evidence,
    /// anchors below half this fraction of the median bridge support are optional.
    pub optional_below: f64,
}

/// A bridge as `bridges` builds it, for a caller that weighs bridges from outside the reads
/// (`open_bridges`).
pub(crate) struct BridgeView {
    pub walk: Vec<Oriented>,
    pub reads: u64,
    pub topological: bool,
    pub weak: bool,
}

/// Junction preferences from outside the target's reads (`pan reconstruct`: a species panel):
/// per bridge walk, the probability that the molecule leaves the walk's first anchor that way.
/// They order only covers the reads rank alike (`pan_posterior::best_matching`).
pub struct BridgePrior {
    pub probability: FxHashMap<Vec<Oriented>, f64>,
}

/// The cover of a component a `BridgePrior` takes among those the reads rank alike: every
/// anchor end a bridge leaves is on one bridge.
#[derive(Clone, Serialize)]
pub struct PriorCover {
    pub cover: CoverReport,
    /// Anchor ends of the component that a bridge leaves, those with more than one bridge to
    /// take, and those of them where the prior tells the bridges apart.
    pub ends: usize,
    pub open_ends: usize,
    pub decided_ends: usize,
    /// Sum over the ends of the prior probability of the bridge taken there: the expected number
    /// of correctly joined ends.
    pub expected_correct_ends: f64,
    /// Covers the reads and the prior rank as well (the one with the fewest molecules is given).
    pub equally_good: u64,
    pub matchings_examined: u64,
    /// The target's own single records (one read through three segments or more) the cover
    /// explains: those through two repeats or more, and all. They were weighed before the prior.
    pub joint_records_explained: u64,
    pub thread_records_explained: u64,
}

fn rc_walk(w: &[Oriented]) -> Vec<Oriented> {
    w.iter().rev().map(|&o| flip(o)).collect()
}

fn parse_label(index: &FxHashMap<&str, usize>, label: &str) -> Result<Oriented> {
    let Some((name, sign)) = label.len().checked_sub(1).map(|i| label.split_at(i)) else {
        bail!("empty segment label");
    };
    let fwd = match sign {
        "+" => true,
        "-" => false,
        _ => bail!("segment label {label:?} must end in + or -"),
    };
    let s = index
        .get(name)
        .copied()
        .with_context(|| format!("evidence names segment {name:?}, which is not in the GFA"))?;
    Ok((s, fwd))
}


type ClosureEdge = (Oriented, Oriented);

fn closure_edge(from: Oriented, to: Oriented) -> ClosureEdge {
    (from, to).min((flip(to), flip(from)))
}

#[derive(Serialize)]
pub struct LocalRepeatClosureLink {
    from: String,
    to: String,
    overlaps_bp: BTreeSet<usize>,
    single_record_crossings: Option<u64>,
}

#[derive(Serialize)]
pub struct LocalRepeatClosure {
    segment: String,
    link: LocalRepeatClosureLink,
    state: &'static str,
    closed_walk_length: Option<usize>,
    competing_links: Vec<LocalRepeatClosureLink>,
}

#[derive(Serialize)]
pub struct LocalRepeatClosureContext {
    source: &'static str,
    path: String,
    reported_records: u64,
    closure_segments: Vec<String>,
}

#[derive(Serialize)]
pub struct LocalRepeatClosures {
    scope: &'static str,
    min_single_record_crossings: u64,
    contexts_truncated: bool,
    closures: Vec<LocalRepeatClosure>,
    contexts: Vec<LocalRepeatClosureContext>,
}

fn local_repeat_closures(
    g: &Graph,
    report: &Value,
    repeats: &[bool],
    graph_bound_records: bool,
    phases: &crate::phase::Constraints,
    min_reads: u64,
    context_budget: u64,
) -> Option<LocalRepeatClosures> {
    let seeds: BTreeSet<_> = g.links.iter()
        .filter_map(|&(a, b)| (a == b && repeats[a.0]).then_some(a.0))
        .collect();
    if seeds.is_empty() {
        return None;
    }
    let index: FxHashMap<_, _> = g.names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let mut edges: BTreeMap<ClosureEdge, BTreeSet<usize>> = BTreeMap::new();
    let mut incident: BTreeMap<usize, BTreeSet<ClosureEdge>> = BTreeMap::new();
    for (&(a, b), &overlap) in g.links.iter().zip(&g.overlaps) {
        let key = closure_edge(a, b);
        edges.entry(key).or_default().insert(overlap);
        for s in [a.0, b.0] {
            if seeds.contains(&s) {
                incident.entry(s).or_default().insert(key);
            }
        }
    }
    let mut record_links = BTreeMap::new();
    if graph_bound_records {
        // ReadModel::load has already checked the graph, schema and uniqueness of
        // these canonical links. Its producer counts each record once per link.
        for row in report["thread_evidence"]["links"].as_array().into_iter().flatten() {
            let from = parse_label(&index, row["from"].as_str().unwrap()).unwrap();
            let to = parse_label(&index, row["to"].as_str().unwrap()).unwrap();
            record_links.insert(closure_edge(from, to), row["records"].as_u64().unwrap());
        }
    }
    let describe_link = |key: ClosureEdge| LocalRepeatClosureLink {
        from: g.label(key.0),
        to: g.label(key.1),
        overlaps_bp: edges[&key].clone(),
        single_record_crossings: graph_bound_records.then(|| record_links.get(&key).copied().unwrap_or(0)),
    };
    let min_reads = min_reads.max(1);
    let closures = seeds.iter().map(|&s| {
        let node = (s, true);
        let key = closure_edge(node, node);
        let link = describe_link(key);
        let overlaps = &link.overlaps_bp;
        let seq = &g.seqs[s];
        let valid_overlap = overlaps.iter().next().copied().filter(|&n| {
            overlaps.len() == 1 && n < seq.len() && seq[..n] == seq[seq.len() - n..]
        });
        let competing_links = incident[&s].iter().copied().filter(|&e| e != key)
            .map(&describe_link).collect::<Vec<_>>();
        let (state, closed_walk_length) = match valid_overlap {
            None => ("unresolved_overlap", None),
            Some(_) if !phases.allows_walk(&[node], true) => ("unresolved_phase_conflict", None),
            Some(_) if link.single_record_crossings.is_none() => ("unresolved_no_graph_bound_crossings", None),
            Some(_) if link.single_record_crossings.unwrap_or(0) < min_reads => ("unresolved_insufficient_crossings", None),
            Some(n) if !competing_links.is_empty() => ("unresolved_competing_context", Some(seq.len() - n)),
            Some(n) => ("unresolved_local_closed_walk", Some(seq.len() - n)),
        };
        LocalRepeatClosure { segment: g.names[s].clone(), link, state, closed_walk_length, competing_links }
    }).collect();

    // Store each canonical context once, even when it touches several closure
    // segments. The separate path-scan budget cannot consume cover-search work.
    let mut context_index = BTreeMap::<(&'static str, Vec<Oriented>), u64>::new();
    let mut remaining = context_budget;
    let mut contexts_truncated = false;
    for (&(entry, segment, exit), &records) in &phases.phases {
        let walk = vec![entry, segment, exit];
        if !walk.iter().any(|o| seeds.contains(&o.0)) {
            continue;
        }
        let Some(next) = remaining.checked_sub(3) else {
            contexts_truncated = true;
            break;
        };
        remaining = next;
        let reverse = rc_walk(&walk);
        let support = context_index.entry(("phase_report", walk.min(reverse))).or_default();
        *support = (*support).max(records);
    }
    let inputs = [
        ("single_record_thread", &report["thread_evidence"]["paths"], graph_bound_records),
        ("observed_walk", &report["observed_walks"], true),
        ("phased_path", &report["phased_paths"], true),
    ];
    'sources: for (source, rows, available) in inputs {
        if contexts_truncated || !available {
            continue;
        }
        for row in rows.as_array().into_iter().flatten() {
            let (Some(path), Some(records)) = (row["path"].as_str(), row["reads"].as_u64()) else {
                continue;
            };
            let mut walk = Vec::new();
            for token in path.split_whitespace() {
                let Some(next) = remaining.checked_sub(1) else {
                    contexts_truncated = true;
                    break 'sources;
                };
                remaining = next;
                let Ok(node) = parse_label(&index, token) else {
                    walk.clear();
                    break;
                };
                walk.push(node);
            }
            if walk.len() < 2 || !walk.iter().any(|o| seeds.contains(&o.0))
                || walk.windows(2).any(|w| !edges.contains_key(&closure_edge(w[0], w[1])))
            {
                continue;
            }
            let reverse = rc_walk(&walk);
            let entry = context_index.entry((source, walk.min(reverse))).or_default();
            *entry = (*entry).max(records);
        }
    }
    let contexts = context_index.into_iter().map(|((source, walk), reported_records)| {
        let touched: BTreeSet<_> = walk.iter().map(|o| o.0).filter(|s| seeds.contains(s)).collect();
        LocalRepeatClosureContext {
            source,
            path: walk.iter().map(|&o| g.label(o)).collect::<Vec<_>>().join(" "),
            reported_records,
            closure_segments: touched.into_iter().map(|s| g.names[s].clone()).collect(),
        }
    }).collect();
    Some(LocalRepeatClosures {
        scope: "local repeat self-loop spellings; no complete molecule or physical topology inference",
        min_single_record_crossings: min_reads,
        contexts_truncated,
        closures,
        contexts,
    })
}

/// Bridges in both directions: phased paths plus direct links between two anchors.
fn bridges(
    g: &Graph,
    report: &Value,
    repeat: &[bool],
    min_reads: u64,
) -> Result<(Vec<Bridge>, Vec<String>)> {
    let index: FxHashMap<&str, usize> = g
        .names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    let mut best: FxHashMap<Vec<Oriented>, u64> = FxHashMap::default();
    // read-supported walks below `min_reads`, the fallback for ends nothing else leaves
    let mut weak: FxHashMap<Vec<Oriented>, u64> = FxHashMap::default();
    let mut self_loops = Vec::new();
    // repeats some read crosses (inside a phased path)
    let mut crossed = vec![false; g.names.len()];
    let mut add = |walk: Vec<Oriented>, reads: u64| {
        for w in [rc_walk(&walk), walk] {
            let e = best.entry(w).or_default();
            *e = (*e).max(reads);
        }
    };
    for p in report["phased_paths"]
        .as_array()
        .context("evidence has no phased_paths")?
    {
        let reads = p["reads"].as_u64().unwrap_or(0);
        let path = p["path"].as_str().context("phased path without a path")?;
        let walk = path
            .split_whitespace()
            .map(|l| parse_label(&index, l))
            .collect::<Result<Vec<_>>>()?;
        if walk.len() < 3 {
            continue;
        }
        if reads < min_reads {
            if reads > 0 && !repeat[walk[0].0] && !repeat[walk[walk.len() - 1].0] {
                for w in [rc_walk(&walk), walk] {
                    let e = weak.entry(w).or_default();
                    *e = (*e).max(reads);
                }
            }
            continue;
        }
        for o in &walk[1..walk.len() - 1] {
            crossed[o.0] = true;
        }
        if walk[0].0 == walk[walk.len() - 1].0 {
            self_loops.push(path.to_string());
        }
        add(walk, reads);
    }
    for l in report["links"]
        .as_array()
        .context("evidence has no links")?
    {
        // Exact dense-edge words establish adjacency between anchors. They never
        // enter phased_paths or authorize an entry/exit pairing through a repeat.
        // max is a conservative lower bound because the two read sets can overlap.
        let reads = l["reads"]
            .as_u64()
            .unwrap_or(0)
            .max(l["adjacency_reads"].as_u64().unwrap_or(0));
        let a = parse_label(&index, l["from"].as_str().context("link without from")?)?;
        let b = parse_label(&index, l["to"].as_str().context("link without to")?)?;
        if repeat[a.0] || repeat[b.0] {
            continue;
        }
        if reads < min_reads {
            if reads > 0 {
                for w in [vec![flip(b), flip(a)], vec![a, b]] {
                    let e = weak.entry(w).or_default();
                    *e = (*e).max(reads);
                }
            }
            continue;
        }
        if a.0 == b.0 {
            self_loops.push(format!("{} {}", g.label(a), g.label(b)));
        }
        add(vec![a, b], reads);
    }
    // Repeats no read crosses (short reads; or a repeat longer than every read): the reads
    // cannot tell which entry pairs with which exit, so every pairing the graph allows is a
    // bridge, weighted alike, and a representative through them comes out non-decisive. Where
    // reads cross a repeat, only their pairings count. An anchor end with no bridge of
    // `min_reads` takes, in order of evidence: the walks reads did take from it below
    // `min_reads` (weak bridges, below), else, if no read took any way from it, the graph's
    // pairings through the repeats beyond it (topological). Either way a molecule through it
    // is not decisive. Without this a long repeat crossed by few reads leaves an end dead and
    // a copy out: Antitrichia's plastid IR (1 read per closing pairing), Holcus's (whose
    // closing pairing has a strong bridge at its other end), Geg-14's 18.9 kb mitochondrial
    // repeat (6 reads over four pairings; the 356,239 bp circle came out 337,301 bp linear).
    // Reads crossing repeats from each end. Below `ABSENCE_READS` a pairing nobody read may
    // just not have been sampled: Sha_Ara-1 thinned to 15x had 9, 5 and 3 reads on three of the
    // four pairings around its 4.2 kb repeat, and the master circle needs the fourth; without it
    // the molecule came out linear with one copy of the repeat, 4,193 bp short.
    let mut crossing: FxHashMap<Oriented, u64> = FxHashMap::default();
    for (w, &r) in &best {
        if w.len() >= 3 {
            *crossing.entry(w[0]).or_default() += r;
        }
    }
    let sparse = |o: Oriented| crossing.get(&o).is_some_and(|&r| r < ABSENCE_READS);
    let mut next: FxHashMap<Oriented, Vec<Oriented>> = FxHashMap::default();
    for &(x, y) in &g.links {
        next.entry(x).or_default().push(y);
        next.entry(flip(y)).or_default().push(flip(x));
    }
    let mut topo: BTreeSet<Vec<Oriented>> = BTreeSet::new();
    for a in (0..g.names.len()).filter(|&a| !repeat[a]) {
        for start in [(a, true), (a, false)] {
            let mut stack = vec![vec![start]];
            while let Some(walk) = stack.pop() {
                let at = walk[walk.len() - 1];
                for &nx in next.get(&at).map(Vec::as_slice).unwrap_or(&[]) {
                    if !repeat[nx.0] {
                        if walk.len() >= 2 {
                            let mut w = walk.clone();
                            w.push(nx);
                            if !best.contains_key(&w) {
                                if nx.0 == a {
                                    self_loops.push(
                                        w.iter().map(|&o| g.label(o)).collect::<Vec<_>>().join(" "),
                                    );
                                }
                                topo.insert(rc_walk(&w));
                                topo.insert(w);
                            }
                        }
                    } else if !crossed[nx.0] && walk.len() <= TOPO_MAX_REPEATS {
                        let mut w = walk.clone();
                        w.push(nx);
                        stack.push(w);
                    }
                }
            }
        }
    }
    let mut out: Vec<Bridge> = best
        .into_iter()
        .map(|(walk, reads)| Bridge {
            walk,
            reads,
            topological: false,
            weak: false,
        })
        .chain(topo.into_iter().map(|walk| Bridge {
            walk,
            reads: TOPO_READS,
            topological: true,
            weak: false,
        }))
        .collect();
    // An anchor end that no bridge leaves can only end a linear molecule, and every repeat
    // beyond it is laid down once too few, or the anchor is left out. Where reads cross a
    // repeat, its other pairings get no topological bridge above, so an end whose own pairings
    // no read (or fewer than `min_reads`) supports is dead although the graph continues: the
    // Antitrichia and Holcus plastids lost an IR copy (the closing pairing read once), and the
    // 117.7 kb Buddleja mitochondrial unitig between both ends of an 18.2 kb repeat no read
    // spans was left out of a 404 kb circle (the genome is 540 kb). The reads cannot tell how
    // such an end pairs, so, as for an uncrossed repeat, every walk the graph allows from it
    // through repeats is a bridge, weighted by the reads that follow it (at least
    // `TOPO_READS`); a molecule through one is not decisive. A direct link between two anchors
    // is different: the evidence measured that junction itself, so it counts only where some
    // read supports it.
    let leaving: BTreeSet<Oriented> = out.iter().map(|b| b.walk[0]).collect();
    let existing: FxHashSet<Vec<Oriented>> = out.iter().map(|b| b.walk.clone()).collect();
    let mut fallback: BTreeMap<Vec<Oriented>, u64> = BTreeMap::new();
    for a in (0..g.names.len()).filter(|&a| !repeat[a]) {
        for start in [(a, true), (a, false)] {
            if leaving.contains(&start) {
                continue;
            }
            let mut stack = vec![vec![start]];
            while let Some(walk) = stack.pop() {
                let at = walk[walk.len() - 1];
                for &nx in next.get(&at).map(Vec::as_slice).unwrap_or(&[]) {
                    let mut w = walk.clone();
                    w.push(nx);
                    if repeat[nx.0] {
                        if walk.len() <= TOPO_MAX_REPEATS {
                            stack.push(w);
                        }
                        continue;
                    }
                    let read = weak.get(&w).copied().unwrap_or(0);
                    if walk.len() == 1 && read == 0 {
                        continue;
                    }
                    let reads = read.max(TOPO_READS);
                    for v in [rc_walk(&w), w] {
                        if !existing.contains(&v) {
                            let e = fallback.entry(v).or_default();
                            *e = (*e).max(reads);
                        }
                    }
                }
            }
        }
    }
    out.extend(fallback.into_iter().map(|(walk, reads)| Bridge {
        topological: !weak.contains_key(&walk),
        walk,
        reads,
        weak: true,
    }));
    // Ends with few crossing reads (see `ABSENCE_READS`): the pairings nobody read through the
    // repeats beyond them may just not have been sampled, so the graph's pairings join as
    // topological bridges. This runs after the dead-end fallback so that it never takes an end
    // away from the read-weighted weak bridges there; it only adds walks not yet bridged.
    let have: FxHashSet<Vec<Oriented>> = out.iter().map(|b| b.walk.clone()).collect();
    let mut extra: BTreeMap<Vec<Oriented>, u64> = BTreeMap::new();
    for a in (0..g.names.len()).filter(|&a| !repeat[a]) {
        for start in [(a, true), (a, false)] {
            if !sparse(start) {
                continue;
            }
            let mut stack = vec![vec![start]];
            while let Some(walk) = stack.pop() {
                let at = walk[walk.len() - 1];
                for &nx in next.get(&at).map(Vec::as_slice).unwrap_or(&[]) {
                    let mut w = walk.clone();
                    w.push(nx);
                    if repeat[nx.0] {
                        if walk.len() <= TOPO_MAX_REPEATS {
                            stack.push(w);
                        }
                    } else if walk.len() >= 2 {
                        let reads = weak.get(&w).copied().unwrap_or(0).max(TOPO_READS);
                        for v in [rc_walk(&w), w] {
                            if !have.contains(&v) {
                                let e = extra.entry(v).or_default();
                                *e = (*e).max(reads);
                            }
                        }
                    }
                }
            }
        }
    }
    out.extend(extra.into_iter().map(|(walk, reads)| {
        let measured = weak.contains_key(&walk);
        Bridge {
            topological: !measured,
            walk,
            reads,
            weak: measured,
        }
    }));
    out.sort_by(|x, y| x.walk.cmp(&y.walk));
    Ok((out, self_loops))
}

/// Most walks `spanned_walks` tries from one anchor end.
const SPANNED_MAX_STEPS: usize = 100_000;

/// Walks through repeats that reads span from anchor to anchor, by length (the evidence's
/// `anchor_spans`): the walks from one anchor end through repeats to the other whose bases
/// between the two anchors (`inner_length`) are what at least `min_reads` reads hold there,
/// within the evidence's `spacing_tol` (at least 50 bp), with those reads. Crossings say which
/// segment follows which, never how often a loop is passed: in a tandem array every copy's
/// junction is alike, so the topological bridges go round it up to `TOPO_MAX_REPEATS` times and
/// depth says how often (`copy_deviation`), and both fall short where a k-mer outspans the unit.
/// Sarcophaga caerulescens (PMAT2 benchmark): its control region's 455 bp unit forms a loop of
/// two segments at 178-212x against 49.7x, depth said 4 copies of each where the 49 reads
/// across the array, 3,489-3,544 bp between its flanks, say 5, and the circle came out
/// 20,787 bp of 21,241.
fn spanned_walks(
    g: &Graph,
    report: &Value,
    repeat: &[bool],
    min_reads: u64,
) -> Result<Vec<(Vec<Oriented>, u64)>> {
    let Some(rows) = report.get("anchor_spans").and_then(Value::as_array) else {
        return Ok(Vec::new());
    };
    let tol = report["params"]["spacing_tol"].as_f64().unwrap_or(0.05);
    let index: FxHashMap<&str, usize> = g
        .names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    let overlaps = overlap_index(g);
    let mut next: FxHashMap<Oriented, Vec<Oriented>> = FxHashMap::default();
    for &(x, y) in &g.links {
        next.entry(x).or_default().push(y);
        next.entry(flip(y)).or_default().push(flip(x));
    }
    for v in next.values_mut() {
        v.sort_unstable();
        v.dedup();
    }
    let max_ov = g.overlaps.iter().copied().max().unwrap_or(0) as i64;
    let fits = |span: i64, inner: i64| (span - inner).abs() as f64 <= (tol * inner.abs() as f64).max(50.0);
    let mut out = Vec::new();
    for row in rows {
        let x = parse_label(&index, row["from"].as_str().context("anchor span without from")?)?;
        let y = parse_label(&index, row["to"].as_str().context("anchor span without to")?)?;
        let spans: Vec<i64> = row["spans"]
            .as_array()
            .context("anchor span without spans")?
            .iter()
            .filter_map(Value::as_i64)
            .collect();
        if (spans.len() as u64) < min_reads || repeat[x.0] || repeat[y.0] {
            continue;
        }
        let longest = *spans.iter().max().expect("non-empty");
        let bound = longest + (tol * longest as f64).max(50.0) as i64;
        // (walk so far, bases after x's end up to the end of its last segment)
        let mut stack: Vec<(Vec<Oriented>, i64)> = vec![(vec![x], 0)];
        let mut steps = 0;
        while let Some((walk, inner)) = stack.pop() {
            steps += 1;
            if steps > SPANNED_MAX_STEPS {
                break;
            }
            let at = walk[walk.len() - 1];
            for &nx in next.get(&at).map(Vec::as_slice).unwrap_or(&[]) {
                let ov = overlaps[&(at, nx)] as i64;
                if !repeat[nx.0] {
                    if nx != y || walk.len() < 2 {
                        continue;
                    }
                    let full = inner - ov;
                    let n = spans.iter().filter(|&&s| fits(s, full)).count() as u64;
                    if n >= min_reads {
                        let mut w = walk.clone();
                        w.push(nx);
                        out.push((w, n));
                    }
                    continue;
                }
                let grown = inner + g.seqs[nx.0].len() as i64 - ov;
                if grown - max_ov <= bound {
                    let mut w = walk.clone();
                    w.push(nx);
                    stack.push((w, grown));
                }
            }
        }
    }
    out.sort();
    Ok(out)
}

/// Bridges with `spanned_walks` in: each spanned walk (both directions) becomes a read-measured
/// bridge, and a topological or weak bridge between the same two anchor ends through repeats
/// that no span fits leaves; returns the canonical spanned walks.
fn merge_spanned(br: &mut Vec<Bridge>, spanned: &[(Vec<Oriented>, u64)]) -> Vec<Vec<Oriented>> {
    if spanned.is_empty() {
        return Vec::new();
    }
    let mut walks: FxHashMap<Vec<Oriented>, u64> = FxHashMap::default();
    let mut ends: FxHashSet<(Oriented, Oriented)> = FxHashSet::default();
    let mut canonical = Vec::new();
    for (w, n) in spanned {
        let rc = rc_walk(w);
        canonical.push(w.clone().min(rc.clone()));
        for v in [rc, w.clone()] {
            ends.insert((v[0], v[v.len() - 1]));
            let e = walks.entry(v).or_default();
            *e = (*e).max(*n);
        }
    }
    br.retain(|b| {
        let unmeasured = b.topological || b.weak;
        !(unmeasured && b.walk.len() >= 3 && ends.contains(&(b.walk[0], b.walk[b.walk.len() - 1])))
            || walks.contains_key(&b.walk)
    });
    for b in br.iter_mut() {
        if let Some(&n) = walks.get(&b.walk) {
            b.reads = b.reads.max(n);
            b.topological = false;
            b.weak = false;
        }
    }
    let have: FxHashSet<Vec<Oriented>> = br.iter().map(|b| b.walk.clone()).collect();
    for (w, &n) in &walks {
        if !have.contains(w) {
            br.push(Bridge {
                walk: w.clone(),
                reads: n,
                topological: false,
                weak: false,
            });
        }
    }
    br.sort_by(|x, y| x.walk.cmp(&y.walk));
    canonical.sort();
    canonical.dedup();
    canonical
}

/// Best path must be at least this much more likely than the runner-up (ln 2: twice as likely)
/// to be reported as decisive; otherwise co-existing configurations are about equally supported.
const DECISIVE_LN: f64 = std::f64::consts::LN_2;

/// Share of each bridge among its rival bridges. Raw read counts mix how common a configuration
/// is with how many reads are long enough to span the bridge and with sequencing depth; among the
/// bridges at one anchor end those factors largely cancel. A bridge competes at both of its ends
/// (leaving its first anchor, entering its last), so the share is the geometric mean of the two,
/// which also makes a path score the same read in either direction.
///
/// Topological bridges never dilute what reads measured: at an end with read bridges those
/// share among themselves as without the graph's pairings, and a topological bridge leaving
/// that end gets the share of one more read. (A topological bridge added for an end without
/// a read bridge of its own reaches, reversed, ends that have them.)
fn exit_fractions(br: &[Bridge]) -> Vec<f64> {
    let mut total: FxHashMap<Oriented, u64> = FxHashMap::default();
    let mut read_total: FxHashMap<Oriented, u64> = FxHashMap::default();
    for b in br {
        *total.entry(b.walk[0]).or_default() += b.reads;
        if !b.topological {
            *read_total.entry(b.walk[0]).or_default() += b.reads;
        }
    }
    let pos: FxHashMap<&[Oriented], usize> = br
        .iter()
        .enumerate()
        .map(|(i, b)| (b.walk.as_slice(), i))
        .collect();
    let out: Vec<f64> = br
        .iter()
        .map(|b| {
            let reads = read_total.get(&b.walk[0]).copied().unwrap_or(0);
            let denominator = match (reads, b.topological) {
                (0, _) => total[&b.walk[0]],
                (r, false) => r,
                (r, true) => r + b.reads,
            };
            b.reads as f64 / denominator.max(1) as f64
        })
        .collect();
    br.iter()
        .enumerate()
        .map(|(i, b)| {
            // `bridges` always adds the reverse complement, which leaves from the far end.
            let j = pos[rc_walk(&b.walk).as_slice()];
            (out[i] * out[j]).sqrt()
        })
        .collect()
}

/// Log-likelihood of a candidate from the exit fractions of its bridges in path order (a circular
/// candidate includes the bridge that closes it). Higher is better.
fn path_score(fractions: &[f64]) -> f64 {
    fractions.iter().map(|f| f.max(1e-9).ln()).sum()
}

/// One molecule of a cover.
#[derive(Clone)]
struct Found {
    start: Oriented,
    /// Bridge indices in path order.
    bridges: Vec<usize>,
    circular: bool,
}

/// Anchors joined by bridges, as connected components (each sorted, ordered by first anchor).
/// Anchors without any bridge are left out.
fn components(n_seg: usize, br: &[Bridge]) -> Vec<Vec<usize>> {
    let mut parent: Vec<usize> = (0..n_seg).collect();
    fn root(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    let mut bridged = vec![false; n_seg];
    for b in br {
        let (x, y) = (b.walk[0].0, b.walk[b.walk.len() - 1].0);
        bridged[x] = true;
        bridged[y] = true;
        let (rx, ry) = (root(&mut parent, x), root(&mut parent, y));
        parent[rx] = ry;
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for a in (0..n_seg).filter(|&a| bridged[a]) {
        let r = root(&mut parent, a);
        groups.entry(r).or_default().push(a);
    }
    let mut out: Vec<Vec<usize>> = groups.into_values().collect();
    out.sort();
    out
}

/// Length-weighted median depth of some anchors (the depth of the base in the middle of them),
/// as `ovasm check evidence` takes it over the whole graph.
fn weighted_median_depth(g: &Graph, depth: impl Fn(usize) -> f64, segs: &[usize]) -> f64 {
    let mut v: Vec<(f64, usize)> = segs.iter().map(|&s| (depth(s), g.seqs[s].len())).collect();
    v.sort_by(|a, b| a.0.total_cmp(&b.0));
    let half = v.iter().map(|x| x.1).sum::<usize>() / 2;
    let mut acc = 0;
    v.iter()
        .find(|x| {
            acc += x.1;
            acc > half
        })
        .map_or(0.0, |x| x.0)
}

/// Scale evidence work with the number of candidate covers on bounded components.
/// A common fixed budget could discard *all* read-based rankings simply because a
/// previously omitted cover size was included. Large components and candidate sets
/// retain the original budget; the expanded budget has a hard upper bound.
fn candidate_work_budget(base: u64, anchors: usize, candidates: usize) -> u64 {
    if anchors <= 64 && candidates <= 100_000 {
        base.saturating_mul(candidates.div_ceil(20_000).clamp(1, 5) as u64)
    } else {
        base
    }
}

fn scoring_work_budget(anchors: usize, candidates: usize) -> u64 {
    candidate_work_budget(crate::read_model::WORK_BUDGET, anchors, candidates)
}

/// All covers of `anchors` (one component) by exactly `m` molecules that together visit every
/// required anchor once (either orientation); optional anchors may be passed through, each at most
/// once. A molecule ends where it closes back to its start or, failing that, as a path.
fn search(
    n_seg: usize,
    anchors: &[usize],
    optional: &[bool],
    br: &[Bridge],
    m: usize,
    max_steps: u64,
) -> (Vec<Vec<Found>>, bool) {
    search_phased(
        n_seg,
        anchors,
        optional,
        br,
        m,
        max_steps,
        &crate::phase::Constraints::default(),
    )
}

fn search_phased(
    n_seg: usize,
    anchors: &[usize],
    optional: &[bool],
    br: &[Bridge],
    m: usize,
    max_steps: u64,
    phases: &crate::phase::Constraints,
) -> (Vec<Vec<Found>>, bool) {
    let mut out_of: FxHashMap<Oriented, Vec<usize>> = FxHashMap::default();
    for (i, b) in br.iter().enumerate() {
        out_of.entry(b.walk[0]).or_default().push(i);
    }
    let (covers, truncated, _) =
        search_phased_indexed(n_seg, anchors, optional, br, m, max_steps, phases, &out_of, None);
    (covers, truncated)
}

fn search_phased_indexed(
    n_seg: usize,
    anchors: &[usize],
    optional: &[bool],
    br: &[Bridge],
    m: usize,
    max_steps: u64,
    phases: &crate::phase::Constraints,
    out_of: &FxHashMap<Oriented, Vec<usize>>,
    thread_ends: Option<&crate::read_model::ThreadEnds>,
) -> (Vec<Vec<Found>>, bool, bool) {
    search_phased_indexed_impl(n_seg, anchors, optional, br, m, max_steps, phases, out_of, thread_ends, true)
}

/// Conservative directed reachability on the anchor endpoints. If the component
/// cannot be represented exactly or exceeds its bounds, disable this shortcut.
fn component_reachability(
    anchors: &[usize], br: &[Bridge], out_of: &FxHashMap<Oriented, Vec<usize>>,
) -> Option<(FxHashMap<Oriented, usize>, Vec<Vec<bool>>)> {
    if anchors.is_empty() || anchors.len() > 256 { return None; }
    let mut index = FxHashMap::default();
    for &a in anchors { for side in [false, true] {
        let next = index.len();
        if index.insert((a, side), next).is_some() { return None; }
    }}
    let n = index.len();
    let mut adj = vec![Vec::new(); n];
    let mut edges = 0usize;
    for (&node, &i) in &index {
        for &j in out_of.get(&node).map(Vec::as_slice).unwrap_or(&[]) {
            let b = br.get(j)?;
            if b.walk.first().copied() != Some(node) { return None; }
            let end = *b.walk.last()?;
            let dst = *index.get(&end)?;
            adj[i].push(dst);
            edges = edges.checked_add(1)?;
            if edges > 65_536 { return None; }
        }
    }
    if n.checked_mul(n.checked_add(edges)?)? > 8_000_000 { return None; }
    let mut reach = vec![vec![false; n]; n];
    for src in 0..n {
        let mut stack = vec![src]; reach[src][src] = true;
        while let Some(x) = stack.pop() { for &y in &adj[x] {
            if !reach[src][y] { reach[src][y] = true; stack.push(y); }
        }}
    }
    Some((index, reach))
}

fn search_phased_indexed_impl(
    n_seg: usize, anchors: &[usize], optional: &[bool], br: &[Bridge], m: usize,
    max_steps: u64, phases: &crate::phase::Constraints,
    out_of: &FxHashMap<Oriented, Vec<usize>>,
    thread_ends: Option<&crate::read_model::ThreadEnds>, optimize: bool,
) -> (Vec<Vec<Found>>, bool, bool) {
    search_phased_indexed_bounded(n_seg, anchors, optional, br, m, max_steps, phases, out_of, thread_ends, optimize, None)
}

fn search_phased_indexed_bounded(
    n_seg: usize,
    anchors: &[usize],
    optional: &[bool],
    br: &[Bridge],
    m: usize,
    max_steps: u64,
    phases: &crate::phase::Constraints,
    out_of: &FxHashMap<Oriented, Vec<usize>>,
    thread_ends: Option<&crate::read_model::ThreadEnds>,
    optimize: bool,
    candidate_limit: Option<usize>,
) -> (Vec<Vec<Found>>, bool, bool) {
    struct Ctx<'a> {
        phases: &'a crate::phase::Constraints,
        out_of: &'a FxHashMap<Oriented, Vec<usize>>,
        br: &'a [Bridge],
        optional: &'a [bool],
        anchors: &'a [usize],
        want: usize,
        m: usize,
        max_steps: u64,
        candidate_limit: Option<usize>,
        topology_neutral: bool,
        optimize: bool,
        reach: Option<&'a (FxHashMap<Oriented, usize>, Vec<Vec<bool>>)>,
    }
    struct State {
        used: Vec<bool>,
        /// Required anchors visited so far.
        covered: usize,
        /// Finished molecules of the cover being built.
        done: Vec<Found>,
        /// Bridges of the molecule being built.
        path: Vec<usize>,
        /// Required anchor the current molecule must visit: the smallest one unused when it
        /// started, so the molecules of a cover are found in one order only.
        must: usize,
        steps: u64,
        truncated: bool,
        found: Vec<Vec<Found>>,
    }
    fn visits(c: &Ctx, start: Oriented, path: &[usize], a: usize) -> bool {
        start.0 == a || path.iter().any(|&i| c.br[i].walk.last().unwrap().0 == a)
    }
    fn dfs(c: &Ctx, s: &mut State, start: Oriented, at: Oriented) {
        s.steps += 1;
        if s.steps > c.max_steps {
            s.truncated = true;
            return;
        }
        // Every later molecule must consume a distinct unused required anchor.
        // Once too few remain, no continuation of this partial cover can succeed.
        if c.optimize && c.want.saturating_sub(s.covered) < c.m.saturating_sub(s.done.len() + 1) {
            return;
        }
        // Ignore used vertices and phase restrictions in this upper-bound graph:
        // only provably unreachable mandatory anchors are pruned.
        if !s.used[s.must] {
            if let Some((index, reach)) = c.reach {
                if let (Some(&from), Some(&p), Some(&q)) = (index.get(&at), index.get(&(s.must, false)), index.get(&(s.must, true))) {
                    if !reach[from][p] && !reach[from][q] { return; }
                }
            }
        }
        if s.covered < c.want {
            for &i in c.out_of.get(&at).map(Vec::as_slice).unwrap_or(&[]) {
                let next = *c.br[i].walk.last().unwrap();
                if s.used[next.0] {
                    continue;
                }
                // No node can block a walk when no phase constraints are active.
                if c.phases.active_nodes > 0 {
                    let mut walk = vec![start];
                    for &j in &s.path {
                        walk.extend_from_slice(&c.br[j].walk[1..]);
                    }
                    walk.extend_from_slice(&c.br[i].walk[1..]);
                    if !c.phases.allows_walk(&walk, false) {
                        continue;
                    }
                }
                s.used[next.0] = true;
                s.covered += usize::from(!c.optional[next.0]);
                s.path.push(i);
                dfs(c, s, start, next);
                s.path.pop();
                s.covered -= usize::from(!c.optional[next.0]);
                s.used[next.0] = false;
                if s.truncated {
                    return;
                }
            }
        }
        // End the molecule here.
        // `must` was unused when this molecule began; only this molecule can
        // have set its bit. Keep the original scan as an oracle for tests.
        let visited_must = if c.optimize { s.used[s.must] } else { visits(c, start, &s.path, s.must) };
        debug_assert_eq!(visited_must, visits(c, start, &s.path, s.must));
        if !visited_must { return; }
        let mut closings: Vec<usize> = c.out_of.get(&at).map(|v| {
            v.iter()
                .copied()
                .filter(|&i| *c.br[i].walk.last().unwrap() == start)
                .filter(|&i| {
                    if c.phases.active_nodes == 0 {
                        return true;
                    }
                    let mut walk = vec![start];
                    for &j in &s.path {
                        walk.extend_from_slice(&c.br[j].walk[1..]);
                    }
                    walk.extend_from_slice(&c.br[i].walk[1..]);
                    walk.pop();
                    c.phases.allows_walk(&walk, true)
                })
                .collect()
        }).unwrap_or_default();
        if !c.topology_neutral {
            closings = closings.into_iter().max_by_key(|&i| c.br[i].reads).into_iter().collect();
        }
        // A path ends only where it cannot go on: no bridge from either end reaches a required
        // anchor still unused. Cutting it anywhere else just multiplies covers.
        let open_to = |o: Oriented| {
            c.out_of.get(&o).is_some_and(|v| {
                v.iter().any(|&i| {
                    let a = c.br[i].walk.last().unwrap().0;
                    if s.used[a] || c.optional[a] {
                        return false;
                    }
                    if c.phases.active_nodes == 0 {
                        return true;
                    }
                    let mut walk = vec![start];
                    for &j in &s.path {
                        walk.extend_from_slice(&c.br[j].walk[1..]);
                    }
                    if o == at {
                        walk.extend_from_slice(&c.br[i].walk[1..]);
                    } else {
                        walk = rc_walk(&walk);
                        walk.extend_from_slice(&c.br[i].walk[1..]);
                    }
                    c.phases.allows_walk(&walk, false)
                })
            })
        };
        let can_end_open = !open_to(at) && !open_to(flip(start));
        if closings.is_empty() && !can_end_open {
            return;
        }
        let mut endings: Vec<Option<usize>> = closings.iter().copied().map(Some).collect();
        if endings.is_empty() || (c.topology_neutral && can_end_open) { endings.push(None); }
        // Every closing bridge competes, and an admissible open path also competes.
        // The previous maximum-count closing preselection could erase whole-read linkage.
        for closing in endings {
        if c.topology_neutral {
            s.steps += 1;
            if s.steps > c.max_steps { s.truncated = true; return; }
        }
        let mut bridges = s.path.clone();
        bridges.extend(closing);
        let mol = Found {
            start,
            bridges,
            circular: closing.is_some(),
        };
        if s.covered == c.want {
            if s.done.len() + 1 == c.m {
                let mut cover = s.done.clone();
                cover.push(mol);
                if c.candidate_limit.is_some_and(|limit| s.found.len() >= limit) {
                    eprintln!("[ovasm] warning: adaptive search reached its candidate limit");
                    s.truncated = true;
                    return;
                }
                s.found.push(cover);
            }
            continue;
        }
        if s.done.len() + 1 == c.m {
            continue;
        }
        // Start the next molecule; it must take the smallest required anchor still unused.
        let Some(&must) = c.anchors.iter().find(|&&a| !s.used[a] && !c.optional[a]) else {
            return;
        };
        let (saved_path, saved_must) = (std::mem::take(&mut s.path), s.must);
        s.done.push(mol);
        s.must = must;
        let free: Vec<usize> = c.anchors.iter().copied().filter(|&a| !s.used[a]).collect();
        for a in free {
            for st in [(a, true), (a, false)] {
                s.used[a] = true;
                s.covered += usize::from(!c.optional[a]);
                dfs(c, s, st, st);
                s.covered -= usize::from(!c.optional[a]);
                s.used[a] = false;
                if s.truncated {
                    break;
                }
            }
            if s.truncated {
                break;
            }
        }
        s.done.pop();
        s.path = saved_path;
        s.must = saved_must;
        if s.truncated { return; }
        }
    }

    fn charge(budget: &mut u64, work: usize) -> Option<()> {
        *budget = budget.checked_sub(work as u64)?;
        Some(())
    }

    /// Match a record end against existing anchor bridges. Only repeats
    /// occur inside a bridge, so two different bridges cannot both be a
    /// prefix of the same sequence of oriented segments.
    fn terminal_bridges(
        c: &Ctx, start: Oriented, mut tail: &[Oriented], occupied: &FxHashSet<usize>, budget: &mut u64,
    ) -> Option<Vec<usize>> {
        let mut at = start;
        let mut added = FxHashSet::default();
        let mut bridges = Vec::new();
        while !tail.is_empty() {
            let mut matched = None;
            for &id in c.out_of.get(&at).into_iter().flatten() {
                let walk = &c.br[id].walk;
                charge(budget, 1)?;
                if walk.len() - 1 > tail.len() {
                    continue;
                }
                charge(budget, walk.len() - 1)?;
                if tail.starts_with(&walk[1..]) {
                    matched = Some(id);
                    break;
                }
            }
            let Some(id) = matched else { return Some(Vec::new()); };
            let walk = &c.br[id].walk;
            at = *walk.last().unwrap();
            if !c.optional[at.0] || occupied.contains(&at.0) || !added.insert(at.0) {
                return Some(Vec::new());
            }
            bridges.push(id);
            tail = &tail[walk.len() - 1..];
        }
        Some(bridges)
    }

    fn extend_terminal_covers(
        c: &Ctx, threads: &crate::read_model::ThreadEnds, covers: &mut Vec<Vec<Found>>, mut budget: u64,
    ) -> Option<()> {
        // Complete the original search first. The remaining budget can only add
        // open alternatives and cannot erase an original cover.
        for i in 0..covers.len() {
            charge(&mut budget, 1)?;
            let cover = &covers[i];
            let last = cover.last().unwrap();
            if last.circular {
                continue;
            }
            let (at, previous) = last.bridges.last().map_or((last.start, None), |&j| {
                let walk = &c.br[j].walk;
                (walk[walk.len() - 1], Some(walk[walk.len() - 2]))
            });
            if c.optional[at.0] || !threads.has_anchor(at, previous) {
                continue;
            }
            charge(&mut budget, 1 + last.bridges.len())?;
            let length = 1 + last.bridges.iter().map(|&j| c.br[j].walk.len() - 1).sum::<usize>();
            charge(&mut budget, length)?;
            let mut walk = candidate_walk(c.br, last);
            let tails = threads.tails(&walk, &mut budget)?;
            if tails.is_empty() {
                continue;
            }
            charge(&mut budget, cover.len())?;
            let cost = cover.iter().map(|f| f.bridges.len() + 1).sum::<usize>();
            charge(&mut budget, cost)?;
            let occupied = cover.iter().flat_map(|f| {
                std::iter::once(f.start.0).chain(f.bridges.iter().map(|&j| c.br[j].walk.last().unwrap().0))
            }).collect();
            let cover = cover.clone();
            for tail in tails {
                let bridges = terminal_bridges(c, at, tail, &occupied, &mut budget)?;
                if bridges.is_empty() {
                    continue;
                }
                if c.phases.active_nodes > 0 {
                    charge(&mut budget, length + tail.len())?;
                    walk.extend_from_slice(tail);
                    let allowed = c.phases.allows_walk(&walk, false);
                    walk.truncate(length);
                    if !allowed {
                        continue;
                    }
                }
                charge(&mut budget, cost + bridges.len())?;
                let mut extended = cover.clone();
                extended.last_mut().unwrap().bridges.extend(bridges);
                if c.candidate_limit.is_some_and(|limit| covers.len() >= limit) {
                    eprintln!("[ovasm] warning: adaptive terminal extension reached its candidate limit");
                    return None;
                }
                covers.push(extended);
            }
        }
        Some(())
    }

    let reach = optimize.then(|| component_reachability(anchors, br, out_of)).flatten();
    let ctx = Ctx {
        phases,
        out_of,
        br,
        optional,
        anchors,
        want: anchors.iter().filter(|&&a| !optional[a]).count(),
        m,
        max_steps,
        candidate_limit,
        topology_neutral: thread_ends.is_some(),
        optimize,
        reach: reach.as_ref(),
    };
    let mut st = State {
        used: vec![false; n_seg],
        covered: 0,
        done: Vec::new(),
        path: Vec::new(),
        must: anchors
            .iter()
            .copied()
            .find(|&a| !optional[a])
            .unwrap_or(anchors[0]),
        steps: 0,
        truncated: false,
        found: Vec::new(),
    };
    for &a in anchors {
        for start in [(a, true), (a, false)] {
            st.used[a] = true;
            st.covered = usize::from(!optional[a]);
            dfs(&ctx, &mut st, start, start);
            st.used[a] = false;
            if st.truncated {
                break;
            }
        }
        if st.truncated {
            break;
        }
    }
    let terminal_truncated = !st.truncated && thread_ends.is_some_and(|threads| {
        !threads.is_empty() && extend_terminal_covers(&ctx, threads, &mut st.found, max_steps - st.steps).is_none()
    });
    (st.found, st.truncated, terminal_truncated)
}

/// Longest tail `cut_tail` takes off a circle, as a fraction of the circle's anchor sequence.
const TAIL_MAX_FRACTION: f64 = 0.1;

/// A linear molecule whose end leads back into its own interior is a circle with a tail: the
/// fewest-molecules rule opens the circle to thread a side branch in. Luzula sylvatica's
/// mitochondrion came out as one linear molecule of 649,000 bp, u4 (15,644 bp) leading into a
/// path of 633,356 bp whose last anchor is bridged back to its first by 30 reads; the published
/// assembly is that circle and a 17 kb linear record. Such a molecule is split where the
/// bridge from its last anchor re-enters it (or, at the other end, where a bridge enters its
/// first anchor from inside), into the circle and the tail; the bridge between the two is
/// left out. Only a bridge reads support in full (neither weak nor topological) closes the
/// circle, and only a tail of at most `TAIL_MAX_FRACTION` of the circle's anchor sequence is
/// cut off: a longer one is a part of the molecule in its own right, and the path stays whole
/// (Nipponbare, whose 376 kb path would fall into 190 + 186 kb from ONT reads and 226 + 150 kb
/// from short reads). Of several possible cuts the largest circle wins.
///
/// Returns the circle, the tail and the index of the bridge left out.
fn cut_tail(g: &Graph, br: &[Bridge], f: &Found) -> Option<(Found, Found, usize)> {
    cut_tail_phased(g, br, f, &crate::phase::Constraints::default())
}

fn cut_tail_phased(
    g: &Graph,
    br: &[Bridge],
    f: &Found,
    phases: &crate::phase::Constraints,
) -> Option<(Found, Found, usize)> {
    if f.circular || f.bridges.is_empty() {
        return None;
    }
    // anchors in path order
    let mut at = vec![f.start];
    at.extend(f.bridges.iter().map(|&i| *br[i].walk.last().unwrap()));
    let n = at.len() - 1;
    let bp = |r: std::ops::Range<usize>| -> usize { r.map(|i| g.seqs[at[i].0].len()).sum() };
    let closing = |from: Oriented, to: Oriented| {
        (0..br.len())
            .filter(|&i| {
                let b = &br[i];
                !b.weak && !b.topological && b.walk[0] == from && *b.walk.last().unwrap() == to
            })
            .max_by_key(|&i| br[i].reads)
    };
    // (circle bp, closing reads, circle, tail, bridge left out)
    let mut best: Option<(usize, u64, Found, Found, usize)> = None;
    let mut offer =
        |circle_bp: usize, tail_bp: usize, c: usize, circle: Found, tail: Found, cut| {
            if !phases.allows_walk(&candidate_walk(br, &circle), true)
                || !phases.allows_walk(&candidate_walk(br, &tail), false)
            {
                return;
            }
            if tail_bp as f64 <= TAIL_MAX_FRACTION * circle_bp as f64
                && best
                    .as_ref()
                    .map_or(true, |b| (circle_bp, br[c].reads) > (b.0, b.1))
            {
                best = Some((circle_bp, br[c].reads, circle, tail, cut));
            }
        };
    for j in 1..=n {
        // the last anchor leads back to anchor j: anchors j..=n close, 0..j are the tail
        if let Some(c) = closing(at[n], at[j]) {
            let mut bridges = f.bridges[j..].to_vec();
            bridges.push(c);
            let circle = Found {
                start: at[j],
                bridges,
                circular: true,
            };
            let tail = Found {
                start: at[0],
                bridges: f.bridges[..j - 1].to_vec(),
                circular: false,
            };
            offer(bp(j..n + 1), bp(0..j), c, circle, tail, f.bridges[j - 1]);
        }
    }
    for j in 0..n {
        // anchor j leads back to the first: anchors 0..=j close, the rest are the tail
        if let Some(c) = closing(at[j], at[0]) {
            let mut bridges = f.bridges[..j].to_vec();
            bridges.push(c);
            let circle = Found {
                start: at[0],
                bridges,
                circular: true,
            };
            let tail = Found {
                start: at[j + 1],
                bridges: f.bridges[j + 1..].to_vec(),
                circular: false,
            };
            offer(
                bp(0..j + 1),
                bp(j + 1..n + 1),
                c,
                circle,
                tail,
                f.bridges[j],
            );
        }
    }
    best.map(|(_, _, circle, tail, cut)| (circle, tail, cut))
}

/// Paths for a component that no cover of up to `max_molecules` molecules solves: the longest
/// walks its read bridges make, by greedy choice. Without them such a component gave nothing
/// at all (Bromus sterilis, frozen seeds: 32 anchors, five of them bridged at one end only, no
/// cover of 3 or 4 molecules, and 337 kb of anchors that are all in the published
/// mitochondrion went unreported). Bridges are taken best first: between required anchors
/// before those reaching optional ones, fully read before weak, then by their share of the
/// reads at their two ends, then by reads. A bridge is taken when the anchor end it leaves and
/// the one it enters are both free and it joins two different paths; a path whose last anchor
/// is bridged back to its first is then closed. Topological bridges, which no read took, are
/// not used: a partial path claims only what reads show. Optional anchors stay only inside a
/// path. Every required anchor of the component is on exactly one path, alone when none of
/// its bridges is taken.
///
/// These are not a cover: they say which anchors the reads join, not how the component
/// unrings, and a repeat is laid down wherever a taken bridge crosses it.
fn partial_paths(
    n_seg: usize,
    comp: &[usize],
    optional: &[bool],
    br: &[Bridge],
    frac: &[f64],
) -> Vec<Found> {
    partial_paths_phased(
        n_seg,
        comp,
        optional,
        br,
        frac,
        &crate::phase::Constraints::default(),
    )
}

fn partial_paths_phased(
    n_seg: usize,
    comp: &[usize],
    optional: &[bool],
    br: &[Bridge],
    frac: &[f64],
    phases: &crate::phase::Constraints,
) -> Vec<Found> {
    fn root(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    let pos: FxHashMap<&[Oriented], usize> = br
        .iter()
        .enumerate()
        .map(|(i, b)| (b.walk.as_slice(), i))
        .collect();
    let mut inside = vec![false; n_seg];
    for &a in comp {
        inside[a] = true;
    }
    let ends = |i: usize| (br[i].walk[0], *br[i].walk.last().unwrap());
    let mut order: Vec<usize> = (0..br.len())
        .filter(|&i| !br[i].topological && inside[br[i].walk[0].0])
        .collect();
    order.sort_by(|&x, &y| {
        let class = |i: usize| {
            let (a, b) = ends(i);
            (optional[a.0] || optional[b.0], br[i].weak)
        };
        class(x)
            .cmp(&class(y))
            .then(frac[y].total_cmp(&frac[x]))
            .then(br[y].reads.cmp(&br[x].reads))
            .then_with(|| br[x].walk.cmp(&br[y].walk))
    });
    // The bridge by which a path leaves each oriented anchor. Entering b uses the end that
    // flip(b) leaves by, so a taken bridge is entered at both of its ends (reversed at the far
    // one).
    let mut exit: FxHashMap<Oriented, usize> = FxHashMap::default();
    let mut parent: Vec<usize> = (0..n_seg).collect();
    for i in order {
        let (a, b) = ends(i);
        if exit.contains_key(&a) || exit.contains_key(&flip(b)) {
            continue;
        }
        let (ra, rb) = (root(&mut parent, a.0), root(&mut parent, b.0));
        if ra == rb {
            // would close a path on itself; closing comes last, at full length
            continue;
        }
        let follow = |start: Oriented| {
            let mut walk = vec![start];
            while let Some(&j) = exit.get(walk.last().unwrap()) {
                walk.extend_from_slice(&br[j].walk[1..]);
            }
            walk
        };
        let mut joined = rc_walk(&follow(flip(a)));
        joined.extend_from_slice(&br[i].walk[1..]);
        joined.extend_from_slice(&follow(b)[1..]);
        if !phases.allows_walk(&joined, false) {
            continue;
        }
        parent[ra] = rb;
        exit.insert(a, i);
        exit.insert(flip(b), pos[rc_walk(&br[i].walk).as_slice()]);
    }
    let mut closing_candidates: FxHashMap<(Oriented, Oriented), Vec<usize>> = FxHashMap::default();
    for i in 0..br.len() {
        if !br[i].topological {
            closing_candidates.entry(ends(i)).or_default().push(i);
        }
    }
    let closing = |from: Oriented, to: Oriented| {
        closing_candidates.get(&(from, to)).into_iter().flatten().copied()
            .filter(|&i| {
                let mut walk = vec![to];
                while *walk.last().unwrap() != from {
                    let Some(&j) = exit.get(walk.last().unwrap()) else {
                        return false;
                    };
                    walk.extend_from_slice(&br[j].walk[1..]);
                }
                walk.extend_from_slice(&br[i].walk[1..]);
                walk.pop();
                phases.allows_walk(&walk, true)
            })
            .max_by_key(|&i| (!br[i].weak, br[i].reads))
    };
    let mut out = Vec::new();
    for &a in comp {
        // a path starts at an anchor with a free end (nothing was closed, so each path has two)
        let start = if !exit.contains_key(&(a, false)) {
            (a, true)
        } else if !exit.contains_key(&(a, true)) {
            (a, false)
        } else {
            continue;
        };
        let mut at = vec![start];
        let mut bridges = Vec::new();
        while let Some(&i) = exit.get(&at[at.len() - 1]) {
            bridges.push(i);
            at.push(ends(i).1);
        }
        // each path is met from both of its ends: keep it from the lower one
        if at.len() > 1 && at[at.len() - 1].0 < a {
            continue;
        }
        let required: Vec<usize> = (0..at.len()).filter(|&i| !optional[at[i].0]).collect();
        let (Some(&lo), Some(&hi)) = (required.first(), required.last()) else {
            continue;
        };
        if let Some(c) = closing(at[at.len() - 1], at[0]) {
            bridges.push(c);
            out.push(Found {
                start,
                bridges,
                circular: true,
            });
            continue;
        }
        let mut bridges = bridges[lo..hi].to_vec();
        let closed = closing(at[hi], at[lo]);
        bridges.extend(closed);
        out.push(Found {
            start: at[lo],
            bridges,
            circular: closed.is_some(),
        });
    }
    out
}

/// Segment walk of a candidate (for a circular one, without repeating the first anchor).
fn candidate_walk(br: &[Bridge], f: &Found) -> Vec<Oriented> {
    let mut walk = vec![f.start];
    for &i in &f.bridges {
        walk.extend_from_slice(&br[i].walk[1..]);
    }
    if f.circular {
        walk.pop();
    }
    walk
}

/// The bridges `run` joins the anchors of a graph by, and which segments are repeats: what a
/// caller needs to weigh the bridges from outside the reads before `run_with_prior`.
pub(crate) fn open_bridges(gfa: &Path, evidence: &Path, min_reads: u64) -> Result<(Graph, Vec<bool>, Vec<BridgeView>)> {
    let g = parse_gfa(gfa)?;
    let report: Value = serde_json::from_str(
        &fs::read_to_string(evidence).with_context(|| format!("cannot read {}", evidence.display()))?,
    )?;
    let repeat_names: BTreeSet<&str> = report["repeat_segments"]
        .as_array()
        .context("evidence has no repeat_segments")?
        .iter()
        .filter_map(Value::as_str)
        .collect();
    let repeat: Vec<bool> = g.names.iter().map(|n| repeat_names.contains(n.as_str())).collect();
    let (mut br, _) = bridges(&g, &report, &repeat, min_reads)?;
    let spanned = spanned_walks(&g, &report, &repeat, min_reads)?;
    merge_spanned(&mut br, &spanned);
    let views = br
        .into_iter()
        .map(|b| BridgeView { walk: b.walk, reads: b.reads, topological: b.topological, weak: b.weak })
        .collect();
    Ok((g, repeat, views))
}

/// The closed cover of a component (every anchor once, every anchor end on one bridge) that a
/// prior takes among those the reads rank alike: `pan_posterior::best_matching` over the
/// component's bridges, each with its share of the reads and the prior's probability at both of
/// its ends (an end the prior does not know: its bridges alike). An anchor end no bridge leaves
/// is where the graph stops: it ends a linear molecule, as in every cover, and an optional
/// anchor may be left out. The `PriorCover` comes back without its report. None for a
/// component whose bridged ends cannot be joined but for two of a repeat's, or with more
/// matchings than the budget.
fn prior_cover(
    comp: &[usize],
    optional: &[bool],
    br: &[Bridge],
    frac: &[f64],
    expected: &[Option<usize>],
    prior: &BridgePrior,
    model: Option<&crate::read_model::ReadModel>,
) -> Option<(Vec<Found>, PriorCover)> {
    let local: FxHashMap<usize, usize> = comp.iter().enumerate().map(|(k, &a)| (a, k)).collect();
    // an anchor end, as the oriented anchor a walk leaves it by: anchor k has ends 2k and 2k + 1
    let end = |o: Oriented| local.get(&o.0).map(|&k| 2 * k + usize::from(!o.1));
    let pos: FxHashMap<&[Oriented], usize> = br.iter().enumerate().map(|(i, b)| (b.walk.as_slice(), i)).collect();
    let mut leaving: Vec<Vec<usize>> = vec![Vec::new(); 2 * comp.len()];
    for (i, b) in br.iter().enumerate() {
        if let Some(e) = end(b.walk[0]) {
            leaving[e].push(i);
        }
    }
    let probability = |i: usize| -> f64 {
        prior.probability.get(&br[i].walk).copied().unwrap_or_else(|| {
            let e = end(br[i].walk[0]).expect("a bridge of the component");
            1.0 / leaving[e].len().max(1) as f64
        })
    };
    let mut edges = Vec::new();
    // per edge: the bridge leaving its first end and the same bridge leaving its other end
    let mut directions: Vec<(usize, usize)> = Vec::new();
    for (i, b) in br.iter().enumerate() {
        let (Some(a), Some(z)) = (end(b.walk[0]), end(flip(b.walk[b.walk.len() - 1]))) else {
            continue;
        };
        let j = *pos.get(rc_walk(&b.walk).as_slice())?;
        // a walk back into the end it left is no junction of two ends; each bridge once
        if a == z || j <= i {
            continue;
        }
        edges.push(crate::pan_posterior::Edge {
            a,
            b: z,
            read_ln: frac[i].max(1e-9).ln(),
            incomplete: b.weak || b.topological,
            prior: probability(i) + probability(j),
            passes: b.walk[1..b.walk.len() - 1].iter().map(|o| o.0).collect(),
        });
        directions.push((i, j));
    }
    // An optional (shallow) anchor may be left out, as in every cover: its two ends joined to
    // each other, with no bridge. Leaving it out costs no read share, so it is left out unless
    // an end beside it would stay open, as the search does.
    const LEFT_OUT: usize = usize::MAX - 1;
    for (k, _) in comp.iter().enumerate().filter(|&(_, &a)| optional[a]) {
        edges.push(crate::pan_posterior::Edge { a: 2 * k, b: 2 * k + 1, read_ln: 0.0, incomplete: false, prior: 0.0, passes: Vec::new() });
        directions.push((LEFT_OUT, LEFT_OUT));
    }
    // The target's own single records (one read through three segments or more) come before
    // the prior, as they come first where `cover_order` has them (`ReadFit::order`). The anchor
    // ends a record crosses are weighed together, and a matching explains the record when its
    // bridges, joined at the anchors they share, hold the record's walk.
    let threads: Vec<(&[Oriented], u64, bool)> = model.map(|m| m.threads().collect()).unwrap_or_default();
    let coupled: Vec<Vec<usize>> = threads
        .iter()
        .map(|(p, _, _)| {
            let mut crossed = Vec::new();
            for (i, &o) in p.iter().enumerate() {
                if i + 1 < p.len() {
                    crossed.extend(end(o));
                }
                if i > 0 {
                    crossed.extend(end(flip(o)));
                }
            }
            crossed.sort_unstable();
            crossed.dedup();
            crossed
        })
        .filter(|crossed| crossed.len() > 1)
        .collect();
    let records = |taken: &[usize]| -> (u64, u64) {
        if threads.is_empty() {
            return (0, 0);
        }
        // the bridges taken, each in both directions, by the oriented anchor they leave
        let mut leaving_anchor: FxHashMap<Oriented, usize> = FxHashMap::default();
        let mut taken_bridges: Vec<usize> = Vec::new();
        for &e in taken {
            let (i, j) = directions[e];
            if i == LEFT_OUT {
                continue;
            }
            for b in [i, j] {
                leaving_anchor.insert(br[b].walk[0], b);
                taken_bridges.push(b);
            }
        }
        let holds = |p: &[Oriented]| -> bool {
            taken_bridges.iter().any(|&start| {
                (0..br[start].walk.len()).any(|from| {
                    let (mut w, mut k) = (&br[start].walk, from);
                    for (t, &node) in p.iter().enumerate() {
                        if w[k] != node {
                            return false;
                        }
                        if t + 1 == p.len() {
                            return true;
                        }
                        if k + 1 < w.len() {
                            k += 1;
                        } else {
                            // at the bridge's last anchor: on along the bridge taken at its other end
                            match leaving_anchor.get(&w[k]) {
                                Some(&next) => {
                                    w = &br[next].walk;
                                    k = 1;
                                }
                                None => return false,
                            }
                        }
                    }
                    false
                })
            })
        };
        let (mut joint, mut all) = (0u64, 0u64);
        for &(p, n, through_two) in &threads {
            if holds(p) {
                all += n;
                if through_two {
                    joint += n;
                }
            }
        }
        (joint, all)
    };
    let m = crate::pan_posterior::best_matching(2 * comp.len(), &edges, &|r| expected[r], COPY_PENALTY_LN, &coupled, &records, 2_000_000)?;
    // the bridge leaving each end; none where the graph stops (the end of a linear molecule)
    let leave: Vec<usize> = m
        .edge_at
        .iter()
        .enumerate()
        .map(|(x, &e)| match e {
            usize::MAX => usize::MAX,
            e if edges[e].a == x => directions[e].0,
            e => directions[e].1,
        })
        .collect();
    // the molecules: leave an anchor by one end, enter the next, leave it by its other end
    let oriented = |x: usize| (comp[x / 2], x % 2 == 0);
    let mut on_molecule: Vec<bool> = (0..comp.len()).map(|k| leave[2 * k] == LEFT_OUT).collect();
    let left_out = on_molecule.iter().filter(|&&x| x).count();
    let mut cover = Vec::new();
    // linear ones first, each from an end the graph stops at to the next such end
    for stop in 0..2 * comp.len() {
        if leave[stop] != usize::MAX || on_molecule[stop / 2] {
            continue;
        }
        let (mut x, mut taken) = (stop ^ 1, Vec::new());
        on_molecule[x / 2] = true;
        while leave[x] != usize::MAX {
            taken.push(leave[x]);
            let entered = *br[leave[x]].walk.last().expect("a bridge has anchors");
            x = end(entered).expect("a bridge of the component");
            on_molecule[x / 2] = true;
        }
        cover.push(Found { start: oriented(stop ^ 1), bridges: taken, circular: false });
    }
    for k in 0..comp.len() {
        if on_molecule[k] {
            continue;
        }
        let (mut x, mut taken) = (2 * k, Vec::new());
        loop {
            on_molecule[x / 2] = true;
            taken.push(leave[x]);
            let entered = *br[leave[x]].walk.last().expect("a bridge has anchors");
            x = end(entered).expect("a bridge of the component");
            if x == 2 * k {
                break;
            }
        }
        cover.push(Found { start: (comp[k], true), bridges: taken, circular: true });
    }
    debug_assert_eq!(m.molecules, cover.len() + left_out);
    debug_assert_eq!(m.uses_incomplete, cover.iter().flat_map(|f| &f.bridges).any(|&i| br[i].weak || br[i].topological));
    if cover.is_empty() {
        return None;
    }
    let open: Vec<&Vec<usize>> = leaving.iter().filter(|l| l.len() > 1).collect();
    let decided = open
        .iter()
        .filter(|l| {
            let p: Vec<f64> = l.iter().map(|&i| probability(i)).collect();
            p.iter().copied().fold(f64::NEG_INFINITY, f64::max) - p.iter().copied().fold(f64::INFINITY, f64::min) > 1e-9
        })
        .count();
    let stats = PriorCover {
        cover: CoverReport { score: 0.0, molecules: Vec::new(), weak_bridges: 0, topological_bridges: 0, copy_deviation: 0, read_support: None },
        ends: leave.iter().filter(|&&i| i < LEFT_OUT).count(),
        open_ends: open.len(),
        decided_ends: decided,
        expected_correct_ends: m.expected,
        equally_good: m.equally_good,
        matchings_examined: m.examined,
        joint_records_explained: m.records.0,
        thread_records_explained: m.records.1,
    };
    Some((cover, stats))
}

/// One key for a path and its reverse complement, and for every rotation of a circle.
fn canonical_key(walk: &[Oriented], circular: bool) -> Vec<Oriented> {
    let rc = rc_walk(walk);
    let mut variants = vec![walk.to_vec(), rc.clone()];
    if circular {
        for r in 1..walk.len() {
            variants.push([&walk[r..], &walk[..r]].concat());
            variants.push([&rc[r..], &rc[..r]].concat());
        }
    }
    variants.into_iter().min().unwrap()
}

fn link_overlap(g: &Graph, a: Oriented, b: Oriented) -> Option<usize> {
    g.links
        .iter()
        .zip(&g.overlaps)
        .find(|(&(x, y), _)| (x, y) == (a, b) || (flip(y), flip(x)) == (a, b))
        .map(|(_, &ov)| ov)
}

fn overlap_index(g: &Graph) -> FxHashMap<(Oriented, Oriented), usize> {
    let mut overlaps = FxHashMap::default();
    for (&(a, b), &ov) in g.links.iter().zip(&g.overlaps) {
        // Match link_overlap's first-edge semantics, including reverse orientation.
        overlaps.entry((a, b)).or_insert(ov);
        overlaps.entry((flip(b), flip(a))).or_insert(ov);
    }
    overlaps
}

fn spelled_length(
    g: &Graph,
    walk: &[Oriented],
    circular: bool,
    overlaps: &FxHashMap<(Oriented, Oriented), usize>,
) -> Result<usize> {
    let first = *walk.first().context("empty molecule walk")?;
    let mut length = g.seqs[first.0].len();
    for w in walk.windows(2) {
        let ov = *overlaps.get(&(w[0], w[1])).context("missing link in molecule walk")?;
        let added = g.seqs[w[1].0].len().checked_sub(ov).context("overlap exceeds segment length")?;
        length = length.checked_add(added).context("molecule length overflow")?;
    }
    if circular {
        let ov = *overlaps.get(&(*walk.last().unwrap(), first)).context("missing closing link")?;
        length = length.checked_sub(ov).context("closing overlap exceeds molecule length")?;
    }
    Ok(length)
}

/// Sequence of a walk, each link's overlap written once. A circular walk also drops the overlap
/// of the closing link from its end.
pub(crate) fn spell(g: &Graph, walk: &[Oriented], circular: bool) -> Result<Vec<u8>> {
    let mut seq = g.oriented_seq(walk[0]);
    for w in walk.windows(2) {
        let ov = link_overlap(g, w[0], w[1]).with_context(|| {
            format!("no link {} -> {} in the GFA", g.label(w[0]), g.label(w[1]))
        })?;
        seq.extend_from_slice(&g.oriented_seq(w[1])[ov..]);
    }
    if circular {
        let (last, first) = (walk[walk.len() - 1], walk[0]);
        let ov = link_overlap(g, last, first)
            .with_context(|| format!("no closing link {} -> {}", g.label(last), g.label(first)))?;
        seq.truncate(seq.len() - ov);
    }
    Ok(seq)
}

#[derive(Clone, Serialize)]
pub struct Linearization {
    pub path: String,
    pub circular: bool,
    pub length: usize,
    /// Reads of each bridge in path order.
    pub bridge_reads: Vec<u64>,
    /// Each bridge's share of the reads leaving its anchor, in path order.
    pub exit_fractions: Vec<f64>,
    /// Log-likelihood of the path (sum of ln exit fractions).
    pub score: f64,
}

/// A set of molecules covering one component.
#[derive(Clone, Serialize)]
pub struct CoverReport {
    /// Sum of the molecules' scores.
    pub score: f64,
    /// Longest first.
    pub molecules: Vec<Linearization>,
    /// Bridges taken at otherwise dead anchor ends that the cover passes (see `bridges`).
    pub weak_bridges: usize,
    /// Junctions supplied by graph topology because no read phases them. Their synthetic
    /// weights describe possible paths, not observations of a chromosome configuration.
    pub topological_bridges: usize,
    /// Repeat copies laid down more or fewer times than the repeats' depths say (see
    /// `copy_deviation`); each costs `COPY_PENALTY_LN` in the ranking.
    pub copy_deviation: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read_support: Option<crate::read_model::ReadFit>,
}

impl CoverReport {
    fn linear(&self) -> usize {
        self.molecules.iter().filter(|m| !m.circular).count()
    }

    /// Score the ranking uses: the likelihood less the copy-number penalty.
    pub(crate) fn ranked_score(&self) -> f64 {
        self.score - COPY_PENALTY_LN * self.copy_deviation as f64
    }

    fn incomplete_junction_evidence(&self) -> bool {
        self.weak_bridges > 0 || self.topological_bridges > 0
    }
}

/// Repeat copies a cover lays down beyond or short of what the depths say: per repeat, its
/// passes in the cover against max(1, depth / median depth rounded). Read pairings decide which
/// copy goes where; depth says how many there are, which is what the reads cannot tell when few
/// of them span a repeat (low depth): a 2x repeat passed once leaves a copy out, a 1x segment
/// passed twice duplicates it.
fn copy_deviation(expected: &[Option<usize>], br: &[Bridge], cover: &[Found]) -> usize {
    let mut passes: FxHashMap<usize, usize> = FxHashMap::default();
    for f in cover {
        for &i in &f.bridges {
            let w = &br[i].walk;
            for o in &w[1..w.len() - 1] {
                *passes.entry(o.0).or_default() += 1;
            }
        }
    }
    let mut dev = 0;
    for (seg, &want) in expected.iter().enumerate() {
        if let Some(want) = want {
            let got = passes.get(&seg).copied().unwrap_or(0);
            // a repeat this cover never reaches belongs to another component or is not placed
            if got > 0 {
                dev += got.abs_diff(want);
            }
        }
    }
    dev
}

/// One key for a link and its reverse complement.
fn link_key(a: Oriented, b: Oriented) -> (Oriented, Oriented) {
    (a, b).min((flip(b), flip(a)))
}

/// Least times a molecule passes a link, from single reads (evidence `link_copies_in_reads`):
/// the most crossings at least `min_reads` reads show, each read alone. A read that crosses a
/// tandem array's unit junction n times holds n copies of the unit, so its molecule does too,
/// whatever the depth says; `min_reads` keeps one chimeric or misjoined read from setting it.
/// With it, the median read distance between consecutive crossings (the unit's length).
fn read_copies(
    g: &Graph,
    report: &Value,
    min_reads: u64,
) -> Result<FxHashMap<(Oriented, Oriented), (usize, usize)>> {
    let index: FxHashMap<&str, usize> =
        g.names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let mut out = FxHashMap::default();
    // older evidence reports have no such field: no read-measured copies
    let Some(rows) = report["link_copies_in_reads"].as_array() else {
        return Ok(out);
    };
    for row in rows {
        let spacing = row["median_spacing"].as_i64().unwrap_or(0).max(0) as usize;
        let a = parse_label(&index, row["from"].as_str().context("link copies without from")?)?;
        let b = parse_label(&index, row["to"].as_str().context("link copies without to")?)?;
        let mut by: Vec<(u64, u64)> = row["reads_by_crossings"]
            .as_array()
            .context("link copies without reads_by_crossings")?
            .iter()
            .filter_map(|x| Some((x[0].as_u64()?, x[1].as_u64()?)))
            .collect();
        by.sort_unstable_by(|x, y| y.cmp(x));
        let mut at_least = 0;
        for (n, reads) in by {
            at_least += reads;
            if at_least >= min_reads.max(1) {
                out.insert(link_key(a, b), (n as usize, spacing));
                break;
            }
        }
    }
    Ok(out)
}

/// Upper bound on the unit copies `lay_read_copies` adds to one molecule.
const MAX_READ_COPIES_ADDED: usize = 1_000;
/// Longest tandem unit (in segments) `lay_read_copies` looks for in the graph.
const MAX_UNIT_SEGMENTS: usize = 16;
/// How far (fraction, at least 50 bp) a unit's length may differ from the read distance between
/// consecutive crossings of its junction, as `--spacing-tol` for read walks.
const UNIT_SPACING_TOL: f64 = 0.05;

/// Shortest walk in the graph from `from` to `to` (both included; `[from]` when they are the
/// same), of at most `MAX_UNIT_SEGMENTS` segments.
fn shortest_return(
    next: &FxHashMap<Oriented, Vec<Oriented>>,
    from: Oriented,
    to: Oriented,
) -> Option<Vec<Oriented>> {
    if from == to {
        return Some(vec![from]);
    }
    let mut prev: FxHashMap<Oriented, Oriented> = FxHashMap::default();
    let mut layer = vec![from];
    for _ in 1..MAX_UNIT_SEGMENTS {
        let mut further = Vec::new();
        for &o in &layer {
            for &nx in next.get(&o).map(Vec::as_slice).unwrap_or(&[]) {
                if nx == from || prev.contains_key(&nx) {
                    continue;
                }
                prev.insert(nx, o);
                if nx == to {
                    let mut path = vec![to];
                    while *path.last().unwrap() != from {
                        path.push(prev[path.last().unwrap()]);
                    }
                    path.reverse();
                    return Some(path);
                }
                further.push(nx);
            }
        }
        if further.is_empty() {
            break;
        }
        layer = further;
    }
    None
}

/// A molecule with every link passed at least as often as single reads cross it (`read_copies`).
/// Where the walk passes a link x -> y too seldom, the shortest return y .. x in the graph (the
/// tandem unit; `[x]` for a self-loop) is laid between x and y, once per missing pass, when the
/// unit adds as many bases as the reads' crossings lie apart (`UNIT_SPACING_TOL`). A link with no
/// such return, a unit of another length, or a unit that is the whole circle is left as it is: a
/// read longer than a circular genome crosses its links twice, a turn apart. None when nothing is
/// added; else the new walk, per unit the copies added, and whether a linear walk was closed.
///
/// A linear walk whose missing unit, closed by its link, passes every segment of the walk is that
/// circle: the reads crossed the link a whole genome apart. Nipponbare ONT plastid (981c9d0): the
/// inverted repeat came out at 1.18x the LSC depth, the bridges laid it once in the linear
/// u0+ u1- u4- u2- u3- (LSC, IR, SSC), and reads crossing u0+ -> u1- twice, 134 kb apart, laid
/// the whole circle twice more: a 382.6 kb linear molecule, 2.85 genomes. Closed, it is the
/// 134.4 kb circle u0+ u1- u4- u2- u3- u2+ u4+ u1+.
///
/// The bridges lay a tandem unit as often as its depth says (`copy_deviation`), and an anchor once.
/// Where the reads are unevenly deep along the genome that undercounts: Haliclona indistincta
/// (ERR12666651, PMAT2 panel, HiFi) has read depth 150-175x over most of its 24.9 kb
/// mitochondrion and 13x at its control region, a ~5 copy array of a 986 bp unit (u1 793 bp,
/// u2 652 bp, 459 bp overlap at k=501). The unit's segments came out at 1.36x and 1.05x of the
/// median depth, u2 an anchor, and the 20,232 bp circle laid the unit twice (reference coverage
/// 0.81) although single reads crossed its junction up to 4 times (9 reads; 10 reads 3 times).
/// From another draw of its reads the same graph gave the linear u0+ u1- u2+, which passes the
/// unit junction once and holds no cycle of its own: the return u2+ -> u1- comes from the graph.
///
/// A unit of mostly single-copy sequence (`repeat` false for more than half its bases) is not
/// laid: those segments stand at the depth of one copy, so a read that crosses the link again a
/// unit later went round a subcircle, it did not read a second copy. Nipponbare ONT mitochondrion
/// (ont-sr-fix read pool, 11-15): seven units of 59-131 kb, 85-100% of them single-copy segments
/// at 0.67-1.30x the median depth, each from 3-5 reads with two crossings, were laid once each:
/// 997 kb of molecules for a 490 kb genome. Haliclona's 986 bp unit is 45% single copy (u2 652 bp;
/// u1 at 1.36x is a repeat) and is still laid. A linear walk that a unit closes is closed first.
pub(crate) fn lay_read_copies(
    g: &Graph,
    walk: &[Oriented],
    circular: bool,
    need: &FxHashMap<(Oriented, Oriented), (usize, usize)>,
    next: &FxHashMap<Oriented, Vec<Oriented>>,
    overlaps: &FxHashMap<(Oriented, Oriented), usize>,
    repeat: &[bool],
) -> Option<(Vec<Oriented>, Vec<(Vec<Oriented>, usize)>, bool)> {
    if need.is_empty() || walk.is_empty() {
        return None;
    }
    let mut w = walk.to_vec();
    let mut circular = circular;
    let mut closed = false;
    let mut added: Vec<(Vec<Oriented>, usize)> = Vec::new();
    let mut total = 0;
    'grow: while total < MAX_READ_COPIES_ADDED {
        let n = w.len();
        let steps = if circular { n } else { n - 1 };
        let mut passes: FxHashMap<(Oriented, Oriented), usize> = FxHashMap::default();
        for t in 0..steps {
            *passes.entry(link_key(w[t], w[(t + 1) % n])).or_default() += 1;
        }
        for t in 0..steps {
            let (x, y) = (w[t], w[(t + 1) % n]);
            let key = link_key(x, y);
            let Some(&(want, spacing)) = need.get(&key) else {
                continue;
            };
            if passes[&key] >= want {
                continue;
            }
            let Some(unit) = shortest_return(next, y, x) else {
                continue;
            };
            // bases one more unit adds: each of its segments less its overlap with the one before
            let mut prev = x;
            let mut adds = 0usize;
            for &o in &unit {
                let ov = overlaps.get(&(prev, o)).copied().unwrap_or(0);
                adds += g.seqs[o.0].len().saturating_sub(ov);
                prev = o;
            }
            adds = adds.saturating_sub(overlaps.get(&(prev, y)).copied().unwrap_or(0))
                + overlaps.get(&(x, y)).copied().unwrap_or(0);
            if (adds as f64 - spacing as f64).abs() > (UNIT_SPACING_TOL * spacing as f64).max(50.0) {
                continue;
            }
            // a circle that is its own unit is the molecule (a 60 bp circle of a short tandem
            // repeat stays 60 bp)
            if circular && unit.len() == n && canonical_key(&unit, true) == canonical_key(&w, true) {
                continue;
            }
            if !circular {
                let ring: FxHashSet<usize> = unit.iter().map(|o| o.0).collect();
                if w.iter().all(|o| ring.contains(&o.0)) {
                    // the unit runs y .. x and the link x -> y closes it; start where the walk did
                    let mut c = unit;
                    if let Some(i) = c.iter().position(|&o| o == w[0]) {
                        c.rotate_left(i);
                    }
                    w = c;
                    circular = true;
                    closed = true;
                    continue 'grow;
                }
            }
            let single: usize = unit
                .iter()
                .filter(|o| !repeat.get(o.0).copied().unwrap_or(true))
                .map(|o| g.seqs[o.0].len())
                .sum();
            let all: usize = unit.iter().map(|o| g.seqs[o.0].len()).sum();
            if 2 * single > all {
                continue;
            }
            w.splice(t + 1..t + 1, unit.iter().copied());
            match added.iter_mut().find(|(u, _)| *u == unit) {
                Some((_, c)) => *c += 1,
                None => added.push((unit, 1)),
            }
            total += 1;
            continue 'grow;
        }
        break;
    }
    (!added.is_empty() || closed).then_some((w, added, closed))
}

/// Order of two covers of one component, across cover sizes; `Less` ranks first.
/// Fewer linear molecules first: a linear molecule lacks its closing bridge, whose ln share is
/// never positive, so scores alone would favour opening circles. Prefer covers whose every
/// junction has sufficient read support. Where all covers need unphased or weak junctions,
/// their synthetic weights cannot resolve molecule count: choose the smallest feasible
/// representative, while retaining the larger covers as ambiguous alternatives. Fully
/// read-supported covers compete by score across molecule counts.
fn cover_order(a: &CoverReport, b: &CoverReport) -> Ordering {
    if let (Some(x),Some(y))=(&a.read_support,&b.read_support) {
        return x.order(y).then(a.copy_deviation.cmp(&b.copy_deviation))
            // Representation-only tiebreak, never a biological evidence advantage.
            .then(a.molecules.len().cmp(&b.molecules.len()));
    }
    a.linear()
        .cmp(&b.linear())
        .then(
            a.incomplete_junction_evidence()
                .cmp(&b.incomplete_junction_evidence()),
        )
        .then_with(|| {
            if a.incomplete_junction_evidence() && b.incomplete_junction_evidence() {
                a.molecules.len().cmp(&b.molecules.len())
            } else {
                Ordering::Equal
            }
        })
        .then(b.ranked_score().total_cmp(&a.ranked_score()))
        // A deterministic representative when the evidence ties, not extra evidence.
        .then(a.molecules.len().cmp(&b.molecules.len()))
}

/// Anchors joined by bridges, and the molecules unringing finds for them.
#[derive(Serialize)]
pub struct ComponentReport {
    pub ranking_model: &'static str,
    /// Ties on measured record/link support, ignoring topology and copy-number preferences.
    pub evidence_tied_candidates: usize,
    pub read_model_budget_exhausted: bool,
    pub read_path_cover: Option<crate::read_model::ReadPathCover>,
    /// An explanatory family of candidate configurations, not certified physical isoforms.
    pub witness_configurations: Vec<CoverReport>,
    pub ambiguity_configurations: Vec<CoverReport>,
    pub anchors: Vec<String>,
    /// Length-weighted median depth of the anchors, when the evidence records depth.
    pub depth: Option<f64>,
    /// Molecules in the selected cover, plus the tails cut off (0 when unsolved or
    /// skipped).
    pub molecules: usize,
    /// Required anchors bridged at one end only: each must end a linear molecule, so the
    /// component needs at least half as many molecules (`min_molecules`). Where the graph is
    /// broken (a gap, an unresolved error structure) these are the places to look.
    pub path_ends: Vec<String>,
    pub min_molecules: usize,
    /// Every anchor is optional (shallow): a low-frequency configuration left out.
    pub skipped: bool,
    pub candidates: usize,
    pub search_truncated: bool,
    /// Whether a bounded retry was performed in automatic budget mode.
    pub adaptive_retry: bool,
    /// Largest cover size enumerated for this component (zero when skipped).
    pub searched_max_molecules: usize,
    /// Larger cover sizes, up to one molecule per anchor, were excluded by the configured
    /// molecule-count limit. Even a complete expansion search cannot rule those covers out.
    pub cover_size_limited: bool,
    /// The selected cover is the read-joined paths (`partial_paths`) of at least
    /// `MIN_UNBRIDGED_BP`, taken because the search ran out of budget before finding any cover of
    /// the sizes it tried. Those paths visit every
    /// required anchor once through read-backed bridges, so they are a cover of their own size;
    /// a smaller one may exist unseen, and the component is not decisive.
    pub read_join_cover: bool,
    pub best: Option<CoverReport>,
    /// For a component without a cover: the paths its read bridges make (see `partial_paths`),
    /// so that its sequence is reported. Not a solution: `best` stays empty, `molecules` 0, and
    /// the component counts among `unsolved_components`. With `read_join_cover`, the read-joined
    /// paths under `MIN_UNBRIDGED_BP` (the longer ones are `best`).
    pub partial: Option<CoverReport>,
    /// Bridge/copy score lead of the selected representative over its strongest competitor
    /// with as many linear molecules, across all searched cover sizes. May be negative when
    /// incomplete junction evidence requires a parsimonious representative.
    pub margin_ln: Option<f64>,
    pub decisive: bool,
    /// Bridges left out where the best cover's linear molecule led back into itself and is
    /// reported as a circle and its tail (see `cut_tail`); reads support these bridges too, so
    /// a component with one is not decisive. `alternatives` and `margin_ln` refer to the
    /// covers before the cut.
    pub tails_cut: Vec<String>,
    /// Tandem units the best cover lays down more often than its bridges do, because single
    /// reads cross the unit's junction more often (see `lay_read_copies`): "unit walk xN" per
    /// extra copy count. `alternatives` keep the bridges' copies.
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub read_copies_added: Vec<String>,
    pub alternatives: Vec<CoverReport>,
    /// With a `BridgePrior`: the closed cover it takes among those the reads rank alike. Not
    /// among `candidates`: it is built, not searched, so a truncated search does not miss it.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub prior_cover: Option<PriorCover>,
}

/// Summarize only ranked components; unplaced/skipped fragments have no ranking.
fn effective_ranking_model<'a>(models: impl Iterator<Item = &'a str>) -> &'static str {
    let mut native = false;
    let mut legacy = false;
    let mut fallback = false;
    for model in models {
        match model {
            "legacy_after_read_budget_exhaustion" => fallback = true,
            "legacy_bridge_copy" => legacy = true,
            _ => native = true,
        }
    }
    match (native, legacy, fallback) {
        (false, false, false) => "no_ranked_components",
        (true, false, false) => crate::read_model::MODEL,
        (false, true, false) => "legacy_bridge_copy",
        (false, false, true) => "legacy_after_read_budget_exhaustion",
        _ => "mixed_component_models",
    }
}

#[derive(Serialize)]
pub struct LinearizeReport {
    /// Requested/available model; see effective_ranking_model for actual component use.
    pub ranking_model: &'static str,
    /// Actual ranking used across non-skipped components with ranked candidates.
    pub effective_ranking_model: &'static str,
    /// True if any component abandoned read-based scoring because its budget ran out.
    pub read_model_budget_exhausted: bool,
    pub physical_topology_identified: bool,
    pub configurations_fasta: Option<String>,
    /// Requested FASTA of ranked component covers.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub candidate_fasta: Option<String>,
    pub graph: String,
    pub evidence: String,
    pub anchors: Vec<String>,
    pub repeats: Vec<String>,
    pub bridges: usize,
    /// Of those, pairings through repeats no read crosses, taken from the graph alone.
    pub topological_bridges: usize,
    /// Of those, pairings taken at anchor ends no bridge of `min_reads` leaves.
    pub weak_bridges: usize,
    /// Self-loops and paths returning to their own anchor (tandem or sub-circle evidence).
    pub self_loops: Vec<String>,
    /// Local repeat self-loop hypotheses do not participate in complete covers.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub local_repeat_closures: Option<LocalRepeatClosures>,
    /// Anchors without any bridge: no path joins them to the rest, so each is reported as its own
    /// linear molecule (unless shallow); the graph is broken there or the reads do not cover it.
    pub unbridged_anchors: Vec<String>,
    /// Shallow or weakly bridged anchors (see `optional_below`): visited only when needed.
    pub optional_anchors: Vec<String>,
    /// Optional anchors neither a molecule nor a partial path holds (low-frequency
    /// configurations).
    pub skipped_anchors: Vec<String>,
    /// One entry per connected component of the anchor graph, largest first.
    pub components: Vec<ComponentReport>,
    /// Components no cover of up to `max_molecules` molecules solves (their sequence is in
    /// `partial_paths`).
    pub unsolved_components: usize,
    /// Molecules over all components (the FASTA records).
    pub molecules: usize,
    /// What no molecule holds but reads still join or the graph requires: the partial paths of
    /// unsolved components and, each alone, the required anchors without any bridge; longest
    /// first, the records of `partial_fasta`. Incomplete by definition, and never counted in
    /// `molecules`.
    pub partial_paths: Vec<Linearization>,
    /// FASTA of `partial_paths` (records `partial.N`, tagged `complete=false`); not written
    /// when there are none.
    pub partial_fasta: Option<String>,
    /// Segments in neither a molecule nor a partial path: optional anchors left out, shallow
    /// components, repeats no chosen bridge crosses.
    pub unplaced_segments: Vec<String>,
    pub candidates: usize,
    pub search_truncated: bool,
    /// At least one searched component omitted cover sizes because of `max_molecules`.
    pub cover_size_limited: bool,
    /// The longest molecule; `alternatives` and `margin_ln` refer to its component.
    pub best: Option<Linearization>,
    /// Bridge/copy score lead of the main component's representative over its strongest
    /// competitor of the same circularity (None with no such competitor).
    pub margin_ln: Option<f64>,
    /// Every searched component has no omitted cover sizes, a complete expansion search,
    /// no weak or topological bridges, and a score lead of
    /// at least ln(2) over a runner-up with the same number of linear molecules (when present).
    /// This is confidence within the bridge/copy score and configured cover-size limit,
    /// not certification of physical chromosome count.
    pub decisive: bool,
    pub alternatives: Vec<Linearization>,
    /// Repeats no molecule lays down (no read-supported path crosses them).
    pub repeats_not_placed: Vec<String>,
}

/// FASTA records `{prefix}` (one record) or `{prefix}.N`, each with its topology, length and
/// path, and `tag` appended to the header.
fn fasta_records(
    g: &Graph,
    records: &[(Vec<Oriented>, Linearization)],
    prefix: &str,
    number_single: bool,
    tag: &str,
) -> Result<String> {
    let mut fa = String::new();
    for (i, (walk, lin)) in records.iter().enumerate() {
        let seq = spell(g, walk, lin.circular)?;
        let name = if records.len() == 1 && !number_single {
            prefix.to_string()
        } else {
            format!("{prefix}.{}", i + 1)
        };
        fa.push_str(&format!(
            ">{name} circular={} length={} path={}{tag}\n",
            lin.circular,
            seq.len(),
            lin.path.replace(' ', ",")
        ));
        for chunk in seq.chunks(80) {
            fa.push_str(std::str::from_utf8(chunk)?);
            fa.push('\n');
        }
    }
    Ok(fa)
}

/// Where the partial paths go unless told otherwise: beside `out_fasta`, `.partial` before
/// its extension (`molecules.fasta` -> `molecules.partial.fasta`).
fn partial_fasta_path(out_fasta: &Path) -> std::path::PathBuf {
    let stem = out_fasta.file_stem().unwrap_or_default().to_string_lossy();
    let name = match out_fasta.extension() {
        Some(ext) => format!("{stem}.partial.{}", ext.to_string_lossy()),
        None => format!("{stem}.partial"),
    };
    out_fasta.with_file_name(name)
}

fn configuration_fasta_path(out_fasta: &Path) -> PathBuf {
    let stem = out_fasta.file_stem().unwrap_or_default().to_string_lossy();
    out_fasta.with_file_name(format!("{stem}.configurations.fasta"))
}

/// Resolve existing symlinks and compare paths even before the output directories exist.
fn candidate_path_identity(path: &Path) -> Result<PathBuf> {
    let mut ancestor = std::path::absolute(path)?;
    let mut suffix = Vec::new();
    loop {
        match fs::canonicalize(&ancestor) {
            Ok(mut resolved) => {
                resolved.extend(suffix.into_iter().rev());
                return Ok(resolved);
            }
            Err(error)
                if error.kind() == std::io::ErrorKind::NotFound
                    && fs::symlink_metadata(&ancestor).is_err() =>
            {
                let name = ancestor.file_name().with_context(|| {
                    format!("cannot resolve output location {}", path.display())
                })?;
                suffix.push(name.to_os_string());
                ancestor.pop();
            }
            Err(error) => return Err(error).with_context(|| format!("resolving {}", path.display())),
        }
    }
}

fn cover_fasta_records(
    g: &Graph,
    cover: &CoverReport,
    label_index: &FxHashMap<&str, usize>,
    prefix: &str,
    tag: &str,
) -> Result<String> {
    let records = cover.molecules.iter().map(|molecule| {
        let walk = molecule.path.split_whitespace()
            .map(|label| parse_label(label_index, label))
            .collect::<Result<Vec<_>>>()?;
        Ok((walk, molecule.clone()))
    }).collect::<Result<Vec<_>>>()?;
    fasta_records(g, &records, prefix, true, tag)
}

fn write_candidate_fasta(
    g: &Graph,
    components: &[ComponentReport],
    label_index: &FxHashMap<&str, usize>,
    path: &Path,
) -> Result<()> {
    static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);
    let parent = path.parent().unwrap_or(Path::new("."));
    fs::create_dir_all(parent)?;
    let pending = parent.join(format!(
        ".ovasm-candidates-{}-{}.tmp",
        std::process::id(),
        NEXT_TEMP.fetch_add(1, AtomicOrdering::Relaxed)
    ));
    let file = fs::OpenOptions::new().write(true).create_new(true).open(&pending)?;
    let result = (|| -> Result<()> {
        let mut writer = BufWriter::new(file);
        for (ci, component) in components.iter().enumerate() {
            let ranked_count = component.alternatives.len()
                + usize::from(component.best.is_some() && !component.read_join_cover);
            let tag = format!(
                " representation=ranked_component_cover physical_topology_identified=false search_truncated={} cover_size_limited={} read_model_budget_exhausted={} alternatives_truncated={}",
                component.search_truncated,
                component.cover_size_limited,
                component.read_model_budget_exhausted,
                ranked_count < component.candidates,
            );
            for (rank, cover) in component.best.iter().chain(&component.alternatives).enumerate() {
                if rank == 0 && component.read_join_cover {
                    continue;
                }
                let prefix = format!("component_candidate.{}.{}", ci + 1, rank + 1);
                let fasta = cover_fasta_records(g, cover, label_index, &prefix, &tag)?;
                writer.write_all(fasta.as_bytes())?;
            }
        }
        writer.flush()?;
        drop(writer);
        fs::rename(&pending, path)?;
        Ok(())
    })();
    if result.is_err() {
        let _ = fs::remove_file(&pending);
    }
    result.with_context(|| format!("writing ranked candidates to {}", path.display()))
}

/// `out_partial`: FASTA for the partial paths (default: see `partial_fasta_path`). A file left
/// there by an earlier run is removed when this run has none, so that it cannot be taken for
/// this run's.
pub fn run(
    gfa: &Path,
    evidence: &Path,
    p: LinearizeParams,
    out_fasta: &Path,
    out_partial: Option<&Path>,
    out_json: &Path,
) -> Result<LinearizeReport> {
    run_inner(gfa, evidence, p, out_fasta, out_partial, out_json, None, None, None)
}

pub fn run_phased(
    gfa: &Path,
    evidence: &Path,
    p: LinearizeParams,
    out_fasta: &Path,
    out_partial: Option<&Path>,
    out_json: &Path,
    phase_report: Option<&Path>,
    out_candidate: Option<&Path>,
) -> Result<LinearizeReport> {
    run_inner(gfa, evidence, p, out_fasta, out_partial, out_json, phase_report, out_candidate, None)
}

/// `run`, and per component the closed cover `prior` takes among those the reads rank alike
/// (`ComponentReport::prior_cover`). The search, its ranking and the files are as `run` gives
/// them.
pub(crate) fn run_with_prior(
    gfa: &Path,
    evidence: &Path,
    p: LinearizeParams,
    out_fasta: &Path,
    out_partial: Option<&Path>,
    out_json: &Path,
    prior: &BridgePrior,
) -> Result<LinearizeReport> {
    run_inner(gfa, evidence, p, out_fasta, out_partial, out_json, None, None, Some(prior))
}

#[allow(clippy::too_many_arguments)]
fn run_inner(
    gfa: &Path,
    evidence: &Path,
    p: LinearizeParams,
    out_fasta: &Path,
    out_partial: Option<&Path>,
    out_json: &Path,
    phase_report: Option<&Path>,
    out_candidate: Option<&Path>,
    prior: Option<&BridgePrior>,
) -> Result<LinearizeReport> {
    let partial_fasta = out_partial.map_or_else(|| partial_fasta_path(out_fasta), Path::to_path_buf);
    let configuration_path = configuration_fasta_path(out_fasta);
    if let Some(candidate) = out_candidate {
        let identity = candidate_path_identity(candidate)?;
        for occupied in [gfa, evidence, out_fasta, out_json, &partial_fasta, &configuration_path]
            .into_iter().chain(phase_report)
        {
            if candidate_path_identity(occupied)? == identity {
                bail!("candidate FASTA must be separate from {}", occupied.display());
            }
        }
    }
    let g = parse_gfa(gfa)?;
    let phases = phase_report
        .map(|p| crate::phase::Constraints::load(&g, p))
        .transpose()?
        .unwrap_or_default();
    let report: Value = serde_json::from_str(
        &fs::read_to_string(evidence)
            .with_context(|| format!("cannot read {}", evidence.display()))?,
    )?;
    crate::repeat_depth::validate_report_graph(gfa, &report)?;
    let repeat_names: BTreeSet<&str> = report["repeat_segments"]
        .as_array()
        .context("evidence has no repeat_segments")?
        .iter()
        .filter_map(Value::as_str)
        .collect();
    let repeat: Vec<bool> = g
        .names
        .iter()
        .map(|n| repeat_names.contains(n.as_str()))
        .collect();
    let anchors: Vec<usize> = (0..g.names.len()).filter(|&s| !repeat[s]).collect();
    let read_model=crate::read_model::ReadModel::load(&g,&report,&repeat,p.min_reads)?;
    if g.names.is_empty() {
        bail!("the assembly graph has no segments; the reads gave nothing to assemble");
    }
    if anchors.is_empty() {
        bail!("every segment is a repeat; nothing anchors a linear path");
    }
    let (mut br, self_loops) = bridges(&g, &report, &repeat, p.min_reads)?;
    br.retain(|b| phases.allows_walk(&b.walk, false));
    let spanned: Vec<(Vec<Oriented>, u64)> = spanned_walks(&g, &report, &repeat, p.min_reads)?
        .into_iter()
        .filter(|(w, _)| phases.allows_walk(w, false))
        .collect();
    let spanned = merge_spanned(&mut br, &spanned);
    // Count each bridge once per distinct endpoint, including self-returning bridges.
    // Scanning all bridges for each anchor is quadratic on fragmented short-read graphs.
    let mut endpoint_support = vec![0u64; g.names.len()];
    for b in &br {
        let (a, z) = (b.walk[0].0, b.walk[b.walk.len() - 1].0);
        endpoint_support[a] += b.reads;
        if z != a {
            endpoint_support[z] += b.reads;
        }
    }
    let support = |a: usize| endpoint_support[a];
    let mut totals: Vec<u64> = anchors
        .iter()
        .map(|&a| support(a))
        .filter(|&t| t > 0)
        .collect();
    totals.sort_unstable();
    let median = totals.get(totals.len() / 2).copied().unwrap_or(0) as f64;
    let mut optional = vec![false; g.names.len()];
    let anchor_components = components(g.names.len(), &br);
    let mut component_depth: Vec<Option<f64>> = vec![None; anchor_components.len()];
    match (
        report["segment_depths"].as_object(),
        report["median_depth"].as_f64(),
    ) {
        // Low-frequency configurations are shallow: judge by depth when the graph records it.
        (Some(depths), Some(md)) => {
            let depth = |a: usize| {
                depths
                    .get(&g.names[a])
                    .and_then(Value::as_f64)
                    .unwrap_or(md)
            };
            for &a in &anchors {
                optional[a] = depth(a) < p.optional_below * md;
            }
            // Shallow against what? Chromosomes differ in copy number, so the graph's median
            // is the wrong yardstick for an anchor of a rarer chromosome: Trifolium dubium's
            // two smaller chromosomes sit at 0.50-0.74 of the largest (182 kb, over half the
            // graph), and a 23 kb anchor at 0.499 of the median was left out as a
            // low-frequency configuration although it is as deep as its neighbours. A
            // configuration is rare relative to the anchors it is bridged with: where a
            // component is shallower than the anchors of the graph as a whole, the threshold
            // comes down by that ratio. (With one component the ratio is 1.) A component
            // shallow as a whole, its ratio below the fraction, keeps the graph's threshold:
            // such parts are left out as before (Bromus sterilis: components of 0.8-3.6 kb
            // at 0.2-0.4 of the median, none of them in the published mitochondrion).
            let bridged: Vec<usize> = anchor_components.iter().flatten().copied().collect();
            let all = weighted_median_depth(&g, depth, &bridged);
            for (comp, own) in anchor_components.iter().zip(component_depth.iter_mut()) {
                let d = weighted_median_depth(&g, depth, comp);
                *own = Some(d);
                if d < all && d >= p.optional_below * all {
                    for &a in comp {
                        optional[a] = depth(a) < p.optional_below * md * (d / all);
                    }
                }
            }
        }
        // Otherwise by bridge support, which also drops for anchors behind long repeats, so
        // only much weaker ones qualify.
        _ => {
            for &a in &anchors {
                optional[a] = (support(a) as f64) < 0.5 * p.optional_below * median;
            }
        }
    }
    if anchors.iter().all(|&a| optional[a]) {
        optional.iter_mut().for_each(|o| *o = false);
    }
    // copies each repeat should be laid down, from its depth, or from the reads that span it
    let mut expected: Vec<Option<usize>> = match (
        report["segment_depths"].as_object(),
        report["median_depth"].as_f64(),
    ) {
        (Some(depths), Some(md)) if md > 0.0 => (0..g.names.len())
            .map(|s| {
                let d = depths.get(&g.names[s]).and_then(Value::as_f64)?;
                repeat[s].then(|| ((d / md).round() as usize).max(1))
            })
            .collect(),
        _ => vec![None; g.names.len()],
    };
    // passes per spanned anchor pair (its walks are alternatives: the most of any), summed over
    // the pairs (each a place the repeat is laid down)
    let mut per_pair: BTreeMap<(Oriented, Oriented), FxHashMap<usize, usize>> = BTreeMap::new();
    for w in &spanned {
        let mut passes: FxHashMap<usize, usize> = FxHashMap::default();
        for o in &w[1..w.len() - 1] {
            *passes.entry(o.0).or_default() += 1;
        }
        let most = per_pair.entry((w[0], w[w.len() - 1])).or_default();
        for (s, n) in passes {
            let e = most.entry(s).or_default();
            *e = (*e).max(n);
        }
    }
    let mut by_span: FxHashMap<usize, usize> = FxHashMap::default();
    for most in per_pair.values() {
        for (&s, &n) in most {
            *by_span.entry(s).or_default() += n;
        }
    }
    for (s, n) in by_span {
        expected[s] = Some(n);
    }
    let read_copies = read_copies(&g, &report, p.min_reads)?;
    let mut graph_next: FxHashMap<Oriented, Vec<Oriented>> = FxHashMap::default();
    for &(x, y) in &g.links {
        graph_next.entry(x).or_default().push(y);
        graph_next.entry(flip(y)).or_default().push(flip(x));
    }
    let frac = exit_fractions(&br);
    let label = |walk: &[Oriented]| {
        walk.iter()
            .map(|&o| g.label(o))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let indexed_overlaps = overlap_index(&g);
    let describe = |f: &Found| -> Result<(Vec<Oriented>, Linearization)> {
        let walk = candidate_walk(&br, f);
        let fractions: Vec<f64> = f.bridges.iter().map(|&i| frac[i]).collect();
        let lin = Linearization {
            path: label(&walk),
            circular: f.circular,
            length: spelled_length(&g, &walk, f.circular, &indexed_overlaps)?,
            score: path_score(&fractions),
            bridge_reads: f.bridges.iter().map(|&i| br[i].reads).collect(),
            exit_fractions: fractions,
        };
        Ok((walk, lin))
    };
    // a set of molecules as reported: longest first, each with its walk
    let cover_report = |cover: &[Found]| -> Result<(Vec<Vec<Oriented>>, CoverReport)> {
        let mut mols = cover.iter().map(&describe).collect::<Result<Vec<_>>>()?;
        mols.sort_by(|(wa, a), (wb, b)| b.length.cmp(&a.length).then_with(|| wa.cmp(wb)));
        let walks = mols.iter().map(|(w, _)| w.clone()).collect();
        let molecules: Vec<Linearization> = mols.into_iter().map(|(_, l)| l).collect();
        let report = CoverReport {
            score: molecules.iter().map(|l| l.score).sum(),
            molecules,
            weak_bridges: cover
                .iter()
                .flat_map(|f| &f.bridges)
                .filter(|&&i| br[i].weak)
                .count(),
            topological_bridges: cover
                .iter()
                .flat_map(|f| &f.bridges)
                .filter(|&&i| br[i].topological)
                .count(),
            copy_deviation: copy_deviation(&expected, &br, cover),
            read_support: None,
        };
        Ok((walks, report))
    };

    let mut comps: Vec<ComponentReport> = Vec::new();
    // every component's best molecules, for spelling and placement checks
    let mut chosen: Vec<(Vec<Oriented>, Linearization)> = Vec::new();
    // what no cover explains, as partial paths
    let mut partial: Vec<(Vec<Oriented>, Linearization)> = Vec::new();
    let bridge_starts: FxHashSet<Oriented> = br.iter().map(|b| b.walk[0]).collect();
    let mut search_out_of: FxHashMap<Oriented, Vec<usize>> = FxHashMap::default();
    for (i, b) in br.iter().enumerate() {
        search_out_of.entry(b.walk[0]).or_default().push(i);
    }
    for (comp, depth) in anchor_components.into_iter().zip(component_depth) {
        let component_model=read_model.as_ref().map(|m|m.for_anchors(&comp,&repeat));
        let thread_ends = component_model.as_ref().map(|model| model.thread_ends(&comp, &optional));
        let mut rep = ComponentReport {
            ranking_model: if read_model.is_some() {crate::read_model::MODEL} else {"legacy_bridge_copy"},
            evidence_tied_candidates: 0,
            read_model_budget_exhausted: false,
            read_path_cover: None,
            witness_configurations: Vec::new(),
            ambiguity_configurations: Vec::new(),
            anchors: comp.iter().map(|&a| g.names[a].clone()).collect(),
            depth,
            molecules: 0,
            path_ends: Vec::new(),
            min_molecules: 0,
            skipped: comp.iter().all(|&a| optional[a]),
            candidates: 0,
            search_truncated: false,
            adaptive_retry: false,
            searched_max_molecules: 0,
            cover_size_limited: false,
            read_join_cover: false,
            best: None,
            partial: None,
            margin_ln: None,
            decisive: false,
            tails_cut: Vec::new(),
            read_copies_added: Vec::new(),
            prior_cover: None,
            alternatives: Vec::new(),
        };
        if rep.skipped {
            comps.push(rep);
            continue;
        }
        let bridged_at = |o: Oriented| bridge_starts.contains(&o);
        let ends: Vec<usize> = comp
            .iter()
            .copied()
            .filter(|&a| !optional[a] && bridged_at((a, true)) != bridged_at((a, false)))
            .collect();
        rep.path_ends = ends.iter().map(|&a| g.names[a].clone()).collect();
        rep.min_molecules = ends.len().div_ceil(2).max(1);
        // Compare every feasible cover size under the same bridge/copy score. A component
        // sharing repeats can hold several chromosomes; finding a single-molecule cover
        // does not establish that it is the representative best supported by the reads.
        let mut ranked: Vec<(Vec<Vec<Oriented>>, CoverReport, Vec<Found>)> = Vec::new();
        let mut seen = BTreeSet::new();
        // An automatic envelope must include the topological lower bound:
        // otherwise the feasible-size loop is empty before any search.
        // Explicit caller limits remain authoritative, and omitted larger
        // sizes still prohibit a decisive claim. No support threshold changes.
        let envelope = if p.max_molecules == 0 {
            if read_model.is_some() && comp.len() <= 12 {comp.len()} else {4.max(rep.min_molecules)}
        } else {
            p.max_molecules
        };
        // An automatic envelope that holds no cover widens one size at a time until one does:
        // a branched component of an incomplete large mitochondrion (Lophophora williamsii,
        // 996 kb in 15 segments, four molecules at least by its path ends) needed six, and all
        // of it went to the partial paths. Components with a cover inside the envelope are
        // searched as before; explicit limits stay authoritative.
        let mut last = envelope.min(comp.len());
        let mut m = rep.min_molecules;
        // A requested larger envelope can exhaust work at several cover sizes.
        // Allow at most one retry per searched size, capped at six in total;
        // default and compact-component searches retain one retry.
        let retry_limit = if comp.len() > 32 && p.max_molecules > 4 {
            p.max_molecules.min(6)
        } else { 1 };
        let mut retries = 0;
        while m <= last {
            let (mut found, mut truncated, mut terminal_truncated) = search_phased_indexed(
                g.names.len(), &comp, &optional, &br, m, p.max_steps, &phases, &search_out_of, thread_ends.as_ref(),
            );
            // Retain the existing compact-component retry. For moderately sized
            // components, retry once with a bounded candidate count so that a larger
            // step budget does not turn into unbounded candidate memory.
            // Explicit --max-steps remains authoritative.
            let compact_retry = comp.len() <= 32 && found.len() <= 1_000_000;
            let bounded_retry = comp.len() > 32 && comp.len() <= 64 && found.len() <= 1_000_000;
            if (truncated || terminal_truncated) && p.adaptive_steps && retries < retry_limit
                && ((compact_retry && br.len() <= 128) || (bounded_retry && br.len() <= 256))
            {
                rep.adaptive_retry = true;
                retries += 1;
                drop(found);
                let factor = if bounded_retry && p.max_molecules > 5 && m >= 6 { 20 }
                    else if bounded_retry && p.max_molecules > 4 && m >= 5 { 10 }
                    else if bounded_retry { 4 } else { 2 };
                let cap = bounded_retry.then_some(if p.max_molecules > 5 && m >= 6 {
                    2_000_000
                } else { 1_000_000 });
                (found, truncated, terminal_truncated) = search_phased_indexed_bounded(
                    g.names.len(), &comp, &optional, &br, m,
                    p.max_steps.saturating_mul(factor), &phases, &search_out_of, thread_ends.as_ref(),
                    true, cap,
                );
            }
            rep.search_truncated |= truncated || terminal_truncated;
            for cover in found {
                let mut key: Vec<(bool, Vec<Oriented>)> = cover.iter()
                    .map(|f| (f.circular, canonical_key(&candidate_walk(&br, f), f.circular)))
                    .collect();
                key.sort();
                if seen.insert(key) {
                    let (walks, report) = cover_report(&cover)?;
                    ranked.push((walks, report, cover));
                }
            }
            // A truncated search cannot establish that this cover size has no solution.
            // Do not multiply an exhausted budget by trying every larger cover size.
            if truncated {
                last = m;
                break;
            }
            if m == last && ranked.is_empty() && p.max_molecules == 0 && last < comp.len() {
                last += 1;
            }
            m += 1;
        }
        rep.searched_max_molecules = last;
        rep.cover_size_limited = rep.searched_max_molecules < comp.len();
        if ranked.is_empty() {
            let paths = partial_paths_phased(g.names.len(), &comp, &optional, &br, &frac, &phases);
            let (walks, report) = cover_report(&paths)?;
            // A search that ran out of budget has not shown that no cover exists, and the
            // read-joined paths are one: every required anchor once, joined only by bridges
            // reads took. With no molecule limit set, they are the cover found. Otherwise a
            // better-linked graph would come out with less than a more broken one, whose pieces
            // are covers of single anchors (rice SAMEA115396709's mitochondrion, 34 anchors
            // holding 252 of its 316 kb: truncated at four molecules, all of it partial, 54 kb
            // of molecules against 287 kb from a graph in more pieces).
            // Only paths of at least `MIN_UNBRIDGED_BP` are molecules, as for an anchor no bridge
            // reaches: shorter ones are the stubs short-read graphs leave, and stay partial (rice
            // SAMEA115396710 from 150 bp reads: a nuclear component of 1,686 anchors at 0.58 times
            // the mitochondrion's depth came out as 853 molecules, 122.5 kb, all under 1 kb).
            let (long, short): (Vec<Found>, Vec<Found>) = paths.iter().cloned().partition(|f| {
                describe(f).is_ok_and(|(_, l)| l.length >= MIN_UNBRIDGED_BP)
            });
            if rep.search_truncated && p.max_molecules == 0 && !long.is_empty() {
                rep.read_join_cover = true;
                if !short.is_empty() {
                    let (walks, report) = cover_report(&short)?;
                    partial.extend(walks.into_iter().zip(report.molecules.iter().cloned()));
                    rep.partial = Some(report);
                }
                let (walks, report) = cover_report(&long)?;
                ranked.push((walks, report, long));
            } else {
                partial.extend(walks.into_iter().zip(report.molecules.iter().cloned()));
                rep.partial = Some(report);
            }
        }
        if let Some(model)=&component_model {
            let mut budget=scoring_work_budget(comp.len(), ranked.len());
            for (walks,cover,_) in &mut ranked {
                let circular:Vec<_>=cover.molecules.iter().map(|m|m.circular).collect();
                cover.read_support=model.score(walks,&circular,&mut budget);
                if cover.read_support.is_none() {rep.read_model_budget_exhausted=true;break;}
            }
            if rep.read_model_budget_exhausted {
                for (_,cover,_) in &mut ranked {cover.read_support=None;}
                rep.ranking_model="legacy_after_read_budget_exhaustion";
            }
        }
        ranked.sort_by(|(wa, a, _), (wb, b, _)| cover_order(a, b).then_with(|| wa.cmp(wb)));
        rep.candidates = ranked.len();
        rep.margin_ln = ranked.first().and_then(|(_, best, _)| {
            ranked
                .iter()
                .skip(1)
                .filter(|(_, other, _)| other.linear() == best.linear())
                .map(|(_, other, _)| other.ranked_score())
                .max_by(f64::total_cmp)
                .map(|other| best.ranked_score() - other)
        });
        // A partial search cannot establish the score margin over unseen covers. Omitted
        // cover sizes, weak/topological bridges and ties also prevent a decisive
        // result. Representative parsimony must never manufacture confidence: compare against
        // every competitor of the same circularity, not only the next ranked representative.
        rep.decisive = !rep.search_truncated
            && !rep.cover_size_limited
            && ranked
                .first()
                .is_some_and(|(_, c, _)| !c.incomplete_junction_evidence())
            && rep.margin_ln.map_or(true, |m| m >= DECISIVE_LN);
        if let Some((_,best,_))=ranked.first() {
            if let Some(fit)=&best.read_support {
                rep.margin_ln=None; // A record-support objective is not a log likelihood.
                rep.evidence_tied_candidates=ranked.iter().filter(|(_,c,_)|c.read_support.as_ref().is_some_and(|s|
                    s.joint_records_explained==fit.joint_records_explained
                    && s.thread_records_explained==fit.thread_records_explained
                    && s.adjacency_record_occurrences_explained==fit.adjacency_record_occurrences_explained)).count();
                rep.decisive=!rep.search_truncated && !rep.cover_size_limited
                    && !best.incomplete_junction_evidence() && fit.complete()
                    && rep.evidence_tied_candidates==1;
            }
        }
        if rep.read_model_budget_exhausted {rep.decisive=false;}
        if let Some(model)=&component_model {
            if !ranked.is_empty() && !rep.read_model_budget_exhausted {
                let fits:Vec<_>=ranked.iter().map(|(_,c,_)|c.read_support.as_ref().unwrap()).collect();
                let path_cover=model.path_cover(&fits,candidate_work_budget(1_000_000,comp.len(),ranked.len()));
                rep.witness_configurations=path_cover.candidate_indices.iter().map(|&i|ranked[i].1.clone()).collect();
                rep.ambiguity_configurations=path_cover.equivalent_candidate_indices.iter().map(|&i|ranked[i].1.clone()).collect();
                rep.read_path_cover=Some(path_cover);
            }
        }
        let mut it = ranked.into_iter();
        if let Some((mut walks, mut best, cover)) = it.next() {
            rep.molecules = best.molecules.len();
            // a linear molecule that leads back into itself: the circle and its tail
            let mut cut: Vec<Found> = Vec::new();
            for f in &cover {
                match if read_model.is_none() {cut_tail_phased(&g, &br, f, &phases)} else {None} {
                    Some((circle, tail, left_out)) => {
                        rep.tails_cut.push(label(&br[left_out].walk));
                        cut.push(circle);
                        // a tail of optional anchors alone is a skipped configuration
                        if candidate_walk(&br, &tail)
                            .iter()
                            .any(|o| !repeat[o.0] && !optional[o.0])
                        {
                            cut.push(tail);
                        }
                    }
                    None => cut.push(f.clone()),
                }
            }
            if !rep.tails_cut.is_empty() {
                (walks, best) = cover_report(&cut)?;
                rep.molecules = cut.len();
                rep.decisive = false;
            }
            for (walk, molecule) in walks.iter_mut().zip(best.molecules.iter_mut()) {
                if let Some((laid, added, closed)) = lay_read_copies(&g, walk, molecule.circular, &read_copies, &graph_next, &indexed_overlaps, &repeat)
                {
                    if closed {
                        molecule.circular = true;
                        rep.read_copies_added.push(format!("{} closed into a circle", label(&laid)));
                    }
                    molecule.length = spelled_length(&g, &laid, molecule.circular, &indexed_overlaps)?;
                    molecule.path = label(&laid);
                    *walk = laid;
                    rep.read_copies_added
                        .extend(added.iter().map(|(unit, n)| format!("{} x{n}", label(unit))));
                }
            }
            chosen.extend(walks.into_iter().zip(best.molecules.iter().cloned()));
            rep.best = Some(best);
        }
        rep.alternatives = it.take(p.max_alternatives).map(|(_, c, _)| c).collect();
        if let Some(prior) = prior {
            if let Some((cover, mut taken)) = prior_cover(&comp, &optional, &br, &frac, &expected, prior, component_model.as_ref()) {
                let (mut walks, mut report) = cover_report(&cover)?;
                // tandem units as often as single reads cross them, as for `best`
                for (walk, molecule) in walks.iter_mut().zip(report.molecules.iter_mut()) {
                    if let Some((laid, _, closed)) = lay_read_copies(&g, walk, molecule.circular, &read_copies, &graph_next, &indexed_overlaps, &repeat) {
                        molecule.circular |= closed;
                        molecule.length = spelled_length(&g, &laid, molecule.circular, &indexed_overlaps)?;
                        molecule.path = label(&laid);
                        *walk = laid;
                    }
                }
                // how the cover fits the target's records, as the searched covers are scored
                if let Some(model) = &component_model {
                    let circular: Vec<bool> = report.molecules.iter().map(|m| m.circular).collect();
                    let mut budget = crate::read_model::WORK_BUDGET;
                    report.read_support = model.score(&walks, &circular, &mut budget);
                }
                taken.cover = report;
                rep.prior_cover = Some(taken);
            }
        }
        comps.push(rep);
    }
    // An anchor no bridge reaches is still assembled sequence of the genome: where the reads
    // are too few to close a circle or to join pieces (Thuidium mitochondrion at ~14x: one
    // 100.7 kb unitig of a 103.1 kb genome; Eleocharis: two pieces of a 2.03 Mb genome), leaving
    // it out reported no sequence at all. It is its own linear molecule, not decisive. Shallow
    // (optional) anchors stay out, as low-frequency configurations do everywhere else, and so
    // do pieces shorter than `MIN_UNBRIDGED_BP`: short-read graphs leave stubs of a few hundred
    // bases at unresolved branches (108 bp in the Salvia and Zal-1 panels), while every genome
    // piece among the archived reports was at least 2.7 kb (HiFi unitigs are >= 1 kb anyway).
    let bridged: FxHashSet<usize> = br
        .iter()
        .flat_map(|b| [b.walk[0].0, b.walk[b.walk.len() - 1].0])
        .collect();
    // single segments that are molecules of their own, and why
    let mut lone: Vec<(usize, &'static str)> = Vec::new();
    for &a in anchors
        .iter()
        .filter(|&&a| !bridged.contains(&a) && !optional[a])
    {
        if g.seqs[a].len() < MIN_UNBRIDGED_BP {
            partial.push(describe(&Found {
                start: (a, true),
                bridges: Vec::new(),
                circular: false,
            })?);
            continue;
        }
        lone.push((a, "unbridged_sequence"));
    }
    // A repeat no chosen bridge crosses is laid down nowhere, yet its sequence is in the genome at
    // least once. Where the reads already leave the genome in linear pieces it is one more piece,
    // as an anchor no bridge reaches is (rice SAMEA115396709's mitochondrion from 100 bp reads:
    // 10.6 kb of its 316 kb in repeats of 1.1-3.1 kb at 1.5-3 times its depth, in no molecule).
    // A result of circles alone stays as it is: there an unplaced repeat is a configuration the
    // reads do not take, reported in `repeats_not_placed`. Pieces under `MIN_UNBRIDGED_BP` stay
    // out, as stubs do above.
    if !lone.is_empty() || chosen.iter().any(|(_, l)| !l.circular) {
        let laid: FxHashSet<usize> = chosen
            .iter()
            .flat_map(|(w, _)| w.iter().map(|o| o.0))
            .collect();
        for r in (0..g.names.len()).filter(|&r| repeat[r] && !laid.contains(&r)) {
            if g.seqs[r].len() >= MIN_UNBRIDGED_BP {
                lone.push((r, "unplaced_repeat"));
            }
        }
    }
    for (a, model) in lone {
        let walk = vec![(a, true)];
        let lin = Linearization {
            path: label(&walk),
            circular: false,
            length: spell(&g, &walk, false)?.len(),
            bridge_reads: Vec::new(),
            exit_fractions: Vec::new(),
            score: 0.0,
        };
        chosen.push((walk, lin.clone()));
        comps.push(ComponentReport {
            ranking_model: model,
            evidence_tied_candidates: 0,
            read_model_budget_exhausted: false,
            read_path_cover: None,
            witness_configurations: Vec::new(),
            ambiguity_configurations: Vec::new(),
            anchors: vec![g.names[a].clone()],
            depth: None,
            molecules: 1,
            path_ends: vec![g.names[a].clone()],
            min_molecules: 1,
            skipped: false,
            candidates: 1,
            search_truncated: false,
            adaptive_retry: false,
            searched_max_molecules: 1,
            cover_size_limited: false,
            read_join_cover: false,
            best: Some(CoverReport {
                score: 0.0,
                molecules: vec![lin],
                weak_bridges: 0,
                topological_bridges: 0,
                copy_deviation: 0,
                read_support: None,
            }),
            partial: None,
            margin_ln: None,
            decisive: false,
            tails_cut: Vec::new(),
            read_copies_added: Vec::new(),
            prior_cover: None,
            alternatives: Vec::new(),
        });
    }
    partial.sort_by(|(wa, a), (wb, b)| b.length.cmp(&a.length).then_with(|| wa.cmp(wb)));
    // largest component (by its longest molecule) first; unsolved ones last
    let longest = |c: &ComponentReport| {
        c.best
            .as_ref()
            .and_then(|b| b.molecules.first())
            .map_or(0, |m| m.length)
    };
    comps.sort_by(|a, b| {
        longest(b)
            .cmp(&longest(a))
            .then_with(|| a.anchors.cmp(&b.anchors))
    });
    chosen.sort_by(|(wa, a), (wb, b)| b.length.cmp(&a.length).then_with(|| wa.cmp(wb)));

    let placed = |a: usize| chosen.iter().any(|(w, _)| w.iter().any(|o| o.0 == a));
    let in_partial = |a: usize| partial.iter().any(|(w, _)| w.iter().any(|o| o.0 == a));
    let repeats_not_placed = if chosen.is_empty() {
        Vec::new()
    } else {
        (0..g.names.len())
            .filter(|&s| repeat[s] && !placed(s))
            .map(|s| g.names[s].clone())
            .collect()
    };
    if !chosen.is_empty() {
        fs::create_dir_all(out_fasta.parent().unwrap_or(Path::new(".")))?;
        fs::write(out_fasta, fasta_records(&g, &chosen, "linear", false, "")?)?;
    }
    if partial.is_empty() {
        // none in this run: one left by an earlier run would pass for this run's
        fs::remove_file(&partial_fasta).ok();
    } else {
        fs::create_dir_all(partial_fasta.parent().unwrap_or(Path::new(".")))?;
        fs::write(
            &partial_fasta,
            fasta_records(&g, &partial, "partial", true, " complete=false")?,
        )?;
    }
    let main = comps.first().filter(|c| c.best.is_some());
    let mut configuration_fasta=String::new();
    let label_index:FxHashMap<_,_>=g.names.iter().enumerate().map(|(i,n)|(n.as_str(),i)).collect();
    for (ci,component) in comps.iter().enumerate() {
        for (fi,cover) in component.witness_configurations.iter().enumerate() {
            let prefix = format!("configuration.{}.{}", ci + 1, fi + 1);
            configuration_fasta.push_str(&cover_fasta_records(&g, cover, &label_index, &prefix,
                " representation=read_path_cover physical_topology_identified=false")?);
        }
    }
    let has_configurations=!configuration_fasta.is_empty();
    if !has_configurations {fs::remove_file(&configuration_path).ok();}
    else {fs::write(&configuration_path,configuration_fasta)?;}
    if let Some(path) = out_candidate {
        write_candidate_fasta(&g, &comps, &label_index, path)?;
    }
    let out = LinearizeReport {
        ranking_model: if read_model.is_some() {crate::read_model::MODEL} else {"legacy_bridge_copy"},
        effective_ranking_model: effective_ranking_model(comps.iter().filter(|c| !c.skipped && c.candidates > 0).map(|c| c.ranking_model)),
        read_model_budget_exhausted: comps.iter().any(|c| c.read_model_budget_exhausted),
        physical_topology_identified: false,
        configurations_fasta: has_configurations.then(||configuration_path.display().to_string()),
        candidate_fasta: out_candidate.map(|path| path.display().to_string()),
        graph: gfa.display().to_string(),
        evidence: evidence.display().to_string(),
        anchors: anchors.iter().map(|&s| g.names[s].clone()).collect(),
        repeats: (0..g.names.len())
            .filter(|&s| repeat[s])
            .map(|s| g.names[s].clone())
            .collect(),
        bridges: br.len() / 2,
        topological_bridges: br.iter().filter(|b| b.topological).count() / 2,
        weak_bridges: br.iter().filter(|b| b.weak).count() / 2,
        optional_anchors: anchors
            .iter()
            .filter(|&&a| optional[a])
            .map(|&a| g.names[a].clone())
            .collect(),
        skipped_anchors: anchors
            .iter()
            .filter(|&&a| optional[a] && !placed(a) && !in_partial(a))
            .map(|&a| g.names[a].clone())
            .collect(),
        unbridged_anchors: anchors
            .iter()
            .filter(|&&s| !br.iter().any(|b| b.walk[0].0 == s))
            .map(|&s| g.names[s].clone())
            .collect(),
        local_repeat_closures: local_repeat_closures(&g, &report, &repeat, read_model.is_some(), &phases, p.min_reads, p.max_steps),
        self_loops,
        unsolved_components: comps
            .iter()
            .filter(|c| !c.skipped && c.best.is_none())
            .count(),
        molecules: chosen.len(),
        partial_fasta: (!partial.is_empty()).then(|| partial_fasta.display().to_string()),
        unplaced_segments: (0..g.names.len())
            .filter(|&s| !placed(s) && !in_partial(s))
            .map(|s| g.names[s].clone())
            .collect(),
        partial_paths: partial.iter().map(|(_, l)| l.clone()).collect(),
        candidates: comps.iter().map(|c| c.candidates).sum(),
        search_truncated: comps.iter().any(|c| c.search_truncated),
        cover_size_limited: comps.iter().any(|c| c.cover_size_limited),
        best: chosen.first().map(|(_, l)| l.clone()),
        margin_ln: main.and_then(|c| c.margin_ln),
        decisive: !chosen.is_empty() && comps.iter().all(|c| c.skipped || c.decisive),
        // the main component's alternatives that are a single molecule, as before
        alternatives: main
            .map(|c| {
                c.alternatives
                    .iter()
                    .filter(|a| a.molecules.len() == 1)
                    .map(|a| a.molecules[0].clone())
                    .collect()
            })
            .unwrap_or_default(),
        components: comps,
        repeats_not_placed,
    };
    fs::create_dir_all(out_json.parent().unwrap_or(Path::new(".")))?;
    fs::write(out_json, serde_json::to_string_pretty(&out)?)?;
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn seq(n: usize, seed: u8) -> Vec<u8> {
        (0..n)
            .map(|i| b"ACGT"[(i * 7 + seed as usize * 13 + i / 3) % 4])
            .collect()
    }

    /// A graph of segments of the given lengths and blunt links, its successor map and overlaps.
    fn blunt(
        lens: &[usize],
        links: &[(Oriented, Oriented)],
    ) -> (Graph, FxHashMap<Oriented, Vec<Oriented>>, FxHashMap<(Oriented, Oriented), usize>) {
        let g = Graph {
            names: (0..lens.len()).map(|i| format!("u{i}")).collect(),
            seqs: lens.iter().enumerate().map(|(i, &n)| seq(n, i as u8)).collect(),
            links: links.to_vec(),
            overlaps: vec![0; links.len()],
        };
        let mut next: FxHashMap<Oriented, Vec<Oriented>> = FxHashMap::default();
        for &(x, y) in links {
            next.entry(x).or_default().push(y);
            next.entry(flip(y)).or_default().push(flip(x));
        }
        let ov = overlap_index(&g);
        (g, next, ov)
    }

    /// The Haliclona shape: unique u0, tandem unit u1-u2. The circle u0+ u1- u2+ u1- laid the
    /// unit once; reads crossing u1- -> u2+ four times, 500 bp apart, need it four times.
    #[test]
    fn read_copies_lay_a_tandem_unit_as_often_as_one_read_crosses_it() {
        let (u0, u1m, u2) = ((0, true), (1, false), (2, true));
        let (g, next, ov) = blunt(&[1000, 300, 200], &[(u0, u1m), (u1m, u2), (u2, u1m), (u1m, u0)]);
        let mut need = FxHashMap::default();
        need.insert(link_key(u1m, u2), (4, 500));
        let walk = vec![u0, u1m, u2, u1m];
        let (laid, added, _) = lay_read_copies(&g, &walk, true, &need, &next, &ov, &[true; 8]).unwrap();
        assert_eq!(laid, vec![u0, u1m, u2, u1m, u2, u1m, u2, u1m, u2, u1m]);
        assert_eq!(added, vec![(vec![u2, u1m], 3)]);
        // a linear walk through the unit once, with no cycle of its own: the graph's return
        let (laid, _, _) = lay_read_copies(&g, &[u0, u1m, u2], false, &need, &next, &ov, &[true; 8]).unwrap();
        assert_eq!(laid, vec![u0, u1m, u2, u1m, u2, u1m, u2, u1m, u2]);
        // the reverse complement of the link counts as the same link
        let mut rc = FxHashMap::default();
        rc.insert(link_key(flip(u2), flip(u1m)), (2, 500));
        assert_eq!(lay_read_copies(&g, &walk, true, &rc, &next, &ov, &[true; 8]).unwrap().0.len(), 6);
        // crossings a whole circle apart (a read longer than the genome) are no tandem unit
        let mut wrap = FxHashMap::default();
        wrap.insert(link_key(u1m, u2), (2, 1500));
        assert!(lay_read_copies(&g, &walk, true, &wrap, &next, &ov, &[true; 8]).is_none());
        // already enough passes, or no return to the link in the graph: unchanged
        need.insert(link_key(u1m, u2), (1, 500));
        assert!(lay_read_copies(&g, &walk, true, &need, &next, &ov, &[true; 8]).is_none());
        let mut open = FxHashMap::default();
        open.insert(link_key(u0, u1m), (3, 500));
        let (g2, line, ov2) = blunt(&[1000, 300, 200], &[(u0, u1m), (u1m, u2)]);
        assert!(lay_read_copies(&g2, &[u0, u1m, u2], false, &open, &line, &ov2, &[true; 8]).is_none());
    }

    /// The plastid shape with its inverted repeat laid once: LSC u0, IR u1, SSC u2, the linear
    /// u0+ u1+ u2+. Reads crossing u0+ -> u1+ twice a genome (1,800 bp) apart close it into the
    /// circle u0+ u1+ u2+ u1-, not a second genome laid on a line.
    #[test]
    fn read_copies_a_genome_apart_close_a_linear_walk_into_its_circle() {
        let (u0, u1, u2) = ((0, true), (1, true), (2, true));
        let (g, next, ov) = blunt(&[1000, 300, 200], &[(u0, u1), (u1, u2), (u2, flip(u1)), (flip(u1), u0)]);
        let mut need = FxHashMap::default();
        need.insert(link_key(u0, u1), (2, 1800));
        let (laid, added, closed) = lay_read_copies(&g, &[u0, u1, u2], false, &need, &next, &ov, &[true; 8]).unwrap();
        assert!(closed);
        assert_eq!(laid, vec![u0, u1, u2, flip(u1)]);
        assert!(added.is_empty());
        assert_eq!(spelled_length(&g, &laid, true, &ov).unwrap(), 1800);
        // crossings a unit apart that is not the whole walk stay a tandem copy on the line
        need.insert(link_key(u0, u1), (2, 900));
        assert!(lay_read_copies(&g, &[u0, u1, u2], false, &need, &next, &ov, &[true; 8]).is_none());
    }

    /// u0 u1 u2 linear; reads cross u1+ -> u2+ twice, 5.5 kb apart, and the graph returns
    /// u2 -> u3 -> u1. With u3 (5 kb) single copy the unit u2 u3 u1 is mostly single copy: a
    /// subcircle read round, nothing is laid. With every segment a repeat it is a tandem copy.
    #[test]
    fn a_unit_of_mostly_single_copy_sequence_is_a_subcircle_not_a_tandem_copy() {
        let (u0, u1, u2, u3) = ((0, true), (1, true), (2, true), (3, true));
        let (g, next, ov) = blunt(&[1000, 300, 200, 5000], &[(u0, u1), (u1, u2), (u2, u3), (u3, u1)]);
        let mut need = FxHashMap::default();
        need.insert(link_key(u1, u2), (2, 5500));
        let single = [false, false, false, false];
        assert!(lay_read_copies(&g, &[u0, u1, u2], false, &need, &next, &ov, &single).is_none());
        let repeats = [false, true, true, true];
        let (laid, added, closed) = lay_read_copies(&g, &[u0, u1, u2], false, &need, &next, &ov, &repeats).unwrap();
        assert!(!closed);
        assert_eq!(laid, vec![u0, u1, u2, u3, u1, u2]);
        assert_eq!(added.len(), 1);
    }

    /// A self-loop's unit is its segment; a circle that is its unit alone stays as it is.
    #[test]
    fn read_copies_lay_a_self_loop_but_not_a_circle_of_its_unit() {
        let (a, r, b) = ((0, true), (1, true), (2, true));
        let (g, next, ov) = blunt(&[1000, 100, 1000], &[(a, r), (r, r), (r, b), (b, a)]);
        let mut need = FxHashMap::default();
        need.insert(link_key(r, r), (3, 100));
        let (laid, added, _) = lay_read_copies(&g, &[a, r, r, b], true, &need, &next, &ov, &[true; 8]).unwrap();
        assert_eq!(laid, vec![a, r, r, r, r, b]);
        assert_eq!(added, vec![(vec![r], 2)]);
        // a link the walk does not pass is not laid
        assert!(lay_read_copies(&g, &[a, r, b], true, &need, &next, &ov, &[true; 8]).is_none());
        assert!(lay_read_copies(&g, &[r], true, &need, &next, &ov, &[true; 8]).is_none());
    }

    /// Anchors A, B, C, D around one repeat R: A->R, C->R, R->B, R->D, plus B->C and D->C.
    fn graph() -> Graph {
        Graph {
            names: ["A", "B", "C", "D", "R"].map(String::from).to_vec(),
            seqs: (0..5).map(|i| seq(40, i)).collect(),
            links: vec![
                ((0, true), (4, true)),
                ((2, true), (4, true)),
                ((4, true), (1, true)),
                ((4, true), (3, true)),
                ((1, true), (2, true)),
                ((3, true), (2, true)),
            ],
            overlaps: vec![0; 6],
        }
    }

    fn evidence(phased: &[(&str, u64)], links: &[(&str, &str, u64)]) -> Value {
        json!({
            "repeat_segments": ["R"],
            "phased_paths": phased.iter().map(|(p, n)| json!({"path": p, "reads": n})).collect::<Vec<_>>(),
            "links": links.iter().map(|(a, b, n)| json!({"from": a, "to": b, "reads": n})).collect::<Vec<_>>(),
        })
    }

    type Solution = (Vec<Oriented>, bool, Vec<u64>, Vec<f64>);

    #[test]
    fn local_adjacency_joins_anchors_but_never_phases_a_repeat() {
        let g = graph();
        let repeat = vec![false, false, false, false, true];
        let ev = json!({
            "phased_paths": [],
            "links": [
                {"from": "B+", "to": "C+", "reads": 0, "adjacency_reads": 15},
                {"from": "A+", "to": "R+", "reads": 0, "adjacency_reads": 30},
                {"from": "R+", "to": "B+", "reads": 0, "adjacency_reads": 30},
            ],
        });
        let (br, _) = bridges(&g, &ev, &repeat, 2).unwrap();
        let bc = br
            .iter()
            .find(|b| b.walk == vec![(1, true), (2, true)])
            .unwrap();
        assert_eq!(bc.reads, 15);
        assert!(!bc.topological && !bc.weak);
        let through_repeat: Vec<_> = br.iter().filter(|b| b.walk.contains(&(4, true))).collect();
        assert!(!through_repeat.is_empty());
        assert!(through_repeat.iter().all(|b| b.topological));
    }

    fn solve(g: &Graph, ev: &Value) -> Vec<Solution> {
        let repeat = vec![false, false, false, false, true];
        let (br, _) = bridges(g, ev, &repeat, 1).unwrap();
        let (found, truncated) = search(5, &[0, 1, 2, 3], &[false; 5], &br, 1, 1_000_000);
        assert!(!truncated);
        let frac = exit_fractions(&br);
        let mut seen = BTreeSet::new();
        found
            .iter()
            .map(|cover| &cover[0])
            .filter_map(|f| {
                let w = candidate_walk(&br, f);
                seen.insert(canonical_key(&w, f.circular)).then(|| {
                    (
                        w,
                        f.circular,
                        f.bridges.iter().map(|&i| br[i].reads).collect(),
                        f.bridges.iter().map(|&i| frac[i]).collect(),
                    )
                })
            })
            .collect()
    }

    #[test]
    fn lays_the_repeat_once_per_crossing() {
        let g = graph();
        let ev = evidence(&[("A+ R+ B+", 30), ("C+ R+ D+", 20)], &[("B+", "C+", 10)]);
        let sols = solve(&g, &ev);
        assert_eq!(
            sols.len(),
            1,
            "a path and its reverse complement are one candidate"
        );
        let (walk, circular, reads, _) = &sols[0];
        let fw = vec![
            (0, true),
            (4, true),
            (1, true),
            (2, true),
            (4, true),
            (3, true),
        ];
        assert!(*walk == fw || *walk == rc_walk(&fw));
        assert!(!circular);
        assert_eq!(spell(&g, &fw, false).unwrap().len(), 6 * 40);
        let mut r = reads.clone();
        r.sort();
        assert_eq!(r, vec![10, 20, 30]);
    }

    #[test]
    fn closes_the_circle_without_repeating_the_first_anchor() {
        let mut g = graph();
        g.links.push(((3, true), (0, true)));
        g.overlaps.push(0);
        let ev = evidence(
            &[("A+ R+ B+", 30), ("C+ R+ D+", 20)],
            &[("B+", "C+", 10), ("D+", "A+", 15)],
        );
        let sols = solve(&g, &ev);
        assert_eq!(sols.len(), 1, "rotations of one circle are one candidate");
        let (walk, circular, reads, _) = &sols[0];
        assert!(*circular);
        assert_eq!(walk.len(), 6);
        assert_eq!(reads.len(), 4, "the closing bridge is scored too");
        assert_eq!(spell(&g, walk, true).unwrap().len(), 6 * 40);
    }

    #[test]
    fn repeats_or_anchors_without_bridges_leave_no_path() {
        let g = graph();
        let ev = evidence(&[("A+ R+ B+", 30), ("C+ R+ D+", 20)], &[]);
        assert!(solve(&g, &ev).is_empty());
    }

    #[test]
    fn overlaps_are_written_once() {
        let x = seq(50, 1);
        let y: Vec<u8> = [&x[40..], &seq(30, 2)[..]].concat();
        let g = Graph {
            names: vec!["X".into(), "Y".into()],
            seqs: vec![x.clone(), y.clone()],
            links: vec![((0, true), (1, true)), ((1, true), (0, true))],
            overlaps: vec![10, 0],
        };
        let lin = spell(&g, &[(0, true), (1, true)], false).unwrap();
        assert_eq!(lin, [&x[..], &y[10..]].concat());
        // Reverse traversal uses the same link's overlap.
        let back = spell(&g, &[(1, false), (0, false)], false).unwrap();
        assert_eq!(back.len(), lin.len());
        // Circular: Y -> X has no overlap, so nothing is trimmed from the end.
        assert_eq!(spell(&g, &[(0, true), (1, true)], true).unwrap(), lin);
    }

    #[test]
    fn the_better_supported_configuration_ranks_first() {
        let g = graph();
        // Strong: A R D, C R B. Weak: A R B, C R D (lexically first, so ties would pick it).
        let ev = evidence(
            &[
                ("A+ R+ D+", 50),
                ("C+ R+ B+", 40),
                ("A+ R+ B+", 5),
                ("C+ R+ D+", 4),
            ],
            &[("B+", "C+", 20), ("D+", "C+", 20)],
        );
        let sols = solve(&g, &ev);
        assert_eq!(sols.len(), 2);
        let best = sols
            .iter()
            .max_by(|a, b| path_score(&a.3).total_cmp(&path_score(&b.3)))
            .unwrap();
        let mut r = best.2.clone();
        r.sort();
        assert_eq!(r, vec![20, 40, 50]);
    }

    #[test]
    fn exit_fractions_are_depth_free_and_direction_free() {
        let g = graph();
        let repeat = vec![false, false, false, false, true];
        // Same configuration at 10x the depth gives the same fractions.
        let f = |scale: u64| {
            let ev = evidence(
                &[("A+ R+ D+", 3 * scale), ("A+ R+ B+", scale)],
                &[("B+", "C+", 7 * scale)],
            );
            // read-supported bridges only (C+ is a dead end: its graph pairings join too)
            let (br, _) = bridges(&g, &ev, &repeat, 1).unwrap();
            let br: Vec<Bridge> = br.into_iter().filter(|b| !b.weak).collect();
            let fr = exit_fractions(&br);
            let at = |w: &[Oriented]| fr[br.iter().position(|b| b.walk == w).unwrap()];
            (
                at(&[(0, true), (4, true), (3, true)]),
                at(&[(3, false), (4, false), (0, false)]),
                at(&[(0, true), (4, true), (1, true)]),
            )
        };
        let (ard, dra, arb) = f(1);
        // A+ leaves 3:1 towards D; D- and B- each have a single exit.
        assert!((ard - 0.75f64.sqrt()).abs() < 1e-12 && (arb - 0.5).abs() < 1e-12);
        assert_eq!(
            ard, dra,
            "a bridge and its reverse complement share one score"
        );
        assert_eq!(f(10), f(1));
    }

    #[test]
    fn a_weakly_bridged_anchor_is_optional_and_skipped() {
        // A R B C R D as before, plus E hanging off D by a weak link: no path can take E in the
        // middle, so requiring it would leave no linearization at all.
        let mut g = graph();
        g.names.push("E".into());
        g.seqs.push(seq(40, 9));
        g.links.push(((3, true), (5, true)));
        g.overlaps.push(0);
        g.links.push(((5, true), (0, true)));
        g.overlaps.push(0);
        let ev = evidence(
            &[("A+ R+ B+", 30), ("C+ R+ D+", 20)],
            &[("B+", "C+", 10), ("D+", "E+", 1), ("E+", "A+", 1)],
        );
        let repeat = vec![false, false, false, false, true, false];
        let (br, _) = bridges(&g, &ev, &repeat, 1).unwrap();
        let anchors = [0, 1, 2, 3, 5];
        let (all, _) = search(6, &anchors, &[false; 6], &br, 1, 1_000_000);
        assert!(
            all.iter().all(|c| c[0].circular),
            "only the closed circle through E visits all"
        );
        let mut optional = [false; 6];
        optional[5] = true;
        let (found, _) = search(6, &anchors, &optional, &br, 1, 1_000_000);
        assert!(
            found.iter().map(|c| &c[0]).any(|f| {
                let w = candidate_walk(&br, f);
                !w.iter().any(|o| o.0 == 5) && !f.circular
            }),
            "the linear path without E is now a candidate"
        );
    }

    /// Anchors only, joined by the given links (all 0M).
    fn anchors_graph(names: &[&str], links: &[(Oriented, Oriented)]) -> Graph {
        Graph {
            names: names.iter().map(|n| n.to_string()).collect(),
            seqs: (0..names.len()).map(|i| seq(40, i as u8)).collect(),
            links: links.to_vec(),
            overlaps: vec![0; links.len()],
        }
    }

    #[test]
    fn separate_chromosomes_unring_apart_and_a_self_link_closes_a_circle() {
        // A <-> B is one circle; C with a self-link is a second one (a one-segment chromosome).
        let g = anchors_graph(
            &["A", "B", "C"],
            &[
                ((0, true), (1, true)),
                ((1, true), (0, true)),
                ((2, true), (2, true)),
            ],
        );
        let ev = json!({
            "repeat_segments": [],
            "phased_paths": [],
            "links": [
                {"from": "A+", "to": "B+", "reads": 10},
                {"from": "B+", "to": "A+", "reads": 8},
                {"from": "C+", "to": "C+", "reads": 12},
            ],
        });
        let (br, loops) = bridges(&g, &ev, &[false; 3], 1).unwrap();
        assert_eq!(loops, vec!["C+ C+".to_string()]);
        let comps = components(3, &br);
        assert_eq!(comps, vec![vec![0, 1], vec![2]]);
        for comp in &comps {
            let (found, _) = search(3, comp, &[false; 3], &br, 1, 1_000);
            assert!(!found.is_empty());
            assert!(found.iter().all(|c| c.len() == 1 && c[0].circular));
        }
        let (c, _) = search(3, &[2], &[false; 3], &br, 1, 1_000);
        let w = candidate_walk(&br, &c[0][0]);
        assert_eq!(w.len(), 1, "one segment, closed by its own link");
        assert_eq!(spell(&g, &w, true).unwrap().len(), 40);
    }

    #[test]
    fn a_component_no_single_path_covers_splits_into_the_fewest_molecules() {
        // Star: X -> A, X -> B, X -> C. Any molecule holds X and one leaf at most, so the
        // fewest molecules is three (X with one leaf, and the other two leaves alone).
        let g = anchors_graph(
            &["X", "A", "B", "C"],
            &[
                ((0, true), (1, true)),
                ((0, true), (2, true)),
                ((0, true), (3, true)),
            ],
        );
        let ev = json!({
            "repeat_segments": [],
            "phased_paths": [],
            "links": [
                {"from": "X+", "to": "A+", "reads": 10},
                {"from": "X+", "to": "B+", "reads": 10},
                {"from": "X+", "to": "C+", "reads": 10},
            ],
        });
        let (br, _) = bridges(&g, &ev, &[false; 4], 1).unwrap();
        let comp = [0, 1, 2, 3];
        for m in 1..=2 {
            assert!(search(4, &comp, &[false; 4], &br, m, 100_000).0.is_empty());
        }
        let (found, truncated) = search(4, &comp, &[false; 4], &br, 3, 100_000);
        assert!(!truncated && !found.is_empty());
        for cover in &found {
            let mut seen: Vec<usize> = cover
                .iter()
                .flat_map(|f| candidate_walk(&br, f))
                .map(|o| o.0)
                .collect();
            seen.sort_unstable();
            assert_eq!(seen, vec![0, 1, 2, 3], "every anchor exactly once");
            assert!(cover.iter().all(|f| !f.circular));
        }
    }

    #[test]
    fn circles_rank_before_likelier_linear_covers() {
        let mol = |circular, score| Linearization {
            path: String::new(),
            circular,
            length: 1,
            bridge_reads: vec![],
            exit_fractions: vec![],
            score,
        };
        let cover = |m: Vec<Linearization>| CoverReport {
            score: m.iter().map(|l| l.score).sum(),
            molecules: m,
            weak_bridges: 0,
            topological_bridges: 0,
            copy_deviation: 0,
            read_support: None,
        };
        let circles = cover(vec![mol(true, -1.0), mol(true, -1.0)]);
        let open = cover(vec![mol(true, -0.1), mol(false, -0.1)]);
        let better_circles = cover(vec![mol(true, -0.5), mol(true, -0.5)]);
        assert_eq!(cover_order(&circles, &open), Ordering::Less);
        assert_eq!(cover_order(&better_circles, &circles), Ordering::Less);
        // a likelier cover that lays a repeat down once too often ranks after one that does not
        let mut off_by_one = better_circles.clone();
        off_by_one.copy_deviation = 1;
        assert_eq!(cover_order(&circles, &off_by_one), Ordering::Less);
    }

    /// Plastid LSC (u0), SSC (u1) and one IR copy (u2), as `ovasm assemble graph` collapses them,
    /// with the junction reads of the Antitrichia holdout: LSC-IR-SSC read twice, the pairing
    /// that closes the circle and its flip-flop rival once each.
    fn run_plastid(phased: Value, dir_name: &str) -> LinearizeReport {
        let dir = std::env::temp_dir().join(format!("{dir_name}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let s = |n, seed| String::from_utf8(seq(n, seed)).unwrap();
        let gfa = format!(
            "H\tVN:Z:1.0\nS\tu0\t{}\nS\tu1\t{}\nS\tu2\t{}\n\
             L\tu0\t-\tu2\t+\t0M\nL\tu0\t+\tu2\t+\t0M\nL\tu1\t-\tu2\t-\t0M\nL\tu1\t+\tu2\t-\t0M\n",
            s(400, 1),
            s(100, 2),
            s(60, 3)
        );
        std::fs::write(dir.join("g.gfa"), gfa).unwrap();
        let ev = json!({
            "repeat_segments": ["u2"],
            "segment_depths": {"u0": 32.8, "u1": 39.1, "u2": 104.5},
            "median_depth": 32.8,
            "phased_paths": phased,
            "links": [
                {"from": "u0-", "to": "u2+", "reads": 82},
                {"from": "u0+", "to": "u2+", "reads": 84},
                {"from": "u1-", "to": "u2-", "reads": 90},
                {"from": "u1+", "to": "u2-", "reads": 80},
            ],
        });
        std::fs::write(dir.join("e.json"), ev.to_string()).unwrap();
        let params = LinearizeParams {
            min_reads: 2,
            max_steps: 1_000_000, adaptive_steps: false,
            max_molecules: 4,
            max_alternatives: 10,
            optional_below: 0.5,
        };
        let r = run(
            &dir.join("g.gfa"),
            &dir.join("e.json"),
            params,
            &dir.join("m.fa"),
            None,
            &dir.join("l.json"),
        )
        .unwrap();
        std::fs::remove_dir_all(&dir).ok();
        r
    }

    /// `run` on a GFA and an evidence report written to a scratch directory.
    fn run_on(gfa: &str, ev: &Value, dir_name: &str) -> LinearizeReport {
        run_files(gfa, ev, dir_name, 4).0
    }

    #[test]
    fn automatic_cover_envelope_includes_five_path_end_pairs_and_keeps_explicit_limits() {
        let mut gfa = String::from("H\tVN:Z:1.0\nS\tR\tACGTACGTACGTACGTACGT\n");
        let mut links = Vec::new();
        let mut paths = Vec::new();
        for i in 0..5 {
            gfa.push_str(&format!(
                "S\tA{i}\t{}\nS\tB{i}\t{}\nL\tA{i}\t+\tR\t+\t0M\nL\tR\t+\tB{i}\t+\t0M\n",
                String::from_utf8(seq(100, i + 1)).unwrap(),
                String::from_utf8(seq(100, i + 10)).unwrap()
            ));
            links.push(json!({"from":format!("A{i}+"),"to":"R+","reads":50}));
            links.push(json!({"from":"R+","to":format!("B{i}+"),"reads":50}));
            for j in 0..5 {
                paths.push(json!({"path":format!("A{i}+ R+ B{j}+"),"reads":10}));
            }
        }
        let ev = json!({"repeat_segments":["R"],"links":links,"phased_paths":paths});
        let explicit = run_files(&gfa, &ev, "ovasm-explicit-cover-four", 4).0;
        assert_eq!(explicit.components[0].min_molecules, 5);
        assert_eq!(explicit.unsolved_components, 1);
        let automatic = run_files(&gfa, &ev, "ovasm-auto-cover-five", 0).0;
        assert_eq!(automatic.components[0].searched_max_molecules, 5);
        assert_eq!(automatic.unsolved_components, 0);
        assert_eq!(automatic.molecules, 5);
        assert!(automatic.cover_size_limited && !automatic.decisive);
    }

    #[test]
    fn an_automatic_envelope_without_a_cover_widens_until_one_holds() {
        // Three leaves into H1, H1 to H2, three leaves out of H2: six path ends make three
        // molecules at least and an automatic envelope of four, but a cover needs five (one
        // leaf-H1-H2-leaf molecule, the four other leaves alone).
        let mut gfa = String::from("H\tVN:Z:1.0\n");
        let mut links = vec![json!({"from":"H1+","to":"H2+","reads":50})];
        gfa.push_str(&format!(
            "S\tH1\t{}\nS\tH2\t{}\nL\tH1\t+\tH2\t+\t0M\n",
            String::from_utf8(seq(100, 31)).unwrap(),
            String::from_utf8(seq(100, 32)).unwrap()
        ));
        for (i, leaf) in ["A", "B", "C"].iter().enumerate() {
            gfa.push_str(&format!(
                "S\t{leaf}1\t{}\nS\t{leaf}2\t{}\nL\t{leaf}1\t+\tH1\t+\t0M\nL\tH2\t+\t{leaf}2\t+\t0M\n",
                String::from_utf8(seq(100, 40 + i as u8)).unwrap(),
                String::from_utf8(seq(100, 50 + i as u8)).unwrap()
            ));
            links.push(json!({"from":format!("{leaf}1+"),"to":"H1+","reads":50}));
            links.push(json!({"from":"H2+","to":format!("{leaf}2+"),"reads":50}));
        }
        let ev = json!({"repeat_segments":[],"links":links,"phased_paths":[]});
        let explicit = run_files(&gfa, &ev, "ovasm-explicit-cover-widen", 4).0;
        assert_eq!(explicit.components[0].min_molecules, 3);
        assert_eq!(explicit.unsolved_components, 1, "an explicit limit stays");
        let automatic = run_files(&gfa, &ev, "ovasm-auto-cover-widen", 0).0;
        assert_eq!(automatic.unsolved_components, 0);
        assert_eq!(automatic.molecules, 5);
        assert_eq!(automatic.components[0].searched_max_molecules, 5);
        assert!(automatic.cover_size_limited && !automatic.decisive);
        assert!(automatic.partial_paths.is_empty());
    }

    /// As `run_on`, with the molecule FASTA and the partial FASTA as written (`None`: no such
    /// file). A partial FASTA of an earlier run is in place beforehand.
    fn run_files(
        gfa: &str,
        ev: &Value,
        dir_name: &str,
        max_molecules: usize,
    ) -> (LinearizeReport, Option<String>, Option<String>) {
        let dir = std::env::temp_dir().join(format!("{dir_name}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("g.gfa"), gfa).unwrap();
        std::fs::write(dir.join("e.json"), ev.to_string()).unwrap();
        std::fs::write(dir.join("m.partial.fa"), ">partial.1 stale\nACGT\n").unwrap();
        let params = LinearizeParams {
            min_reads: 2,
            max_steps: 1_000_000, adaptive_steps: false,
            max_molecules,
            max_alternatives: 10,
            optional_below: 0.5,
        };
        let r = run(
            &dir.join("g.gfa"),
            &dir.join("e.json"),
            params,
            &dir.join("m.fa"),
            None,
            &dir.join("l.json"),
        )
        .unwrap();
        let read = |name: &str| std::fs::read_to_string(dir.join(name)).ok();
        let files = (read("m.fa"), read("m.partial.fa"));
        std::fs::remove_dir_all(&dir).ok();
        (r, files.0, files.1)
    }

    /// Segments `(name, length)` with distinct sequences, and 0M links `(from, to)`.
    fn gfa_text(segs: &[(&str, usize)], links: &[(&str, &str)]) -> String {
        let mut out = String::from("H\tVN:Z:1.0\n");
        for (i, (name, len)) in segs.iter().enumerate() {
            let s = String::from_utf8(seq(*len, i as u8 + 1)).unwrap();
            out.push_str(&format!("S\t{name}\t{s}\n"));
        }
        for (a, b) in links {
            let (an, ao) = a.split_at(a.len() - 1);
            let (bn, bo) = b.split_at(b.len() - 1);
            out.push_str(&format!("L\t{an}\t{ao}\t{bn}\t{bo}\t0M\n"));
        }
        out
    }

    fn link_reads(links: &[(&str, &str, u64)]) -> Vec<Value> {
        links
            .iter()
            .map(|(a, b, n)| json!({"from": a, "to": b, "reads": n}))
            .collect()
    }

    fn names(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    /// Segments of every solved molecule.
    fn placed_segments(r: &LinearizeReport) -> BTreeSet<String> {
        r.components
            .iter()
            .filter_map(|c| c.best.as_ref())
            .flat_map(|b| &b.molecules)
            .flat_map(|m| m.path.split_whitespace())
            .map(|l| l[..l.len() - 1].to_string())
            .collect()
    }

    #[test]
    fn an_anchor_as_deep_as_its_own_chromosome_is_not_optional() {
        // Trifolium dubium (frozen seeds): the largest chromosome is one circular unitig at
        // 44.7x and holds the graph's median; the rest is one component at 32.9, 22.3 and
        // 23.5x. u2 at 0.499 of the graph's median is no rarer than its neighbours.
        let gfa = gfa_text(
            &[("u0", 500), ("u1", 1800), ("u2", 230), ("u3", 190)],
            &[("u0-", "u3+"), ("u1-", "u1-"), ("u2-", "u3+")],
        );
        let ev = json!({
            "repeat_segments": [],
            "segment_depths": {"u0": 32.9, "u1": 44.7, "u2": 22.3, "u3": 23.5},
            "median_depth": 44.7,
            "phased_paths": [],
            "links": link_reads(&[("u0-", "u3+", 29), ("u1-", "u1-", 94), ("u2-", "u3+", 49)]),
        });
        let r = run_on(&gfa, &ev, "ovasm-own-depth");
        assert!(r.optional_anchors.is_empty(), "{:?}", r.optional_anchors);
        assert!(r.skipped_anchors.is_empty());
        let all: BTreeSet<String> = names(&["u0", "u1", "u2", "u3"]).into_iter().collect();
        assert_eq!(placed_segments(&r), all);
        let small = r.components.iter().find(|c| c.anchors.len() == 3).unwrap();
        assert_eq!(small.depth, Some(32.9));
        assert_eq!(small.molecules, 2, "u3 is entered from u0 and from u2");
        assert_eq!(r.unsolved_components, 0);
    }

    #[test]
    fn depth_within_a_component_still_marks_the_rare_anchor_and_the_rare_component() {
        // Main circle A B at 20x. B also leaves by V (4x, a rare arrangement inside the main
        // component): optional and skipped. X Y is a separate component at 7 and 8x, shallow
        // as a whole (0.35-0.4 of the median): left out, as those of Bromus sterilis are.
        let gfa = gfa_text(
            &[("A", 900), ("B", 700), ("V", 120), ("X", 200), ("Y", 150)],
            &[("A+", "B+"), ("B+", "A+"), ("B+", "V+"), ("X+", "Y+")],
        );
        let ev = json!({
            "repeat_segments": [],
            "segment_depths": {"A": 20.0, "B": 21.0, "V": 4.0, "X": 7.0, "Y": 8.0},
            "median_depth": 20.0,
            "phased_paths": [],
            "links": link_reads(&[
                ("A+", "B+", 20),
                ("B+", "A+", 18),
                ("B+", "V+", 3),
                ("X+", "Y+", 6),
            ]),
        });
        let r = run_on(&gfa, &ev, "ovasm-rare-parts");
        assert_eq!(r.optional_anchors, names(&["V", "X", "Y"]));
        assert_eq!(r.skipped_anchors, names(&["V", "X", "Y"]));
        assert_eq!(r.molecules, 1);
        let best = r.best.as_ref().unwrap();
        assert!(best.circular && best.length == 1600, "{}", best.path);
        let rare = r
            .components
            .iter()
            .find(|c| c.anchors == names(&["X", "Y"]))
            .unwrap();
        assert!(rare.skipped && rare.depth == Some(7.0));
        assert_eq!(r.unsolved_components, 0);
    }

    #[test]
    fn a_path_that_leads_back_into_itself_is_a_circle_and_its_tail() {
        // Luzula sylvatica: T leads into X, and the path X Y ends with a bridge back to X.
        // One molecule covers all three only as the linear T X Y; reported are the circle X Y
        // and the tail T.
        let gfa = gfa_text(
            &[("T", 100), ("X", 600), ("Y", 500)],
            &[("T+", "X+"), ("X+", "Y+"), ("Y+", "X+")],
        );
        let ev = json!({
            "repeat_segments": [],
            "segment_depths": {"T": 18.8, "X": 27.0, "Y": 27.0},
            "median_depth": 27.0,
            "phased_paths": [],
            "links": link_reads(&[("T+", "X+", 32), ("X+", "Y+", 40), ("Y+", "X+", 30)]),
        });
        let r = run_on(&gfa, &ev, "ovasm-rho");
        let c = &r.components[0];
        assert_eq!(c.tails_cut, names(&["T+ X+"]));
        assert_eq!((c.molecules, r.molecules, r.unsolved_components), (2, 2, 0));
        let mols = &c.best.as_ref().unwrap().molecules;
        assert!(
            mols[0].circular && mols[0].length == 1100,
            "{}",
            mols[0].path
        );
        assert_eq!(mols[0].bridge_reads, vec![40, 30]);
        assert!(
            !mols[1].circular && mols[1].path == "T+",
            "{}",
            mols[1].path
        );
        assert!(!r.decisive, "reads support the bridge left out as well");
        let all: BTreeSet<String> = names(&["T", "X", "Y"]).into_iter().collect();
        assert_eq!(placed_segments(&r), all);

        // the same at the path's other end: Y leads on to T
        let gfa = gfa_text(
            &[("T", 100), ("X", 600), ("Y", 500)],
            &[("X+", "Y+"), ("Y+", "X+"), ("Y+", "T+")],
        );
        let ev = json!({
            "repeat_segments": [],
            "phased_paths": [],
            "links": link_reads(&[("X+", "Y+", 40), ("Y+", "X+", 30), ("Y+", "T+", 32)]),
        });
        let r = run_on(&gfa, &ev, "ovasm-rho-end");
        let c = &r.components[0];
        assert_eq!(c.tails_cut.len(), 1);
        let mols = &c.best.as_ref().unwrap().molecules;
        assert!(mols[0].circular && mols[0].length == 1100 && mols[1].length == 100);
    }

    #[test]
    fn a_long_tail_a_small_loop_or_a_weak_closure_is_not_cut() {
        // A leads into B, and B closes on itself: A is nine times the loop, so the path stays
        // whole.
        let gfa = gfa_text(&[("A", 900), ("B", 100)], &[("A+", "B+"), ("B+", "B+")]);
        let ev = json!({
            "repeat_segments": [],
            "phased_paths": [],
            "links": link_reads(&[("A+", "B+", 20), ("B+", "B+", 15)]),
        });
        let r = run_on(&gfa, &ev, "ovasm-rho-small");
        let c = &r.components[0];
        assert!(c.tails_cut.is_empty());
        let best = r.best.as_ref().unwrap();
        assert!(!best.circular && best.length == 1000 && r.molecules == 1);

        // T X Y with a tail a quarter of the circle: a part of the molecule, not a side branch
        let gfa = gfa_text(
            &[("T", 275), ("X", 600), ("Y", 500)],
            &[("T+", "X+"), ("X+", "Y+"), ("Y+", "X+")],
        );
        let ev = json!({
            "repeat_segments": [],
            "phased_paths": [],
            "links": link_reads(&[("T+", "X+", 32), ("X+", "Y+", 40), ("Y+", "X+", 30)]),
        });
        let r = run_on(&gfa, &ev, "ovasm-rho-long-tail");
        assert!(r.components[0].tails_cut.is_empty());
        assert_eq!(r.molecules, 1);
        assert_eq!(r.best.as_ref().unwrap().length, 1375);

        // the short tail again, but the bridge back to X was read once (below min_reads):
        // weak, kept only because Y has no other exit, and no ground for cutting T off.
        let gfa = gfa_text(
            &[("T", 100), ("X", 600), ("Y", 500)],
            &[("T+", "X+"), ("X+", "Y+"), ("Y+", "X+")],
        );
        let ev = json!({
            "repeat_segments": [],
            "phased_paths": [],
            "links": link_reads(&[("T+", "X+", 32), ("X+", "Y+", 40), ("Y+", "X+", 1)]),
        });
        let r = run_on(&gfa, &ev, "ovasm-rho-weak");
        assert!(r.components[0].tails_cut.is_empty());
        assert_eq!(r.molecules, 1);
        assert_eq!(r.best.as_ref().unwrap().length, 1200);
    }

    #[test]
    fn component_reachability_is_directed_and_fails_closed_on_external_endpoints() {
        let make = |walk| Bridge { walk, reads: 3, topological: false, weak: false };
        let br = vec![make(vec![(0, true), (1, true)]), make(vec![(1, true), (2, false)])];
        let mut out = FxHashMap::default();
        for (i, b) in br.iter().enumerate() { out.entry(b.walk[0]).or_insert_with(Vec::new).push(i); }
        let (idx, r) = component_reachability(&[0, 1, 2], &br, &out).unwrap();
        assert!(r[idx[&(0, true)]][idx[&(2, false)]]);
        assert!(!r[idx[&(2, false)]][idx[&(0, true)]]);
        assert!(!r[idx[&(0, false)]][idx[&(2, false)]]);
        assert!(component_reachability(&[0, 1], &br, &out).is_none());
        assert!(component_reachability(&(0..257).collect::<Vec<_>>(), &br, &out).is_none());
    }

    #[test]
    fn a_component_without_a_cover_comes_out_as_partial_paths_and_stays_unsolved() {
        // X leaves by A (30 reads) or B (10); A leads on to C, C to the shallow O. One molecule
        // cannot hold both A and B. U has no link at all. With one molecule allowed there is no
        // cover: the read bridges still join X A C, and B and U are reported alone.
        let gfa = gfa_text(
            &[
                ("X", 300),
                ("A", 200),
                ("B", 150),
                ("C", 400),
                ("O", 90),
                ("U", 120),
            ],
            &[("X+", "A+"), ("X+", "B+"), ("A+", "C+"), ("C+", "O+")],
        );
        let ev = json!({
            "repeat_segments": [],
            "segment_depths": {"X": 20.0, "A": 21.0, "B": 19.0, "C": 20.0, "O": 2.0, "U": 20.0},
            "median_depth": 20.0,
            "phased_paths": [],
            "links": link_reads(&[
                ("X+", "A+", 30),
                ("X+", "B+", 10),
                ("A+", "C+", 20),
                ("C+", "O+", 3),
            ]),
        });
        let (r, molecules, partial) = run_files(&gfa, &ev, "ovasm-partial", 1);
        assert_eq!((r.molecules, r.unsolved_components), (0, 1));
        assert!(r.best.is_none() && molecules.is_none() && !r.decisive);
        let c = &r.components[0];
        assert!(c.best.is_none() && c.molecules == 0 && !c.skipped);
        let paths: Vec<&str> = r.partial_paths.iter().map(|m| m.path.as_str()).collect();
        assert_eq!(paths, vec!["X+ A+ C+", "B+", "U+"]);
        assert_eq!(r.partial_paths[0].bridge_reads, vec![30, 20]);
        assert!(r.partial_paths.iter().all(|m| !m.circular));
        // the component's own paths; U is in no component
        assert_eq!(c.partial.as_ref().unwrap().molecules.len(), 2);
        // O is optional: not at a path's end, and reported as left out
        assert_eq!(r.skipped_anchors, names(&["O"]));
        assert_eq!(r.unplaced_segments, names(&["O"]));
        let fa = partial.unwrap();
        let headers: Vec<&str> = fa.lines().filter(|l| l.starts_with('>')).collect();
        assert_eq!(
            headers,
            vec![
                ">partial.1 circular=false length=900 path=X+,A+,C+ complete=false",
                ">partial.2 circular=false length=150 path=B+ complete=false",
                ">partial.3 circular=false length=120 path=U+ complete=false",
            ]
        );
        assert!(r.partial_fasta.as_ref().unwrap().ends_with("m.partial.fa"));

        // with two molecules there is a cover: nothing partial but U, and the component is
        // solved
        let (r, molecules, partial) = run_files(&gfa, &ev, "ovasm-partial-solved", 2);
        assert_eq!((r.molecules, r.unsolved_components), (2, 0));
        assert!(r.components[0].partial.is_none() && molecules.is_some());
        assert_eq!(r.partial_paths.len(), 1);
        assert!(partial
            .unwrap()
            .starts_with(">partial.1 circular=false length=120 path=U+ "));
    }

    #[test]
    fn an_unresolved_result_keeps_each_unplaced_repeat_once() {
        // X A is a linear piece; R is a repeat no bridge crosses. Its sequence is still in the
        // genome: with linear pieces it comes out once, as a piece of its own.
        let gfa = gfa_text(&[("X", 1500), ("A", 1400), ("R", 1200), ("S", 300)], &[("X+", "A+")]);
        let ev = json!({
            "repeat_segments": ["R", "S"],
            "segment_depths": {"X": 20.0, "A": 21.0, "R": 45.0, "S": 44.0},
            "median_depth": 20.0,
            "phased_paths": [],
            "links": link_reads(&[("X+", "A+", 30)]),
        });
        let (r, molecules, _) = run_files(&gfa, &ev, "ovasm-unplaced-repeat", 0);
        let paths: Vec<String> = r.components.iter()
            .flat_map(|c| c.best.iter().flat_map(|b| b.molecules.iter().map(|m| m.path.clone())))
            .collect();
        assert!(paths.contains(&"X+ A+".to_string()) && paths.contains(&"R+".to_string()), "{paths:?}");
        // S is a stub under MIN_UNBRIDGED_BP
        assert!(!paths.iter().any(|p| p.contains("S")));
        assert!(r.components.iter().any(|c| c.ranking_model == "unplaced_repeat"));
        assert_eq!(r.repeats_not_placed, names(&["S"]));
        assert!(molecules.unwrap().contains("path=R+"));

        // a circle: the unplaced repeat is a configuration the reads do not take
        let gfa = gfa_text(&[("X", 1500), ("A", 1400), ("R", 1200)], &[("X+", "A+"), ("A+", "X+")]);
        let ev = json!({
            "repeat_segments": ["R"],
            "segment_depths": {"X": 20.0, "A": 21.0, "R": 45.0},
            "median_depth": 20.0,
            "phased_paths": [],
            "links": link_reads(&[("X+", "A+", 30), ("A+", "X+", 30)]),
        });
        let (r, _, _) = run_files(&gfa, &ev, "ovasm-unplaced-repeat-circle", 0);
        assert!(r.best.as_ref().unwrap().circular);
        assert_eq!(r.molecules, 1);
        assert_eq!(r.repeats_not_placed, names(&["R"]));
    }


    fn repeat_closure_fixture() -> Graph {
        Graph {
            names: ["A", "R", "B", "S"].map(String::from).to_vec(),
            seqs: vec![b"ACGT".repeat(6); 4],
            links: vec![
                ((1, true), (1, true)),
                ((1, false), (1, false)),
                ((0, true), (1, true)),
                ((1, true), (2, true)),
                ((3, true), (3, false)),
            ],
            overlaps: vec![0; 5],
        }
    }

    fn repeat_closure_evidence(g: &Graph, direct: u64) -> Value {
        json!({
            "repeat_segments": ["R", "S"],
            "links": link_reads(&[("R-", "R-", direct), ("R+", "B+", 1), ("A+", "R+", 2)]),
            "observed_walks": [{"path":"R+ R+ B+", "reads":1}, {"path":"B- R- R-", "reads":1}],
            "phased_paths": [{"path":"A+ R+ B+", "reads":2}],
            "thread_evidence": {
                "schema": crate::thread_evidence::SCHEMA,
                "graph_fingerprint": crate::unzip::graph_fingerprint(g),
                "accepted_records": 1,
                "links": [
                    {"from":"R-", "to":"R-", "records":direct},
                    {"from":"R+", "to":"B+", "records":1},
                    {"from":"A+", "to":"R+", "records":2}
                ],
                "paths": [{"path":"R+ R+ B+", "reads":1}]
            }
        })
    }

    #[test]
    fn local_repeat_closures_keep_competing_links_and_distinct_contexts() {
        let g = repeat_closure_fixture();
        let ev = repeat_closure_evidence(&g, 3);
        let repeat = [false, true, false, true];
        assert!(crate::read_model::ReadModel::load(&g, &ev, &repeat, 2).unwrap().is_some());
        let report = local_repeat_closures(&g, &ev, &repeat, true, &Default::default(), 2, 100).unwrap();
        let [closure] = report.closures.as_slice() else { panic!("expected one same-strand loop") };
        assert_eq!((closure.segment.as_str(), closure.state, closure.closed_walk_length),
                   ("R", "unresolved_competing_context", Some(24)));
        assert_eq!(closure.link.single_record_crossings, Some(3));
        assert_eq!(closure.competing_links.len(), 2);
        assert!(closure.competing_links.iter().any(|e| e.single_record_crossings == Some(1)));
        assert_eq!(report.contexts.len(), 3);
        assert!(!report.contexts_truncated);
        assert_eq!(report.contexts.iter().find(|r| r.source == "observed_walk").unwrap().reported_records, 1);
        let limited = local_repeat_closures(&g, &ev, &repeat, true, &Default::default(), 2, 1).unwrap();
        assert!(limited.contexts_truncated);
        assert_eq!(limited.closures[0].competing_links.len(), 2);
        assert!(local_repeat_closures(&g, &ev, &[false; 4], true, &Default::default(), 2, 100).is_none());
    }

    #[test]
    fn local_repeat_closures_withhold_unverified_or_insufficient_crossings_and_ambiguous_overlap() {
        let mut g = repeat_closure_fixture();
        let repeat = [false, true, false, true];
        let ev = repeat_closure_evidence(&g, 1);
        let get = |g: &Graph, bound| local_repeat_closures(g, &ev, &repeat, bound, &Default::default(), 2, 100).unwrap();
        assert_eq!(get(&g, true).closures[0].state, "unresolved_insufficient_crossings");
        assert_eq!(get(&g, false).closures[0].state, "unresolved_no_graph_bound_crossings");
        assert!(get(&g, false).closures[0].closed_walk_length.is_none());
        g.overlaps[1] = 4;
        let ambiguous = get(&g, true);
        assert_eq!(ambiguous.closures[0].state, "unresolved_overlap");
        assert_eq!(ambiguous.closures[0].link.overlaps_bp, BTreeSet::from([0, 4]));
    }

    #[test]
    fn local_repeat_closures_respect_conflicting_hard_phase() {
        let g = repeat_closure_fixture();
        let repeat = [false, true, false, true];
        let ev = repeat_closure_evidence(&g, 3);
        let path = std::env::temp_dir().join(format!("ovasm-local-closure-phase-{}.json", std::process::id()));
        fs::write(&path, json!({
            "graph_fingerprint": crate::unzip::graph_fingerprint(&g),
            "min_count": 2,
            "unique_records": 10,
            "repeat_nodes": [{
                "segment":"R", "classification":"resolved",
                "entries":[{"segment":"A+","adjacency_reads":3},{"segment":"R+","adjacency_reads":3}],
                "exits":[{"segment":"R+","adjacency_reads":3},{"segment":"B+","adjacency_reads":3}],
                "phases":[{"entry":"A+","exit":"R+","reads":3,"supported":true},
                          {"entry":"R+","exit":"B+","reads":3,"supported":true}]
            }]
        }).to_string()).unwrap();
        let phase = crate::phase::Constraints::load(&g, &path).unwrap();
        fs::remove_file(path).unwrap();
        assert!(!phase.allows_walk(&[(1, true)], true));
        let report = local_repeat_closures(&g, &ev, &repeat, true, &phase, 2, 100).unwrap();
        assert_eq!(report.closures[0].state, "unresolved_phase_conflict");
        assert_eq!(report.closures[0].closed_walk_length, None);
        assert_eq!(report.closures[0].competing_links.len(), 2);
        assert_eq!(report.contexts.iter().filter(|row| row.source == "phase_report").count(), 2);
    }

    #[test]
    fn scoring_budget_scales_only_for_bounded_components() {
        let base = crate::read_model::WORK_BUDGET;
        assert_eq!(scoring_work_budget(53, 0), base);
        assert_eq!(scoring_work_budget(53, 20_000), base);
        assert_eq!(scoring_work_budget(53, 44_432), base * 3);
        assert_eq!(scoring_work_budget(53, 100_000), base * 5);
        assert_eq!(scoring_work_budget(53, 100_001), base);
        assert_eq!(scoring_work_budget(65, 44_432), base);
        assert_eq!(candidate_work_budget(1_000_000, 53, 44_432), 3_000_000);
        assert_eq!(candidate_work_budget(1_000_000, 65, 44_432), 1_000_000);
    }

    #[test]
    fn bounded_retry_marks_candidate_limit_as_truncated() {
        let anchors = [0usize];
        let optional = [false];
        let br = vec![Bridge { walk: vec![(0, true), (0, true)], reads: 2, weak: false, topological: false }];
        let mut out = FxHashMap::default();
        out.insert((0, true), vec![0]);
        let phase = crate::phase::Constraints::default();
        let unbounded = search_phased_indexed_impl(1, &anchors, &optional, &br, 1, 1_000, &phase, &out, None, true);
        assert!(!unbounded.0.is_empty() && !unbounded.1);
        let capped = search_phased_indexed_bounded(1, &anchors, &optional, &br, 1, 1_000, &phase, &out, None, true, Some(0));
        assert!(capped.0.is_empty() && capped.1);
        let enough = search_phased_indexed_bounded(1, &anchors, &optional, &br, 1, 1_000, &phase, &out, None, true, Some(unbounded.0.len()));
        assert!(!enough.1 && !enough.2 && enough.0.len() == unbounded.0.len());
    }

    #[test]
    fn optimized_cover_search_matches_unpruned_oracle() {
        // Deterministic directed four-anchor graphs, including self closures,
        // optional anchors, both orientations, and active phase checks.
        let anchors = [0usize, 1, 2, 3];
        let mut seed = 0x1234_5678u64;
        for case in 0..64 {
            let mut br = Vec::new();
            for a in 0..4 { for b in 0..4 {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                if seed >> 62 == 0 {
                    for forward in [true, false] {
                        br.push(Bridge { walk: vec![(a, forward), (b, forward)], reads: (case % 7 + 1) as u64, weak: false, topological: false });
                    }
                }
            }}
            let mut out = FxHashMap::default();
            for (i, b) in br.iter().enumerate() { out.entry(b.walk[0]).or_insert_with(Vec::new).push(i); }
            let mut phase = crate::phase::Constraints::default();
            if case % 2 == 0 { phase.active_nodes = 1; }
            for mask in 0..16 {
                let optional: Vec<bool> = (0..4).map(|i| mask & (1 << i) != 0).collect();
                for m in 1..=4 {
                    let a = search_phased_indexed_impl(4, &anchors, &optional, &br, m, 2_000_000, &phase, &out, None, true);
                    let b = search_phased_indexed_impl(4, &anchors, &optional, &br, m, 2_000_000, &phase, &out, None, false);
                    let repr = |c: &[Vec<Found>]| c.iter().map(|p| p.iter().map(|x| (x.start, x.bridges.clone(), x.circular)).collect::<Vec<_>>()).collect::<Vec<_>>();
                    assert!(!a.1 && !b.1, "fixture budget case={case} mask={mask} m={m}");
                    assert_eq!(repr(&a.0), repr(&b.0), "case={case} mask={mask} m={m}");
                    assert_eq!(a.2, b.2);
                }
            }
        }
    }

    #[test]
    fn optimized_cover_search_matches_phase_and_terminal_oracle() {
        let gfa = gfa_text(&[("A",100),("C",100),("R",100),("B",100),("D",100)],
            &[("A+","R+"),("C+","R+"),("R+","B+"),("R+","D+")]);
        let dir = std::env::temp_dir().join(format!("ovasm-dfs-phase-oracle-{}",std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let gp = dir.join("graph.gfa"); std::fs::write(&gp,&gfa).unwrap();
        let graph = parse_gfa(&gp).unwrap();
        let report = json!({"graph_fingerprint":crate::unzip::graph_fingerprint(&graph), "min_count":80, "unique_records":400,
            "repeat_nodes":[{"segment":"R","classification":"resolved",
            "entries":[{"segment":"A+","adjacency_reads":200},{"segment":"C+","adjacency_reads":200}],
            "exits":[{"segment":"B+","adjacency_reads":200},{"segment":"D+","adjacency_reads":200}],
            "phases":[{"entry":"A+","exit":"B+","reads":80,"supported":true},{"entry":"C+","exit":"D+","reads":80,"supported":true}]}]});
        let pp=dir.join("phase.json"); std::fs::write(&pp,report.to_string()).unwrap();
        let phases=crate::phase::Constraints::load(&graph,&pp).unwrap();
        assert!(!phases.allows((0,true),(2,true),(4,true)));
        assert!(!phases.allows((4,false),(2,false),(0,false)));
        let links=[("A+","R+",2),("C+","R+",2),("R+","B+",2),("R+","D+",2)];
        let ev=terminal_evidence(&gfa,&links,&[("A+ R+ B+",2),("C+ R+ D+",2)]);
        let anchors: Vec<_>=(0..5).collect();let optional=vec![false,false,false,true,true];let repeats=vec![false;5];
        let threads=crate::read_model::ReadModel::load(&graph,&ev,&repeats,2).unwrap().unwrap().thread_ends(&anchors,&optional);
        assert!(!threads.is_empty());
        let (br,_)=bridges(&graph,&ev,&repeats,2).unwrap();let mut out=FxHashMap::default();
        for (i,b) in br.iter().enumerate(){out.entry(b.walk[0]).or_insert_with(Vec::new).push(i);}
        for m in 1..=3 { for with_threads in [false,true] {
            let t=with_threads.then_some(&threads);
            let a=search_phased_indexed_impl(5,&anchors,&optional,&br,m,2_000_000,&phases,&out,t,true);
            let b=search_phased_indexed_impl(5,&anchors,&optional,&br,m,2_000_000,&phases,&out,t,false);
            let repr=|c:&[Vec<Found>]| c.iter().map(|p|p.iter().map(|f|(f.start,f.bridges.clone(),f.circular)).collect::<Vec<_>>()).collect::<Vec<_>>();
            assert!(!a.1 && !b.1 && !a.2 && !b.2);
            assert_eq!(repr(&a.0),repr(&b.0),"m={m} threads={with_threads}");
        }}
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn partial_paths_take_the_best_bridge_at_each_end_and_close_what_closes() {
        // A -> B -> C; A also leaves for D on fewer reads: A's end goes to B, D stays alone.
        let g = anchors_graph(
            &["A", "B", "C", "D"],
            &[
                ((0, true), (1, true)),
                ((1, true), (2, true)),
                ((2, true), (0, true)),
                ((0, true), (3, true)),
            ],
        );
        let ev = json!({
            "repeat_segments": [],
            "phased_paths": [],
            "links": link_reads(&[("A+", "B+", 20), ("B+", "C+", 18), ("A+", "D+", 4)]),
        });
        let (br, _) = bridges(&g, &ev, &[false; 4], 2).unwrap();
        let frac = exit_fractions(&br);
        let paths = partial_paths(4, &[0, 1, 2, 3], &[false; 4], &br, &frac);
        let walks: Vec<Vec<Oriented>> = paths.iter().map(|f| candidate_walk(&br, f)).collect();
        assert_eq!(
            walks,
            vec![vec![(0, true), (1, true), (2, true)], vec![(3, true)]]
        );
        assert!(paths.iter().all(|f| !f.circular));

        // A -> B -> C and back to A, all read: the path closes
        let ev = json!({
            "repeat_segments": [],
            "phased_paths": [],
            "links": link_reads(&[("A+", "B+", 20), ("B+", "C+", 18), ("C+", "A+", 15)]),
        });
        let (br, _) = bridges(&g, &ev, &[false; 4], 2).unwrap();
        let frac = exit_fractions(&br);
        let paths = partial_paths(4, &[0, 1, 2], &[false; 4], &br, &frac);
        assert_eq!(paths.len(), 1);
        assert!(paths[0].circular && paths[0].bridges.len() == 3);
        assert_eq!(
            candidate_walk(&br, &paths[0]).len(),
            3,
            "the first anchor once"
        );

        // an optional anchor inside a path stays, at its end it goes; a path of optional
        // anchors alone is not reported
        let mut optional = [false; 4];
        optional[1] = true;
        optional[3] = true;
        let ev = json!({
            "repeat_segments": [],
            "phased_paths": [],
            "links": link_reads(&[("A+", "B+", 20), ("B+", "C+", 18)]),
        });
        let (br, _) = bridges(&g, &ev, &[false; 4], 2).unwrap();
        let frac = exit_fractions(&br);
        let paths = partial_paths(4, &[0, 1, 2], &optional, &br, &frac);
        assert_eq!(paths.len(), 1);
        assert_eq!(candidate_walk(&br, &paths[0]).len(), 3);
        optional[2] = true;
        let paths = partial_paths(4, &[0, 1, 2], &optional, &br, &frac);
        assert_eq!(paths.len(), 1);
        assert_eq!(candidate_walk(&br, &paths[0]), vec![(0, true)]);
    }

    #[test]
    fn a_solved_graph_leaves_no_partial_file() {
        // the plastid of `run_plastid`: one circle, everything placed
        let gfa = gfa_text(
            &[("u0", 400), ("u1", 100), ("u2", 60)],
            &[
                ("u0-", "u2+"),
                ("u0+", "u2+"),
                ("u1-", "u2-"),
                ("u1+", "u2-"),
            ],
        );
        let ev = json!({
            "repeat_segments": ["u2"],
            "segment_depths": {"u0": 32.8, "u1": 39.1, "u2": 104.5},
            "median_depth": 32.8,
            "phased_paths": [
                {"path": "u0+ u2+ u1-", "reads": 12},
                {"path": "u0- u2+ u1+", "reads": 13},
            ],
            "links": link_reads(&[
                ("u0-", "u2+", 82),
                ("u0+", "u2+", 84),
                ("u1-", "u2-", 90),
                ("u1+", "u2-", 80),
            ]),
        });
        let (r, molecules, partial) = run_files(&gfa, &ev, "ovasm-no-partial", 4);
        assert!(r.best.as_ref().unwrap().circular);
        assert!(molecules
            .unwrap()
            .starts_with(">linear circular=true length=620 "));
        assert!(partial.is_none(), "the earlier run's file is removed");
        assert!(r.partial_paths.is_empty() && r.partial_fasta.is_none());
        assert!(r.unplaced_segments.is_empty() && r.unsolved_components == 0);
    }

    #[test]
    fn a_once_read_pairing_closes_an_end_nothing_else_leaves() {
        let r = run_plastid(
            json!([
                {"path": "u0+ u2+ u1-", "reads": 2},
                {"path": "u0- u2+ u1+", "reads": 1},
                {"path": "u0- u2+ u1-", "reads": 1},
            ]),
            "ovasm-weak",
        );
        let best = r.best.as_ref().unwrap();
        assert!(best.circular, "{}", best.path);
        assert_eq!(best.length, 400 + 100 + 2 * 60, "both IR copies");
        assert!(r.weak_bridges >= 1);
        assert!(
            !r.decisive,
            "a once-read pairing does not decide the configuration"
        );
        assert_eq!(r.components[0].best.as_ref().unwrap().weak_bridges, 1);
    }

    #[test]
    fn a_dead_end_closes_even_where_the_far_end_is_bridged() {
        // Holcus plastid: both flip-flop isomers (4 and 3 reads) leave the same LSC end, so the
        // other LSC end has only a once-read pairing, whose SSC end the isomers already leave.
        let r = run_plastid(
            json!([
                {"path": "u0+ u2+ u1-", "reads": 4},
                {"path": "u0+ u2+ u1+", "reads": 3},
                {"path": "u0- u2+ u1+", "reads": 1},
            ]),
            "ovasm-holcus",
        );
        let best = r.best.as_ref().unwrap();
        assert!(best.circular, "{}", best.path);
        assert_eq!(best.length, 400 + 100 + 2 * 60, "both IR copies");
        assert!(!r.decisive);
    }

    /// Linearize a GFA and evidence report given as text and JSON.
    /// A and C enter repeat R, B and D leave it, B leads on to C and D back to A: one circle
    /// A R B C R D, or two, A R D and C R B. `phased`: reads across R.
    fn prior_run(phased: &[(&str, u64)], prior: &[(&str, f64)], dir_name: &str) -> LinearizeReport {
        prior_run_with(phased, &[], prior, dir_name)
    }

    /// `threads`: single records of the target (one read through three segments or more), as
    /// short reads give them: no phased path, so the bridges through R stay the graph's pairings.
    fn prior_run_with(phased: &[(&str, u64)], threads: &[(&str, u64)], prior: &[(&str, f64)], dir_name: &str) -> LinearizeReport {
        let s = |n, seed| String::from_utf8(seq(n, seed)).unwrap();
        let gfa = format!(
            "H\tVN:Z:1.0\nS\tA\t{}\nS\tB\t{}\nS\tC\t{}\nS\tD\t{}\nS\tR\t{}\n\
             L\tA\t+\tR\t+\t0M\nL\tC\t+\tR\t+\t0M\nL\tR\t+\tB\t+\t0M\nL\tR\t+\tD\t+\t0M\n\
             L\tB\t+\tC\t+\t0M\nL\tD\t+\tA\t+\t0M\n",
            s(400, 1),
            s(300, 2),
            s(350, 3),
            s(250, 4),
            s(120, 5)
        );
        let mut ev = json!({
            "repeat_segments": ["R"],
            "segment_depths": {"A": 50.0, "B": 50.0, "C": 50.0, "D": 50.0, "R": 100.0},
            "median_depth": 50.0,
            "phased_paths": phased.iter().map(|(p, n)| json!({"path": p, "reads": n})).collect::<Vec<_>>(),
            "links": [{"from": "B+", "to": "C+", "reads": 20}, {"from": "D+", "to": "A+", "reads": 20}],
        });
        let dir = std::env::temp_dir().join(format!("{dir_name}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("g.gfa"), &gfa).unwrap();
        let g = parse_gfa(&dir.join("g.gfa")).unwrap();
        if !threads.is_empty() {
            let links = [("A+", "R+"), ("C+", "R+"), ("R+", "B+"), ("R+", "D+"), ("B+", "C+"), ("D+", "A+")];
            ev["thread_evidence"] = json!({"schema": crate::thread_evidence::SCHEMA,
                "graph_fingerprint": crate::unzip::graph_fingerprint(&g),
                "accepted_records": threads.iter().map(|(_, n)| n).sum::<u64>(),
                "links": links.iter().map(|(a, b)| json!({"from": a, "to": b, "records": 20})).collect::<Vec<_>>(),
                "paths": threads.iter().map(|(p, n)| json!({"path": p, "reads": n})).collect::<Vec<_>>()});
        }
        std::fs::write(dir.join("e.json"), ev.to_string()).unwrap();
        let index: FxHashMap<&str, usize> = g.names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
        let prior = BridgePrior {
            probability: prior
                .iter()
                .map(|(labels, p)| (labels.split(' ').map(|l| parse_label(&index, l).unwrap()).collect(), *p))
                .collect(),
        };
        let params = LinearizeParams { min_reads: 2, max_steps: 1_000_000, adaptive_steps: false, max_molecules: 4, max_alternatives: 10, optional_below: 0.5 };
        let r = run_with_prior(&dir.join("g.gfa"), &dir.join("e.json"), params, &dir.join("m.fa"), None, &dir.join("l.json"), &prior).unwrap();
        // the bridges a caller weighs are the solver's own
        let (_, repeat, views) = open_bridges(&dir.join("g.gfa"), &dir.join("e.json"), 2).unwrap();
        assert_eq!(repeat, vec![false, false, false, false, true]);
        assert!(views.iter().any(|b| b.walk.len() == 2 && !b.topological && !b.weak && b.reads == 20));
        assert!(views.iter().all(|b| !repeat[b.walk[0].0] && !repeat[b.walk[b.walk.len() - 1].0]));
        std::fs::remove_dir_all(&dir).ok();
        r
    }

    #[test]
    fn a_prior_takes_one_of_the_covers_reads_leave_open() {
        // no read crosses R: the search's representative is the cover with the fewest molecules
        let two_circles = [("A+ R+ B+", 0.1), ("A+ R+ D+", 0.9), ("C+ R+ B+", 0.9), ("C+ R+ D+", 0.1)];
        let r = prior_run(&[], &two_circles, "ovasm-prior-open");
        let c = &r.components[0];
        assert_eq!(c.best.as_ref().unwrap().molecules.len(), 1);
        assert!(!c.decisive);
        let p = c.prior_cover.as_ref().unwrap();
        assert_eq!(p.cover.molecules.len(), 2);
        assert!(p.cover.molecules.iter().all(|m| m.circular && m.length == 770), "{:?}", p.cover.molecules.iter().map(|m| &m.path).collect::<Vec<_>>());
        assert_eq!((p.ends, p.open_ends, p.decided_ends, p.equally_good), (8, 4, 2, 1));
        // A+ and C+ at 0.9, their bridges' other ends unknown to the prior (1/2), the two read links 1 at each end
        assert!((p.expected_correct_ends - 6.8).abs() < 1e-9, "{}", p.expected_correct_ends);
        assert_eq!(p.cover.topological_bridges, 2);
        // a prior that tells nothing apart: the fewest molecules of the tied covers, as the search
        let r = prior_run(&[], &[], "ovasm-prior-flat");
        let p = r.components[0].prior_cover.as_ref().unwrap();
        assert_eq!((p.cover.molecules.len(), p.decided_ends, p.equally_good), (1, 0, 2));
        assert_eq!(p.cover.molecules[0].length, 400 + 120 + 300 + 350 + 120 + 250);
    }

    #[test]
    fn a_prior_never_outranks_what_reads_measured() {
        let two_circles = [("A+ R+ B+", 0.1), ("A+ R+ D+", 0.9), ("C+ R+ B+", 0.9), ("C+ R+ D+", 0.1)];
        // few reads across R: the unread pairings stay possible, but a cover of read bridges comes first
        let r = prior_run(&[("A+ R+ B+", 3), ("C+ R+ D+", 4)], &two_circles, "ovasm-prior-few");
        let p = r.components[0].prior_cover.as_ref().unwrap();
        assert_eq!((p.cover.molecules.len(), p.cover.topological_bridges), (1, 0));
        // many reads: the unread pairings are no bridges at all
        let r = prior_run(&[("A+ R+ B+", 30), ("C+ R+ D+", 20)], &two_circles, "ovasm-prior-many");
        let p = r.components[0].prior_cover.as_ref().unwrap();
        assert_eq!((p.cover.molecules.len(), p.open_ends, p.equally_good), (1, 0, 1));
        assert_eq!(p.cover.molecules[0].path, r.components[0].best.as_ref().unwrap().molecules[0].path);
    }

    #[test]
    fn the_targets_single_records_come_before_a_prior() {
        let two_circles = [("A+ R+ B+", 0.1), ("A+ R+ D+", 0.9), ("C+ R+ B+", 0.9), ("C+ R+ D+", 0.1)];
        // five reads of the target run A-R-B in one piece: the pairing the prior likes less
        let r = prior_run_with(&[], &[("A+ R+ B+", 5)], &two_circles, "ovasm-prior-thread");
        let p = r.components[0].prior_cover.as_ref().unwrap();
        assert_eq!((p.cover.molecules.len(), p.thread_records_explained, p.joint_records_explained), (1, 5, 0));
        assert_eq!(p.cover.molecules[0].length, 400 + 120 + 300 + 350 + 120 + 250);
        // the prior still tells the bridges apart, and the cover's fit to the records is reported
        assert_eq!((p.open_ends, p.decided_ends), (4, 2));
        assert_eq!(p.cover.read_support.as_ref().map(|f| f.thread_records_explained), Some(5));
        // one read is no record to go by (`min_reads` 2): the prior has its way
        let r = prior_run_with(&[], &[("A+ R+ B+", 1)], &two_circles, "ovasm-prior-thread-one");
        let p = r.components[0].prior_cover.as_ref().unwrap();
        assert_eq!((p.cover.molecules.len(), p.thread_records_explained), (2, 0));
    }

    #[test]
    fn a_prior_cover_leaves_a_shallow_anchor_out_and_ends_where_the_graph_stops() {
        let s = |n, seed| String::from_utf8(seq(n, seed)).unwrap();
        let params = LinearizeParams { min_reads: 2, max_steps: 1_000_000, adaptive_steps: false, max_molecules: 4, max_alternatives: 10, optional_below: 0.5 };
        let solve = |gfa: &str, ev: Value, name: &str| {
            let dir = std::env::temp_dir().join(format!("{name}-{}", std::process::id()));
            std::fs::create_dir_all(&dir).unwrap();
            std::fs::write(dir.join("g.gfa"), gfa).unwrap();
            std::fs::write(dir.join("e.json"), ev.to_string()).unwrap();
            let prior = BridgePrior { probability: FxHashMap::default() };
            let r = run_with_prior(&dir.join("g.gfa"), &dir.join("e.json"), params, &dir.join("m.fa"), None, &dir.join("l.json"), &prior).unwrap();
            std::fs::remove_dir_all(&dir).ok();
            r
        };
        // A and B in a circle, with a shallow segment o between them on a second way
        let gfa = format!(
            "H\tVN:Z:1.0\nS\tA\t{}\nS\tB\t{}\nS\to\t{}\nL\tA\t+\to\t+\t0M\nL\to\t+\tB\t+\t0M\nL\tA\t+\tB\t+\t0M\nL\tB\t+\tA\t+\t0M\n",
            s(400, 1),
            s(300, 2),
            s(200, 3)
        );
        let ev = json!({"repeat_segments": [], "segment_depths": {"A": 50.0, "B": 50.0, "o": 10.0}, "median_depth": 50.0, "phased_paths": [],
            "links": [{"from": "A+", "to": "o+", "reads": 5}, {"from": "o+", "to": "B+", "reads": 5},
                {"from": "A+", "to": "B+", "reads": 40}, {"from": "B+", "to": "A+", "reads": 45}]});
        let r = solve(&gfa, ev, "ovasm-prior-optional");
        assert_eq!(r.optional_anchors, vec!["o".to_string()]);
        let p = r.components[0].prior_cover.as_ref().unwrap();
        assert_eq!(p.cover.molecules.len(), 1);
        assert_eq!((p.cover.molecules[0].length, p.cover.molecules[0].circular, p.ends), (700, true, 4));
        // A, B and C in a row: the graph stops at both ends, and the cover is that one linear molecule
        let gfa = format!(
            "H\tVN:Z:1.0\nS\tA\t{}\nS\tB\t{}\nS\tC\t{}\nL\tA\t+\tB\t+\t0M\nL\tB\t+\tC\t+\t0M\n",
            s(400, 1),
            s(300, 2),
            s(200, 3)
        );
        let ev = json!({"repeat_segments": [], "segment_depths": {"A": 50.0, "B": 50.0, "C": 50.0}, "median_depth": 50.0, "phased_paths": [],
            "links": [{"from": "A+", "to": "B+", "reads": 40}, {"from": "B+", "to": "C+", "reads": 45}]});
        let r = solve(&gfa, ev, "ovasm-prior-linear");
        let p = r.components[0].prior_cover.as_ref().unwrap();
        assert_eq!(p.cover.molecules.len(), 1);
        assert_eq!((p.cover.molecules[0].length, p.cover.molecules[0].circular, p.ends), (900, false, 4));
        assert!(["A+ B+ C+", "C- B- A-"].contains(&p.cover.molecules[0].path.as_str()), "{}", p.cover.molecules[0].path);
    }

    fn run_graph(gfa: &str, ev: Value, dir_name: &str) -> LinearizeReport {
        let dir = std::env::temp_dir().join(format!("{dir_name}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("g.gfa"), gfa).unwrap();
        std::fs::write(dir.join("e.json"), ev.to_string()).unwrap();
        let params = LinearizeParams {
            min_reads: 2,
            max_steps: 1_000_000, adaptive_steps: false,
            max_molecules: 4,
            max_alternatives: 10,
            optional_below: 0.5,
        };
        let r = run(
            &dir.join("g.gfa"),
            &dir.join("e.json"),
            params,
            &dir.join("m.fa"),
            None,
            &dir.join("l.json"),
        )
        .unwrap();
        std::fs::remove_dir_all(&dir).ok();
        r
    }

    #[test]
    fn reads_across_a_tandem_array_say_how_often_its_loop_is_laid_down() {
        // A circle of one 1,000 bp anchor and a 300 bp repeat that loops on itself. Depth says
        // two copies (100x against 50x); twelve reads hold 2,990-3,012 bp between the anchor's
        // two ends, ten copies, more than the topological bridges go round (`TOPO_MAX_REPEATS`).
        let gfa = gfa_text(&[("A", 1000), ("R", 300)], &[("A+", "R+"), ("R+", "R+"), ("R+", "A+")]);
        let ev = |spans: Option<Vec<i64>>| {
            let mut v = json!({
                "repeat_segments": ["R"],
                "segment_depths": {"A": 50.0, "R": 100.0},
                "median_depth": 50.0,
                "phased_paths": [],
                "links": link_reads(&[("A+", "R+", 20), ("R+", "R+", 20), ("R+", "A+", 20)]),
                "params": {"spacing_tol": 0.05},
            });
            if let Some(s) = spans {
                v["anchor_spans"] = json!([{"from": "A+", "to": "A+", "spans": s}]);
            }
            v
        };
        let length = |r: &LinearizeReport| r.components[0].best.as_ref().unwrap().molecules[0].length;
        let by_depth = run_graph(&gfa, ev(None), "ovasm-array-depth");
        assert!(length(&by_depth) <= 1000 + TOPO_MAX_REPEATS * 300);
        let spans: Vec<i64> = (0..12).map(|i| 2990 + 2 * i).collect();
        let by_reads = run_graph(&gfa, ev(Some(spans)), "ovasm-array-reads");
        let best = by_reads.components[0].best.as_ref().unwrap();
        assert!(best.molecules[0].circular);
        assert_eq!(best.molecules[0].length, 1000 + 10 * 300);
        assert_eq!(best.topological_bridges, 0, "the spanned walk is read-measured");
        // too few reads (min_reads 2) leave the graph's answer
        let one = run_graph(&gfa, ev(Some(vec![3000])), "ovasm-array-one");
        assert_eq!(length(&one), length(&by_depth));
    }

    #[test]
    fn an_anchor_between_both_ends_of_an_unspanned_repeat_stays_in_the_cover() {
        // Buddleja mitochondrion: u1 sits between both ends of repeat u4, which no read spans
        // into u1; reads cross u4 only from u3 to u2. All anchor sequence and both repeat
        // copies must be reported, but these reads do not establish one physical circle.
        let s = |n, seed| String::from_utf8(seq(n, seed)).unwrap();
        let gfa = format!(
            "H\tVN:Z:1.0\nS\tu0\t{}\nS\tu1\t{}\nS\tu2\t{}\nS\tu3\t{}\nS\tu4\t{}\n\
             L\tu0\t-\tu2\t-\t0M\nL\tu0\t+\tu3\t+\t0M\nL\tu1\t-\tu4\t-\t0M\n\
             L\tu1\t+\tu4\t+\t0M\nL\tu2\t-\tu4\t+\t0M\nL\tu3\t+\tu4\t-\t0M\n",
            s(400, 1),
            s(300, 2),
            s(100, 3),
            s(200, 4),
            s(60, 5)
        );
        let ev = json!({
            "repeat_segments": ["u4"],
            "segment_depths": {"u0": 36.0, "u1": 34.0, "u2": 38.0, "u3": 39.0, "u4": 55.0},
            "median_depth": 36.0,
            "phased_paths": [{"path": "u3+ u4- u2+", "reads": 6}],
            "links": [
                {"from": "u0+", "to": "u3+", "reads": 20},
                {"from": "u0-", "to": "u2-", "reads": 20},
            ],
        });
        let r = run_graph(&gfa, ev, "ovasm-buddleja");
        let cover = r.components[0].best.as_ref().unwrap();
        assert!(cover.molecules.iter().all(|m| m.circular));
        assert_eq!(
            cover.molecules.iter().map(|m| m.length).sum::<usize>(),
            400 + 200 + 60 + 300 + 60 + 100,
            "all anchors and both repeat copies must remain in the cover"
        );
        assert!(cover.molecules.iter().any(|m| m.path.contains("u1")));
        assert!(r.unbridged_anchors.is_empty());
        assert!(!r.decisive);
    }

    #[test]
    fn well_read_pairings_need_no_fallback() {
        // "well read" means at least `ABSENCE_READS` reads crossing from each end; below it a
        // pairing nobody read may just not have been sampled (see
        // `few_reads_across_a_repeat_do_not_rule_out_the_unread_pairings`)
        let r = run_plastid(
            json!([
                {"path": "u0+ u2+ u1-", "reads": 12},
                {"path": "u0- u2+ u1+", "reads": 15},
                {"path": "u0- u2+ u1-", "reads": 6},
            ]),
            "ovasm-strong",
        );
        let best = r.best.as_ref().unwrap();
        assert!(best.circular);
        assert_eq!(best.length, 400 + 100 + 2 * 60);
        assert_eq!(r.weak_bridges, 0);
        assert!(r.decisive);
    }

    #[test]
    fn an_unbridged_contig_is_a_linear_molecule_of_its_own() {
        // Thuidium mitochondrion at ~14x: one unitig, no links; a shallow piece beside it.
        let s = |n, seed| String::from_utf8(seq(n, seed)).unwrap();
        let gfa = format!("H\tVN:Z:1.0\nS\tu0\t{}\nS\tu1\t{}\n", s(1500, 1), s(80, 2));
        let ev = json!({
            "repeat_segments": [],
            "segment_depths": {"u0": 14.0, "u1": 2.0},
            "median_depth": 14.0,
            "phased_paths": [],
            "links": [],
        });
        let r = run_graph(&gfa, ev, "ovasm-unbridged");
        assert_eq!(r.molecules, 1);
        let best = r.best.as_ref().unwrap();
        assert_eq!(
            (best.path.as_str(), best.circular, best.length),
            ("u0+", false, 1500)
        );
        assert!(!r.decisive);
        assert_eq!(r.unbridged_anchors, vec!["u0", "u1"]);
        assert_eq!(r.skipped_anchors, vec!["u1"], "the shallow piece stays out");
    }

    #[test]
    fn labels_parse_and_reject_unknown_segments() {
        let g = graph();
        let index: FxHashMap<&str, usize> = g
            .names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.as_str(), i))
            .collect();
        assert_eq!(parse_label(&index, "R-").unwrap(), (4, false));
        assert!(parse_label(&index, "Q+").is_err());
        assert!(parse_label(&index, "A").is_err());
    }

    #[test]
    fn an_end_without_a_full_bridge_takes_its_weak_walks_before_the_graphs_pairings() {
        // Reads cross R from A (30), but C's only pairing holds 1 read, below min_reads 2:
        // R counts as crossed and C+ has no full bridge (Geg-14's 18.9 kb repeat). The walk a
        // read did take from C+ is its bridge, weak; no topological guess is added.
        let g = graph();
        let repeat = vec![false, false, false, false, true];
        let ev = evidence(&[("A+ R+ B+", 30), ("C+ R+ D+", 1)], &[]);
        let (br, _) = bridges(&g, &ev, &repeat, 2).unwrap();
        let from_c: Vec<&Bridge> = br.iter().filter(|b| b.walk[0] == (2, true)).collect();
        // the walk a read did take stays a weak read bridge, next to the graph's other walks
        let taken = from_c
            .iter()
            .find(|b| b.walk[b.walk.len() - 1] == (3, true))
            .unwrap();
        assert!(taken.weak && !taken.topological);

        // With no read at all from C+, the graph's pairings stand in.
        let ev = evidence(&[("A+ R+ B+", 30)], &[]);
        let (br, _) = bridges(&g, &ev, &repeat, 2).unwrap();
        let from_c: Vec<&Bridge> = br.iter().filter(|b| b.walk[0] == (2, true)).collect();
        assert!(!from_c.is_empty() && from_c.iter().all(|b| b.topological));
        let mut exits: Vec<Oriented> = from_c.iter().map(|b| b.walk[b.walk.len() - 1]).collect();
        exits.sort();
        assert_eq!(
            exits,
            vec![(1, true), (3, true)],
            "every exit the graph allows"
        );
        // A+ keeps its read bridge; the reverse of D's topological pairings reaches it, but
        // the shares reads measured are unchanged
        let read_a: Vec<&Bridge> = br
            .iter()
            .filter(|b| b.walk[0] == (0, true) && !b.topological)
            .collect();
        assert_eq!(read_a.len(), 1);
        assert_eq!(read_a[0].reads, 30);
        let frac = exit_fractions(&br);
        let at = |w: &[Oriented]| frac[br.iter().position(|b| b.walk == w).unwrap()];
        assert!((at(&[(0, true), (4, true), (1, true)]) - 1.0).abs() < 1e-12);
        assert!(at(&[(0, true), (4, true), (3, true)]) < 0.2);
    }

    #[test]
    fn few_reads_across_a_repeat_do_not_rule_out_the_unread_pairings() {
        let g = graph();
        let repeat = vec![false, false, false, false, true];
        let pairs = |br: &[Bridge]| -> Vec<(Oriented, Oriented, bool)> {
            let mut v: Vec<_> = br
                .iter()
                .filter(|b| b.walk.len() == 3 && b.walk[0].1 && b.walk[2].1)
                .map(|b| (b.walk[0], b.walk[2], b.topological))
                .collect();
            v.sort();
            v
        };
        // 3 and 4 reads: A->D and C->B, never read, may just not have been sampled
        let ev = evidence(&[("A+ R+ B+", 3), ("C+ R+ D+", 4)], &[]);
        let (br, _) = bridges(&g, &ev, &repeat, 2).unwrap();
        assert_eq!(
            pairs(&br),
            vec![
                ((0, true), (1, true), false),
                ((0, true), (3, true), true),
                ((2, true), (1, true), true),
                ((2, true), (3, true), false),
            ]
        );
        // 30 and 20 reads: the unread pairings are taken as absent
        let ev = evidence(&[("A+ R+ B+", 30), ("C+ R+ D+", 20)], &[]);
        let (br, _) = bridges(&g, &ev, &repeat, 2).unwrap();
        assert!(pairs(&br).iter().all(|p| !p.2));
    }

    #[test]
    fn unconstrained_fast_path_preserves_search_candidates_and_budget() {
        let g = graph();
        let ev = evidence(&[("A+ R+ B+", 30), ("C+ R+ D+", 20)], &[]);
        let (br, _) = bridges(&g, &ev, &[false, false, false, false, true], 2).unwrap();
        let fast = crate::phase::Constraints::default();
        let mut full_checks = crate::phase::Constraints::default();
        full_checks.active_nodes = 1; // force walk construction, with no forbidden triples
        let keys = |covers: Vec<Vec<Found>>| {
            covers.into_iter().map(|cover| cover.into_iter().map(|f|
                (f.circular, candidate_walk(&br, &f))).collect::<Vec<_>>()).collect::<Vec<_>>()
        };
        for budget in [0, 10, 10000] {
            for m in [1, 2, 3] {
                let (a, at) = search_phased(5, &[0, 1, 2, 3], &[false; 5], &br, m, budget, &fast);
                let (b, bt) = search_phased(5, &[0, 1, 2, 3], &[false; 5], &br, m, budget, &full_checks);
                assert_eq!(at, bt);
                assert_eq!(keys(a), keys(b));
            }
        }
    }

    #[test]
    fn indexed_lengths_match_spelling_in_both_orientations() {
        let g = Graph {
            names: vec!["A".into(), "B".into()],
            seqs: vec![b"ACGTAC".to_vec(), b"ACGGTA".to_vec()],
            links: vec![((0, true), (1, true)), ((1, true), (0, true))],
            overlaps: vec![2, 1],
        };
        let index = overlap_index(&g);
        for walk in [vec![(0, true), (1, true)], vec![(1, false), (0, false)]] {
            for circular in [false, true] {
                assert_eq!(spelled_length(&g, &walk, circular, &index).unwrap(),
                           spell(&g, &walk, circular).unwrap().len());
            }
        }
        assert!(spelled_length(&g, &[(0, true), (0, true)], false, &index).is_err());
        let mut duplicates = g;
        duplicates.links.extend([((0, true), (1, true)), ((1, false), (0, false)), ((0, true), (0, true))]);
        duplicates.overlaps.extend([4, 3, 2]);
        let index = overlap_index(&duplicates);
        for walk in [vec![(0, true), (1, true)], vec![(1, false), (0, false)], vec![(0, true)]] {
            for circular in [false, true] {
                assert_eq!(spelled_length(&duplicates, &walk, circular, &index).unwrap(),
                           spell(&duplicates, &walk, circular).unwrap().len());
            }
        }
    }

    #[test]
    fn exhausted_search_does_not_expand_to_larger_covers() {
        let dir = std::env::temp_dir().join(format!("ovasm-budget-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("g.gfa"),
            gfa_text(&[("A", 600), ("B", 500), ("C", 400)], &[("A+", "B+"), ("B+", "C+"), ("C+", "A+")])).unwrap();
        let ev = json!({"repeat_segments": [], "phased_paths": [], "links": [
            {"from":"A+", "to":"B+", "reads":10},
            {"from":"B+", "to":"C+", "reads":10},
            {"from":"C+", "to":"A+", "reads":10}
        ]});
        std::fs::write(dir.join("e.json"), ev.to_string()).unwrap();
        let params = LinearizeParams { min_reads: 2, max_steps: 0, adaptive_steps: false, max_molecules: 0,
            max_alternatives: 10, optional_below: 0.5 };
        let report = run(&dir.join("g.gfa"), &dir.join("e.json"), params,
            &dir.join("m.fa"), None, &dir.join("l.json")).unwrap();
        assert!(report.search_truncated);
        assert!(report.components.iter().all(|c| c.searched_max_molecules == 1));
        assert!(report.components.iter().all(|c| !c.decisive));
        // no cover was found, but the read-joined paths are one: A B C, closed by C to A
        let c = &report.components[0];
        assert!(c.read_join_cover && c.partial.is_none());
        assert_eq!((report.molecules, report.unsolved_components), (1, 0));
        assert!(report.best.as_ref().unwrap().circular);
        // the same paths of 24 bp are stubs: partial, the component unsolved
        std::fs::write(dir.join("g.gfa"),
            "H\tVN:Z:1.0\nS\tA\tACGTACGT\nS\tB\tGCTAGCTA\nS\tC\tTGCATGCA\nL\tA\t+\tB\t+\t0M\nL\tB\t+\tC\t+\t0M\nL\tC\t+\tA\t+\t0M\n").unwrap();
        let report = run(&dir.join("g.gfa"), &dir.join("e.json"), params,
            &dir.join("m.fa"), None, &dir.join("l.json")).unwrap();
        let c = &report.components[0];
        assert!(!c.read_join_cover && c.partial.is_some());
        assert_eq!((report.molecules, report.unsolved_components), (0, 1));
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn automatic_retry_is_bounded_and_explicit_budget_is_respected() {
        let dir = std::env::temp_dir().join(format!("ovasm-adaptive-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("g.gfa"),
            "H\tVN:Z:1.0\nS\tA\tACGTACGT\nS\tB\tGCTAGCTA\nS\tC\tTGCATGCA\nL\tA\t+\tB\t+\t0M\nL\tB\t+\tC\t+\t0M\nL\tC\t+\tA\t+\t0M\n").unwrap();
        let ev = json!({"repeat_segments": [], "phased_paths": [], "links": [
            {"from":"A+", "to":"B+", "reads":10},
            {"from":"B+", "to":"C+", "reads":10},
            {"from":"C+", "to":"A+", "reads":10}
        ]});
        std::fs::write(dir.join("e.json"), ev.to_string()).unwrap();
        let invoke = |steps, adaptive| run(&dir.join("g.gfa"), &dir.join("e.json"),
            LinearizeParams { min_reads: 2, max_steps: steps, adaptive_steps: adaptive,
                max_molecules: 1, max_alternatives: 10, optional_below: 0.5 },
            &dir.join("m.fa"), None, &dir.join("l.json")).unwrap();
        let mut checked = false;
        for steps in 1..100 {
            let explicit = invoke(steps, false);
            let automatic = invoke(steps, true);
            assert!(!explicit.components[0].adaptive_retry);
            if explicit.search_truncated && !automatic.search_truncated {
                assert!(automatic.components[0].adaptive_retry);
                checked = true;
                break;
            }
        }
        assert!(checked, "a bounded retry must finish a search near its first limit");
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn an_empty_graph_is_reported_as_empty_not_as_all_repeats() {
        let dir = std::env::temp_dir().join(format!("ovasm-empty-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("g.gfa"), "H\tVN:Z:1.0\n").unwrap();
        std::fs::write(
            dir.join("e.json"),
            r#"{"repeat_segments":[],"phased_paths":[],"links":[]}"#,
        )
        .unwrap();
        let params = LinearizeParams {
            min_reads: 2,
            max_steps: 1000, adaptive_steps: false,
            max_molecules: 4,
            max_alternatives: 10,
            optional_below: 0.5,
        };
        let err = run(
            &dir.join("g.gfa"),
            &dir.join("e.json"),
            params,
            &dir.join("m.fa"),
            None,
            &dir.join("l.json"),
        )
        .err()
        .unwrap()
        .to_string();
        std::fs::remove_dir_all(&dir).ok();
        assert!(err.contains("no segments"), "{err}");
    }

    fn model_evidence(gfa: &str, repeats: &[&str], links: &[(&str,&str,u64)], paths: &[(&str,u64)]) -> Value {
        let file=std::env::temp_dir().join(format!("ovasm-model-fingerprint-{}-{:?}.gfa",std::process::id(),std::thread::current().id()));
        std::fs::write(&file,gfa).unwrap();
        let graph=parse_gfa(&file).unwrap();std::fs::remove_file(file).unwrap();
        let phased:Vec<_>=paths.iter().filter(|(p,_)| {
            let tokens:Vec<_>=p.split_whitespace().collect();
            tokens.len()>2 && tokens[1..tokens.len()-1].iter().all(|t|repeats.contains(&&t[..t.len()-1]))
        }).map(|(p,n)|json!({"path":p,"reads":n})).collect();
        json!({"repeat_segments":repeats,"phased_paths":phased,
            "links":links.iter().map(|(a,b,n)|json!({"from":a,"to":b,"reads":n})).collect::<Vec<_>>(),
            "thread_evidence":{"schema":crate::thread_evidence::SCHEMA,
                "graph_fingerprint":crate::unzip::graph_fingerprint(&graph),
                "accepted_records":paths.iter().map(|(_,n)|n).sum::<u64>(),
                "links":links.iter().map(|(a,b,n)|json!({"from":a,"to":b,"records":n})).collect::<Vec<_>>(),
                "paths":paths.iter().map(|(p,n)|json!({"path":p,"reads":n})).collect::<Vec<_>>()}})
    }

    fn terminal_evidence(gfa: &str, links: &[(&str, &str, u64)], paths: &[(&str, u64)]) -> Value {
        let mut ev = model_evidence(gfa, &[], links, paths);
        ev["segment_depths"] = json!({"X":10,"A":100,"B":100,"Y":10});
        ev["median_depth"] = json!(100);
        ev
    }

    fn terminal_covers(report: &LinearizeReport) -> impl Iterator<Item = &CoverReport> {
        let component = &report.components[0];
        component.best.iter().chain(&component.alternatives)
    }

    #[test]
    fn native_model_represents_whole_threads_after_required_anchors() {
        let gfa = gfa_text(&[("X",100),("A",100),("B",100),("Y",100)],
                           &[("X+","A+"),("A+","B+"),("B+","Y+")]);
        let links = [("X+","A+",2),("A+","B+",2),("B+","Y+",2)];
        let cases = [
            ("ovasm-ends-forward", "X+ A+ B+ Y+", 2, true),
            ("ovasm-ends-reverse", "Y- B- A- X-", 2, true),
            ("ovasm-ends-low", "X+ A+ B+ Y+", 1, false),
        ];
        for (name, path, n, expected) in cases {
            let result = run_on(&gfa, &terminal_evidence(&gfa, &links, &[(path,n)]), name);
            let full = terminal_covers(&result).flat_map(|c| &c.molecules)
                .filter(|m| m.length == 400).collect::<Vec<_>>();
            assert_eq!(full.len(), usize::from(expected), "{name}");
            assert!(full.iter().all(|m| !m.circular));
            assert!(!result.physical_topology_identified);
            let model = result.components[0].read_path_cover.as_ref().unwrap();
            assert_eq!(model.represented_thread_records, if expected {2} else {0}, "{name}");
        }
    }

    #[test]
    fn native_model_does_not_promote_a_record_prefix_to_a_terminal_candidate() {
        let gfa = gfa_text(&[("X",100),("A",100),("B",100),("Y",100),("Z",100)],
                           &[("X+","A+"),("A+","B+"),("B+","Y+"),("Y+","Z+")]);
        let links = [("X+","A+",2),("A+","B+",2),("B+","Y+",2),("Y+","Z+",2)];
        let mut ev = terminal_evidence(&gfa, &links, &[("X+ A+ B+ Y+ Z+",2)]);
        ev["segment_depths"]["Z"] = json!(10);
        let mut old = terminal_evidence(&gfa, &links, &[]);
        old["segment_depths"]["Z"] = json!(10);
        let baseline = run_on(&gfa, &old, "ovasm-ends-prefix-baseline");
        let complete = run_on(&gfa, &ev, "ovasm-ends-prefix-complete");
        let model = complete.components[0].read_path_cover.as_ref().unwrap();
        assert_eq!(model.represented_thread_records, 2);
        let paths = |r: &LinearizeReport| terminal_covers(r).flat_map(|c| &c.molecules)
            .map(|m| m.path.clone()).collect::<BTreeSet<_>>();
        let before = paths(&baseline);
        let after = paths(&complete);
        let added: Vec<_> = after.difference(&before).collect();
        assert_eq!(added.len(), 1);
        assert_eq!(added[0].split_whitespace().count(), 5, "{}", added[0]);
    }

    #[test]
    fn native_model_terminal_budget_keeps_the_original_cover_search() {
        let gfa = gfa_text(&[("X",100),("A",100),("B",100),("Y",100)],
                           &[("X+","A+"),("A+","B+"),("B+","Y+")]);
        let links = [("X+","A+",2),("A+","B+",2),("B+","Y+",2)];
        let ev = terminal_evidence(&gfa, &links, &[("X+ A+ B+ Y+",2)]);
        let old = terminal_evidence(&gfa, &links, &[]);
        let file = std::env::temp_dir().join(format!("ovasm-ends-budget-{}.gfa",std::process::id()));
        std::fs::write(&file, &gfa).unwrap();
        let graph = parse_gfa(&file).unwrap();
        std::fs::remove_file(file).unwrap();
        let anchors: Vec<_> = (0..graph.names.len()).collect();
        let optional: Vec<_> = graph.names.iter().map(|n| n=="X" || n=="Y").collect();
        let repeats = vec![false; anchors.len()];
        let before = crate::read_model::ReadModel::load(&graph, &old, &repeats, 2).unwrap()
            .unwrap().thread_ends(&anchors,&optional);
        let after = crate::read_model::ReadModel::load(&graph, &ev, &repeats, 2).unwrap()
            .unwrap().thread_ends(&anchors,&optional);
        let (bridges, _) = bridges(&graph, &ev, &repeats, 2).unwrap();
        let mut out = FxHashMap::<_,Vec<_>>::default();
        for (i,b) in bridges.iter().enumerate() {
            out.entry(b.walk[0]).or_default().push(i);
        }
        let phases = crate::phase::Constraints::default();
        let solve = |budget, threads| search_phased_indexed(
            anchors.len(),&anchors,&optional,&bridges,1,budget,&phases,&out,Some(threads));
        let budget = (1..200).find(|&n| !solve(n,&before).1).unwrap();
        let (original, truncated, terminal_truncated) = solve(budget,&before);
        assert!(!truncated && !terminal_truncated && !original.is_empty());
        let (extended, truncated, terminal_truncated) = solve(budget,&after);
        assert!(!truncated && terminal_truncated);
        let paths = |covers: Vec<Vec<Found>>| covers.into_iter().map(|c| c.into_iter()
            .map(|f| (f.circular, candidate_walk(&bridges,&f))).collect::<Vec<_>>()).collect::<Vec<_>>();
        assert_eq!(paths(original), paths(extended));
    }

    #[test]
    fn native_model_terminal_threads_do_not_add_circles() {
        let gfa = gfa_text(&[("X",100),("A",100),("B",100),("Y",100)],
                           &[("X+","A+"),("A+","B+"),("B+","Y+"),("Y+","X+")]);
        let links = [("X+","A+",2),("A+","B+",2),("B+","Y+",2),("Y+","X+",2)];
        let baseline = run_on(&gfa, &terminal_evidence(&gfa, &links, &[]), "ovasm-ends-circle-base");
        let extended = run_on(&gfa, &terminal_evidence(&gfa, &links, &[("X+ A+ B+ Y+",2)]),
                              "ovasm-ends-circle-read");
        let circles = |r: &LinearizeReport| terminal_covers(r).flat_map(|c| &c.molecules)
            .filter(|m| m.circular).map(|m| m.path.clone()).collect::<BTreeSet<_>>();
        assert_eq!(circles(&baseline), circles(&extended));
        assert!(!extended.physical_topology_identified);
    }

    #[test]
    fn effective_model_distinguishes_fallback_mixed_and_no_candidates() {
        assert_eq!(super::effective_ranking_model([].into_iter()), "no_ranked_components");
        assert_eq!(super::effective_ranking_model(["legacy_after_read_budget_exhaustion"].into_iter()), "legacy_after_read_budget_exhaustion");
        assert_eq!(super::effective_ranking_model([crate::read_model::MODEL, "legacy_after_read_budget_exhaustion"].into_iter()), "mixed_component_models");
        assert_eq!(super::effective_ranking_model([crate::read_model::MODEL].into_iter()), crate::read_model::MODEL);
        assert_eq!(super::effective_ranking_model(["legacy_bridge_copy"].into_iter()), "legacy_bridge_copy");
    }

    #[test]
    fn native_model_uses_observed_closing_linkage() {
        let g=gfa_text(&[("A",100),("B",100),("C",100)],&[("A+","B+"),("B+","C+"),("C+","A+")]);
        let ev=model_evidence(&g,&[],&[("A+","B+",20),("B+","C+",20),("C+","A+",10)],&[("A+ B+ C+",10),("B+ C+ A+",10)]);
        let r=run_on(&g,&ev,"ovasm-native-model-closed");
        assert_eq!(r.ranking_model,crate::read_model::MODEL);
        assert!(r.best.unwrap().circular); assert!(!r.cover_size_limited);
        let c=&r.components[0];assert!(c.best.as_ref().unwrap().read_support.as_ref().unwrap().complete());
        assert!(c.margin_ln.is_none());assert!(!r.physical_topology_identified);
    }

    #[test]
    fn native_model_does_not_turn_weak_closure_into_certainty() {
        let g=gfa_text(&[("A",100),("B",100),("C",100)],&[("A+","B+"),("B+","C+"),("C+","A+")]);
        let ev=model_evidence(&g,&[],&[("A+","B+",10),("B+","C+",10),("C+","A+",1)],&[("A+ B+ C+",10)]);
        let r=run_on(&g,&ev,"ovasm-native-model-weak-closed");
        assert!(!r.best.unwrap().circular);
        assert!(r.components[0].evidence_tied_candidates>=2);
        assert!(!r.decisive);
        assert!(r.components[0].alternatives.iter().any(|c|c.molecules.iter().any(|m|m.circular)));
    }

    #[test]
    fn native_model_does_not_preselect_the_loudest_closing_bridge() {
        let g=gfa_text(&[("A",100),("R",50),("S",50)],&[("A+","R+"),("R+","A+"),("A+","S+"),("S+","A+")]);
        let mut ev=model_evidence(&g,&["R","S"],&[("A+","R+",5),("R+","A+",5),("A+","S+",10),("S+","A+",10)],&[("A+ R+ A+",5),("A+ S+ A+",10)]);
        ev["phased_paths"]=json!([{"path":"A+ R+ A+","reads":100},{"path":"A+ S+ A+","reads":20}]);
        let r=run_on(&g,&ev,"ovasm-native-model-all-closings");
        assert!(r.best.as_ref().unwrap().path.contains('S'));
        assert!(r.components[0].alternatives.iter().any(|c|c.molecules.iter().any(|m|m.path.contains('R'))));
        assert!(!r.decisive,"coexisting read-supported alternatives remain unresolved");
        assert_eq!(r.components[0].witness_configurations.len(),2);
        let family=r.components[0].read_path_cover.as_ref().unwrap();
        assert_eq!((family.supported_thread_records,family.represented_thread_records),(15,15));
        assert!(family.minimum_family_size_proven);
    }
    #[test]
    fn native_model_keeps_component_evidence_scopes_separate() {
        let g=gfa_text(&[("A",100),("B",100),("C",100),("D",100),("E",100),("F",100)],
            &[("A+","B+"),("B+","C+"),("C+","A+"),("D+","E+"),("E+","F+"),("F+","D+")]);
        let ev=model_evidence(&g,&[],&[("A+","B+",20),("B+","C+",20),("C+","A+",10),
            ("D+","E+",20),("E+","F+",20),("F+","D+",10)],
            &[("A+ B+ C+",10),("B+ C+ A+",10),("D+ E+ F+",10),("E+ F+ D+",10)]);
        let r=run_on(&g,&ev,"ovasm-native-model-components");
        assert_eq!(r.components.len(),2);
        for c in &r.components {assert_eq!(c.best.as_ref().unwrap().read_support.as_ref().unwrap().unexplained_thread_records,0);}
    }
}
