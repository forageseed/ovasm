//! Stage 2b: unringing, representative molecules from read-supported paths (`ovasm linearize`).
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
//! its first. The fewest molecules that cover a component win (a master circle or path before
//! sub-circles); among covers of that size, `cover_order` ranks them by circularity and by
//! `path_score`: the likelihood of the paths when each anchor's exit is chosen in proportion to
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
use std::path::Path;

use anyhow::{bail, Context, Result};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use serde_json::Value;

use crate::evidence::{flip, parse_gfa, Graph, Oriented};

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
    /// Most molecules one component may be split into before it is reported unsolved.
    pub max_molecules: usize,
    pub max_alternatives: usize,
    /// An anchor shallower than this fraction of the median depth is optional: a low-frequency
    /// configuration the representative path need not visit. In a component shallower than
    /// the graph's anchors as a whole the threshold comes down in proportion, unless the
    /// component is itself below the fraction (see `run`). Without depth in the evidence,
    /// anchors below half this fraction of the median bridge support are optional.
    pub optional_below: f64,
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
        let reads = l["reads"].as_u64().unwrap_or(0);
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
/// as `ovasm evidence` takes it over the whole graph.
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
    let mut out_of: FxHashMap<Oriented, Vec<usize>> = FxHashMap::default();
    for (i, b) in br.iter().enumerate() {
        out_of.entry(b.walk[0]).or_default().push(i);
    }
    struct Ctx<'a> {
        out_of: &'a FxHashMap<Oriented, Vec<usize>>,
        br: &'a [Bridge],
        optional: &'a [bool],
        anchors: &'a [usize],
        want: usize,
        m: usize,
        max_steps: u64,
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
        if s.covered < c.want {
            for &i in c.out_of.get(&at).map(Vec::as_slice).unwrap_or(&[]) {
                let next = *c.br[i].walk.last().unwrap();
                if s.used[next.0] {
                    continue;
                }
                s.used[next.0] = true;
                s.covered += usize::from(!c.optional[next.0]);
                s.path.push(i);
                dfs(c, s, start, next);
                s.path.pop();
                s.covered -= usize::from(!c.optional[next.0]);
                s.used[next.0] = false;
            }
        }
        // End the molecule here.
        if !visits(c, start, &s.path, s.must) {
            return;
        }
        let closing = c.out_of.get(&at).and_then(|v| {
            v.iter()
                .copied()
                .filter(|&i| *c.br[i].walk.last().unwrap() == start)
                .max_by_key(|&i| c.br[i].reads)
        });
        // A path ends only where it cannot go on: no bridge from either end reaches a required
        // anchor still unused. Cutting it anywhere else just multiplies covers.
        let open_to = |o: Oriented| {
            c.out_of.get(&o).is_some_and(|v| {
                v.iter().any(|&i| {
                    let a = c.br[i].walk.last().unwrap().0;
                    !s.used[a] && !c.optional[a]
                })
            })
        };
        if closing.is_none() && (open_to(at) || open_to(flip(start))) {
            return;
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
                s.found.push(cover);
            }
            return;
        }
        if s.done.len() + 1 == c.m {
            return;
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
            }
        }
        s.done.pop();
        s.path = saved_path;
        s.must = saved_must;
    }

    let ctx = Ctx {
        out_of: &out_of,
        br,
        optional,
        anchors,
        want: anchors.iter().filter(|&&a| !optional[a]).count(),
        m,
        max_steps,
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
        }
    }
    (st.found, st.truncated)
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
        parent[ra] = rb;
        exit.insert(a, i);
        exit.insert(flip(b), pos[rc_walk(&br[i].walk).as_slice()]);
    }
    let closing = |from: Oriented, to: Oriented| {
        (0..br.len())
            .filter(|&i| !br[i].topological && ends(i) == (from, to))
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

/// Sequence of a walk, each link's overlap written once. A circular walk also drops the overlap
/// of the closing link from its end.
fn spell(g: &Graph, walk: &[Oriented], circular: bool) -> Result<Vec<u8>> {
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
    /// Repeat copies laid down more or fewer times than the repeats' depths say (see
    /// `copy_deviation`); each costs `COPY_PENALTY_LN` in the ranking.
    pub copy_deviation: usize,
}

impl CoverReport {
    fn linear(&self) -> usize {
        self.molecules.iter().filter(|m| !m.circular).count()
    }

    /// Score the ranking uses: the likelihood less the copy-number penalty.
    fn ranked_score(&self) -> f64 {
        self.score - COPY_PENALTY_LN * self.copy_deviation as f64
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

/// Order of two covers of one component with the same number of molecules; `Less` ranks first.
/// Fewer linear molecules first: a linear molecule lacks its closing bridge, whose ln share is
/// never positive, so scores alone would favour opening circles. Then the likelier cover.
fn cover_order(a: &CoverReport, b: &CoverReport) -> Ordering {
    a.linear()
        .cmp(&b.linear())
        .then(b.ranked_score().total_cmp(&a.ranked_score()))
}

/// Anchors joined by bridges, and the molecules unringing finds for them.
#[derive(Serialize)]
pub struct ComponentReport {
    pub anchors: Vec<String>,
    /// Length-weighted median depth of the anchors, when the evidence records depth.
    pub depth: Option<f64>,
    /// Fewest molecules that cover the component, plus the tails cut off (0 when unsolved or
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
    pub best: Option<CoverReport>,
    /// For a component without a cover: the paths its read bridges make (see `partial_paths`),
    /// so that its sequence is reported. Not a solution: `best` stays empty, `molecules` 0, and
    /// the component counts among `unsolved_components`.
    pub partial: Option<CoverReport>,
    /// Log-likelihood lead of the best cover over the runner-up with as many linear molecules.
    pub margin_ln: Option<f64>,
    pub decisive: bool,
    /// Bridges left out where the best cover's linear molecule led back into itself and is
    /// reported as a circle and its tail (see `cut_tail`); reads support these bridges too, so
    /// a component with one is not decisive. `alternatives` and `margin_ln` refer to the
    /// covers before the cut.
    pub tails_cut: Vec<String>,
    pub alternatives: Vec<CoverReport>,
}

#[derive(Serialize)]
pub struct LinearizeReport {
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
    /// The longest molecule; `alternatives` and `margin_ln` refer to its component.
    pub best: Option<Linearization>,
    /// Log-likelihood lead of the best path over the runner-up (None with a single candidate).
    pub margin_ln: Option<f64>,
    /// Whether every component's best cover is at least twice as likely as its runner-up. When
    /// false, the reads support several configurations about equally and the molecules are one
    /// representative of them.
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
    let g = parse_gfa(gfa)?;
    let report: Value = serde_json::from_str(
        &fs::read_to_string(evidence)
            .with_context(|| format!("cannot read {}", evidence.display()))?,
    )?;
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
    if g.names.is_empty() {
        bail!("the assembly graph has no segments; the reads gave nothing to assemble");
    }
    if anchors.is_empty() {
        bail!("every segment is a repeat; nothing anchors a linear path");
    }
    let (br, self_loops) = bridges(&g, &report, &repeat, p.min_reads)?;
    let support = |a: usize| -> u64 {
        br.iter()
            .filter(|b| b.walk[0].0 == a || b.walk[b.walk.len() - 1].0 == a)
            .map(|b| b.reads)
            .sum()
    };
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
    // copies each repeat should be laid down, from its depth
    let expected: Vec<Option<usize>> = match (
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
    let frac = exit_fractions(&br);
    let label = |walk: &[Oriented]| {
        walk.iter()
            .map(|&o| g.label(o))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let describe = |f: &Found| -> Result<(Vec<Oriented>, Linearization)> {
        let walk = candidate_walk(&br, f);
        let fractions: Vec<f64> = f.bridges.iter().map(|&i| frac[i]).collect();
        let lin = Linearization {
            path: label(&walk),
            circular: f.circular,
            length: spell(&g, &walk, f.circular)?.len(),
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
            copy_deviation: copy_deviation(&expected, &br, cover),
        };
        Ok((walks, report))
    };

    let mut comps: Vec<ComponentReport> = Vec::new();
    // every component's best molecules, for spelling and placement checks
    let mut chosen: Vec<(Vec<Oriented>, Linearization)> = Vec::new();
    // what no cover explains, as partial paths
    let mut partial: Vec<(Vec<Oriented>, Linearization)> = Vec::new();
    for (comp, depth) in anchor_components.into_iter().zip(component_depth) {
        let mut rep = ComponentReport {
            anchors: comp.iter().map(|&a| g.names[a].clone()).collect(),
            depth,
            molecules: 0,
            path_ends: Vec::new(),
            min_molecules: 0,
            skipped: comp.iter().all(|&a| optional[a]),
            candidates: 0,
            search_truncated: false,
            best: None,
            partial: None,
            margin_ln: None,
            decisive: false,
            tails_cut: Vec::new(),
            alternatives: Vec::new(),
        };
        if rep.skipped {
            comps.push(rep);
            continue;
        }
        let bridged_at = |o: Oriented| br.iter().any(|b| b.walk[0] == o);
        let ends: Vec<usize> = comp
            .iter()
            .copied()
            .filter(|&a| !optional[a] && bridged_at((a, true)) != bridged_at((a, false)))
            .collect();
        rep.path_ends = ends.iter().map(|&a| g.names[a].clone()).collect();
        rep.min_molecules = ends.len().div_ceil(2).max(1);
        // fewest molecules first: stop at the first cover size that has any cover
        let mut ranked: Vec<(Vec<Vec<Oriented>>, CoverReport, Vec<Found>)> = Vec::new();
        for m in rep.min_molecules..=p.max_molecules.max(1) {
            let (found, truncated) = search(g.names.len(), &comp, &optional, &br, m, p.max_steps);
            rep.search_truncated |= truncated;
            let mut seen = BTreeSet::new();
            for cover in found {
                let (walks, report) = cover_report(&cover)?;
                let mut key: Vec<(bool, Vec<Oriented>)> = walks
                    .iter()
                    .zip(&report.molecules)
                    .map(|(w, l)| (l.circular, canonical_key(w, l.circular)))
                    .collect();
                key.sort();
                if seen.insert(key) {
                    ranked.push((walks, report, cover));
                }
            }
            if !ranked.is_empty() {
                rep.molecules = m;
                break;
            }
        }
        if ranked.is_empty() {
            let paths = partial_paths(g.names.len(), &comp, &optional, &br, &frac);
            let (walks, report) = cover_report(&paths)?;
            partial.extend(walks.into_iter().zip(report.molecules.iter().cloned()));
            rep.partial = Some(report);
        }
        ranked.sort_by(|(wa, a, _), (wb, b, _)| cover_order(a, b).then_with(|| wa.cmp(wb)));
        rep.candidates = ranked.len();
        rep.margin_ln = match (ranked.first(), ranked.get(1)) {
            (Some((_, a, _)), Some((_, b, _))) if a.linear() == b.linear() => {
                Some(a.ranked_score() - b.ranked_score())
            }
            _ => None,
        };
        // a cover through a pairing read fewer than `min_reads` times is not decisive
        rep.decisive = ranked.first().is_some_and(|(_, c, _)| c.weak_bridges == 0)
            && rep.margin_ln.map_or(true, |m| m >= DECISIVE_LN);
        let mut it = ranked.into_iter();
        if let Some((mut walks, mut best, cover)) = it.next() {
            // a linear molecule that leads back into itself: the circle and its tail
            let mut cut: Vec<Found> = Vec::new();
            for f in &cover {
                match cut_tail(&g, &br, f) {
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
            chosen.extend(walks.into_iter().zip(best.molecules.iter().cloned()));
            rep.best = Some(best);
        }
        rep.alternatives = it.take(p.max_alternatives).map(|(_, c, _)| c).collect();
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
            anchors: vec![g.names[a].clone()],
            depth: None,
            molecules: 1,
            path_ends: vec![g.names[a].clone()],
            min_molecules: 1,
            skipped: false,
            candidates: 1,
            search_truncated: false,
            best: Some(CoverReport {
                score: 0.0,
                molecules: vec![lin],
                weak_bridges: 0,
                copy_deviation: 0,
            }),
            partial: None,
            margin_ln: None,
            decisive: false,
            tails_cut: Vec::new(),
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
    let partial_fasta =
        out_partial.map_or_else(|| partial_fasta_path(out_fasta), Path::to_path_buf);
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
    let out = LinearizeReport {
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
            copy_deviation: 0,
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

    /// Plastid LSC (u0), SSC (u1) and one IR copy (u2), as `ovasm assemble` collapses them,
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
            max_steps: 1_000_000,
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
            max_steps: 1_000_000,
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
    fn run_graph(gfa: &str, ev: Value, dir_name: &str) -> LinearizeReport {
        let dir = std::env::temp_dir().join(format!("{dir_name}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("g.gfa"), gfa).unwrap();
        std::fs::write(dir.join("e.json"), ev.to_string()).unwrap();
        let params = LinearizeParams {
            min_reads: 2,
            max_steps: 1_000_000,
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
    fn an_anchor_between_both_ends_of_an_unspanned_repeat_joins_the_circle() {
        // Buddleja mitochondrion: u1 sits between both ends of repeat u4, which no read spans
        // into u1; reads cross u4 only from u3 to u2. The genome is u0 u3 u4 u1 u4 u2.
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
        let best = r.best.as_ref().unwrap();
        assert!(best.circular, "{}", best.path);
        assert_eq!(
            best.length,
            400 + 200 + 60 + 300 + 60 + 100,
            "{}",
            best.path
        );
        assert!(best.path.contains("u1"));
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
            max_steps: 1000,
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
}
