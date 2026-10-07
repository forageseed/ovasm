//! Stage 4 (S4.1): HiFi assembly on a closed-syncmer sparse de Bruijn graph (`ovasm assemble`).
//!
//! Only closed-syncmer k-mers are nodes: a k-mer is sampled when the smallest (hashed, canonical)
//! s-mer inside it sits at its first or its last position. Every window of `k - s + 1` k-mer
//! starts holds at least one, so consecutive sampled k-mers of a read are at most `k - s` apart
//! and overlap by at least `s` bases: the graph spells exact sequence without a consensus step,
//! at a density of about `2 / (k - s + 1)`. Nodes are identified by a hash of the whole k-mer, so
//! a sequencing error anywhere in it gives a different, rarely seen node that the count
//! threshold removes. Edges join consecutive sampled k-mers of a read (with their distance);
//! compacting non-branching runs gives unitigs, written as GFA with exact `(k - d)M` overlaps.
//!
//! With k shorter than the organelle's repeats, each repeat collapses into a branching unitig,
//! which is the graph shape `ovasm evidence` and `ovasm linearize` resolve with reads.

use std::collections::{BTreeSet, VecDeque};
use std::fmt::Write as _;
use std::fs;
use std::hash::{Hash, Hasher};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet, FxHasher};
use serde::Serialize;

use crate::evidence::revcomp;
use crate::kmer::CanonicalKmers;
use crate::recruit::{peak_rss_mb, read_batches};

/// Oriented node: (node index, forward?).
type ONode = (usize, bool);
/// Canonical edge between two oriented k-mer hashes, with the distance between their starts.
type EdgeKey = ((u64, bool), (u64, bool), u32);
/// A sampled k-mer in a read: hash, forward (canonical) strand?, start in the read.
type Sample = (u64, bool, u32);

#[derive(Clone, Copy)]
pub struct AssembleParams {
    pub k: usize,
    pub s: usize,
    /// Minimum reads for a node or an edge; `None` picks it from the node-count distribution.
    pub min_count: Option<u32>,
}

#[inline]
fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn hash_bytes(b: &[u8]) -> u64 {
    let mut h = FxHasher::default();
    h.write(b);
    mix64(h.finish())
}

/// Starts of closed-syncmer k-mers in `seq` (uppercase ACGT runs only).
///
/// Closed: the minimum hashed canonical s-mer of the k-mer equals the value at its first or its
/// last s-mer position. Comparing values (not positions) keeps the rule strand-symmetric.
fn closed_syncmers(seq: &[u8], k: usize, s: usize) -> Vec<usize> {
    let w = k - s + 1;
    let mut out = Vec::new();
    let mut run: Vec<u64> = Vec::new();
    let mut run_start = 0usize;
    let flush = |run: &mut Vec<u64>, start: usize, out: &mut Vec<usize>| {
        if run.len() >= w {
            let mut dq: VecDeque<usize> = VecDeque::new();
            for i in 0..run.len() {
                while dq.back().is_some_and(|&j| run[j] >= run[i]) {
                    dq.pop_back();
                }
                dq.push_back(i);
                if i + 1 >= w {
                    let first = i + 1 - w;
                    while dq.front().is_some_and(|&j| j < first) {
                        dq.pop_front();
                    }
                    let min = run[*dq.front().unwrap()];
                    if run[first] == min || run[i] == min {
                        out.push(start + first);
                    }
                }
            }
        }
        run.clear();
    };
    let mut last: Option<usize> = None;
    for (pos, smer) in CanonicalKmers::new(seq, s) {
        if last.is_some_and(|l| pos != l + 1) {
            flush(&mut run, run_start, &mut out);
        }
        if run.is_empty() {
            run_start = pos;
        }
        run.push(mix64(smer));
        last = Some(pos);
    }
    flush(&mut run, run_start, &mut out);
    out
}

/// (canonical k-mer hash, read shows the canonical strand?, start) for each sampled k-mer.
fn sample_read(seq: &[u8], p: &AssembleParams) -> Vec<Sample> {
    if p.s == 0 {
        return sample_dense(seq, p.k);
    }
    closed_syncmers(seq, p.k, p.s)
        .into_iter()
        .map(|pos| {
            let km = &seq[pos..pos + p.k];
            let (hf, hr) = (hash_bytes(km), hash_bytes(&revcomp(km)));
            (hf.min(hr), hf <= hr, pos as u32)
        })
        .collect()
}

/// Default k in dense mode: the median read length (first 10,000 reads of the first input)
/// minus 23, odd, at least 21. Each read then gives 24 k-mers: Nipponbare 150 bp reads at
/// ~1,150x gave k 127 (one molecule); Salvia 100 bp reads at ~150x gave k 77 (every junction),
/// where k 99 and 127 left no graph at all.
pub fn dense_k(inputs: &[PathBuf]) -> Result<(usize, usize)> {
    let first = inputs.first().context("no read file")?;
    let mut parser = needletail::parse_fastx_file(first)
        .with_context(|| format!("cannot read {}", first.display()))?;
    let mut lens = Vec::new();
    while let Some(r) = parser.next() {
        lens.push(r?.num_bases());
        if lens.len() == 10_000 {
            break;
        }
    }
    if lens.is_empty() {
        bail!("{} holds no reads", first.display());
    }
    lens.sort_unstable();
    let median = lens[lens.len() / 2];
    let k = median.saturating_sub(23).max(21);
    Ok((if k % 2 == 0 { k - 1 } else { k }.max(21), median))
}

/// Every k-mer of a read (dense mode, `s` = 0), skipping those with a non-ACGT base. Short
/// reads need it: sampled k-mers become nodes and neighbours in one read become edges, and
/// syncmers leave up to k - s between samples (s <= 31), so a 150 bp read would join samples
/// only for k below about 90.
fn sample_dense(seq: &[u8], k: usize) -> Vec<Sample> {
    let mut out = Vec::with_capacity(seq.len().saturating_sub(k) + 1);
    let mut rc = Vec::with_capacity(k);
    let mut run = 0usize;
    for (i, &b) in seq.iter().enumerate() {
        run = if matches!(b, b'A' | b'C' | b'G' | b'T') {
            run + 1
        } else {
            0
        };
        if run >= k {
            let pos = i + 1 - k;
            let km = &seq[pos..=i];
            rc.clear();
            rc.extend(km.iter().rev().map(|&c| match c {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            }));
            let (hf, hr) = (hash_bytes(km), hash_bytes(&rc));
            out.push((hf.min(hr), hf <= hr, pos as u32));
        }
    }
    out
}

#[derive(Serialize)]
pub struct AssembleReport {
    pub k: usize,
    /// Automatic k: each k tried, largest first, with its node depth (see `run_auto_k`).
    pub k_tried: Vec<(usize, f64)>,
    /// Automatic k found no k deep enough on the reads given, so it self-corrected them and
    /// assembled these (see `run_auto_k`); `k_tried` then lists both ladders.
    pub corrected_reads: Option<String>,
    pub s: usize,
    pub reads: u64,
    pub bases: u64,
    pub sampled_kmers: u64,
    pub distinct_nodes: usize,
    pub min_count: u32,
    pub min_count_auto: bool,
    /// Typical node depth (count-weighted median of nodes seen at least 3 times).
    pub node_depth: f64,
    pub nodes_kept: usize,
    pub edges_kept: usize,
    pub cleaning: CleanStats,
    /// Dead ends reconnected through k-mers the reads agree on (see `rescue_gaps`).
    pub rescue: RescueStats,
    /// Dense mode: dead ends joined by an exact overlap shorter than k - 1 (see
    /// `join_dead_ends`).
    pub overlap_joins: usize,
    /// Dense mode: gaps between dead ends filled by a local k = 31 assembly (see `fill_gaps`).
    pub gaps_filled: usize,
    /// Bases moved into branching unitigs by repeat-boundary normalization.
    pub boundary_bp_moved: usize,
    /// Links whose overlap was trimmed from the side that has no other link (now 0M).
    pub links_blunted: usize,
    /// Links still written with an overlap (both ends branch).
    pub links_with_overlap: usize,
    pub unitigs: usize,
    /// Unitigs whose last node continues into their first (a closed circle).
    pub circular_unitigs: usize,
    pub links: usize,
    pub total_length: usize,
    pub n50: usize,
    pub longest: usize,
    pub elapsed_seconds: f64,
    pub peak_rss_mb: Option<f64>,
}

struct Unitig {
    nodes: Vec<(ONode, u32)>, // oriented node and the distance from the previous one
    circular: bool,
    /// For a circle, the distance from the last node back to the first (else 0).
    close: u32,
}

fn flip(o: ONode) -> ONode {
    (o.0, !o.1)
}

type Adj = [Vec<(ONode, u32)>];
type AdjList = Vec<Vec<(ONode, u32)>>;
/// Canonical link between unitig ends ((from, forward?), (to, forward?)) and its overlap in bp.
type ULink = (((usize, bool), (usize, bool)), usize);
/// Unitig end walked forward (true) or reverse (false).
type OEnd = (usize, bool);
/// Link seen from one end: (oriented unitig at the other end, overlap).
type OLink = ((usize, bool), usize);

fn live_succ(adj: &Adj, alive: &[bool], o: ONode) -> Vec<(ONode, u32)> {
    adj[2 * o.0 + o.1 as usize]
        .iter()
        .copied()
        .filter(|(n, _)| alive[n.0])
        .collect()
}

/// Maximal non-branching runs of live nodes.
fn compact(adj: &Adj, alive: &[bool]) -> Vec<Unitig> {
    let n = alive.len();
    // The unique successor of `o` when it also has `o` as its unique predecessor.
    let extend = |o: ONode| -> Option<(ONode, u32)> {
        match live_succ(adj, alive, o).as_slice() {
            [(m, d)] if m.0 != o.0 && live_succ(adj, alive, flip(*m)).len() == 1 => Some((*m, *d)),
            _ => None,
        }
    };
    let mut visited = vec![false; n];
    let mut unitigs = Vec::new();
    for i in (0..n).filter(|&i| alive[i]) {
        if visited[i] {
            continue;
        }
        let mut start: ONode = (i, true);
        let mut circular = false;
        let mut guard = 0;
        while let Some((pn, _)) = extend(flip(start)) {
            let prev = flip(pn);
            if prev.0 == i {
                circular = true;
                break;
            }
            start = prev;
            guard += 1;
            if guard > n {
                break;
            }
        }
        let mut path = vec![(start, 0u32)];
        visited[start.0] = true;
        let mut cur = start;
        let mut close = 0;
        while let Some((m, d)) = extend(cur) {
            if m == start {
                circular = true;
                close = d;
                break;
            }
            if visited[m.0] {
                break;
            }
            visited[m.0] = true;
            path.push((m, d));
            cur = m;
        }
        unitigs.push(Unitig {
            nodes: path,
            circular,
            close,
        });
    }
    unitigs
}

/// Links between unitig ends, one per adjacency.
fn unitig_links(adj: &Adj, alive: &[bool], unitigs: &[Unitig], k: usize) -> Vec<ULink> {
    // oriented node -> (unitig, unitig orientation) when leaving / entering the unitig there
    let mut leave: FxHashMap<ONode, (usize, bool)> = FxHashMap::default();
    let mut enter: FxHashMap<ONode, (usize, bool)> = FxHashMap::default();
    for (u, t) in unitigs.iter().enumerate() {
        let (head, tail) = (t.nodes[0].0, t.nodes[t.nodes.len() - 1].0);
        leave.insert(tail, (u, true));
        leave.insert(flip(head), (u, false));
        enter.insert(head, (u, true));
        enter.insert(flip(tail), (u, false));
    }
    let mut links = std::collections::BTreeSet::new();
    for (sl, succs) in adj.iter().enumerate() {
        let a: ONode = (sl / 2, sl % 2 == 1);
        if !alive[a.0] {
            continue;
        }
        for &(b, d) in succs.iter().filter(|(b, _)| alive[b.0]) {
            if let (Some(&from), Some(&to)) = (leave.get(&a), enter.get(&b)) {
                let rc = ((to.0, !to.1), (from.0, !from.1));
                let key = if (from, to) <= rc { (from, to) } else { rc };
                links.insert((key, k - d as usize));
            }
        }
    }
    links.into_iter().collect()
}

/// Whether the unitigs leading into `u` rejoin what lies beyond it by another route of at most
/// `budget` bp that avoids `u`: then `u` is a detour (a bubble arm, simple or nested). "Beyond"
/// is everything reachable from `u`'s exits within `budget`, not just the exits: an error arm is
/// often several unitigs (it branches again, e.g. into a second copy of a nearby repeat), and its
/// first unitig's exits are themselves error pieces that rejoin the true path a step later.
fn has_detour(
    u: usize,
    out: &FxHashMap<(usize, bool), Vec<(usize, bool)>>,
    len: &[usize],
    budget: usize,
) -> bool {
    let next = |o: (usize, bool)| out.get(&o).map(Vec::as_slice).unwrap_or(&[]);
    // Oriented unitigs reachable from `from` avoiding u, each with the cheapest cost of getting
    // there: the length of the unitigs passed on the way (not its own, so a long unitig where
    // the arms rejoin is still reached), at most `budget`.
    let reach = |from: &[(usize, bool)]| {
        let mut seen: FxHashMap<(usize, bool), usize> = FxHashMap::default();
        let mut stack: Vec<(usize, bool)> = Vec::new();
        for &o in from.iter().filter(|o| o.0 != u) {
            seen.insert(o, 0);
            stack.push(o);
        }
        while let Some(cur) = stack.pop() {
            let cost = seen[&cur] + len[cur.0];
            if cost > budget {
                continue;
            }
            for &nx in next(cur) {
                if nx.0 != u && seen.get(&nx).map_or(true, |&c| cost < c) {
                    seen.insert(nx, cost);
                    stack.push(nx);
                }
            }
        }
        seen
    };
    let beyond = reach(next((u, true)));
    if beyond.is_empty() {
        return false;
    }
    // entering u forward means leaving a predecessor p with p -> u+, i.e. flip of out[(u, -)]
    let entries: Vec<(usize, bool)> = next((u, false)).iter().map(|&(v, f)| (v, !f)).collect();
    // the entries themselves are where the arms split, not part of the other route
    let around: Vec<(usize, bool)> = entries
        .iter()
        .flat_map(|&e| next(e).iter().copied())
        .filter(|o| o.0 != u)
        .collect();
    reach(&around).keys().any(|o| beyond.contains_key(o))
}

/// Move sequence that every unitig entering a node carries into that node.
///
/// In a sparse graph a repeat's node holds only the sampled k-mers wholly inside the repeat; the
/// unitigs entering it end with up to ~k bases of the repeat each, joined by short overlaps. So
/// the branch point sits inside the repeat and the repeat's start is duplicated in every entry.
/// When all entries of a node end continue only into it, the bases they share just before the
/// node (aligned at the node's start, i.e. leaving out each entry's own overlap) move from each
/// entry into the node: every overlap stays valid, every walk spells the same sequence, and the
/// branch point moves to where the entries really diverge. Both
/// orientations of every unitig are visited, so exits are handled as the entries of the flip.
fn normalize_boundaries(seqs: &mut [Vec<u8>], links: &[ULink]) -> usize {
    let oriented = |s: &[u8], f: bool| if f { s.to_vec() } else { revcomp(s) };
    let mut out: FxHashMap<(usize, bool), Vec<OLink>> = FxHashMap::default();
    // largest overlap of any link touching each unitig: it must stay longer than that
    let mut max_ov = vec![0usize; seqs.len()];
    for &((a, b), ov) in links {
        out.entry(a).or_default().push((b, ov));
        out.entry((b.0, !b.1)).or_default().push(((a.0, !a.1), ov));
        max_ov[a.0] = max_ov[a.0].max(ov);
        max_ov[b.0] = max_ov[b.0].max(ov);
    }
    let mut moved = 0;
    for _ in 0..10 {
        let mut changed = false;
        for v in 0..seqs.len() {
            for vf in [true, false] {
                // entries of (v, vf): p with p -> (v, vf), i.e. flips of the exits of (v, !vf)
                let preds: Vec<OLink> = out
                    .get(&(v, !vf))
                    .map(|x| x.iter().map(|&(q, ov)| ((q.0, !q.1), ov)).collect())
                    .unwrap_or_default();
                if preds.len() < 2 {
                    continue;
                }
                let mut ids: Vec<usize> = preds.iter().map(|(q, _)| q.0).collect();
                ids.sort_unstable();
                ids.dedup();
                if ids.len() != preds.len()
                    || ids.contains(&v)
                    || preds
                        .iter()
                        .any(|(q, _)| out.get(q).map_or(0, Vec::len) != 1)
                {
                    continue;
                }
                // Each entry reaches `ov` bases into the node; align the entries at the node's
                // start by leaving those bases out, then take the suffix they all share.
                let qs: Vec<Vec<u8>> = preds
                    .iter()
                    .map(|&(q, ov)| {
                        let x = oriented(&seqs[q.0], q.1);
                        x[..x.len().saturating_sub(ov)].to_vec()
                    })
                    .collect();
                // keep at least half of every entry, and every entry longer than the overlaps
                // at its other end (an entry may lose bases at both ends)
                let cap = qs
                    .iter()
                    .zip(&preds)
                    .map(|(x, &(q, _))| {
                        (x.len() / 2).min(seqs[q.0].len().saturating_sub(max_ov[q.0] + 1))
                    })
                    .min()
                    .unwrap_or(0);
                let mut shift = 0;
                while shift < cap
                    && qs
                        .iter()
                        .all(|x| x[x.len() - 1 - shift] == qs[0][qs[0].len() - 1 - shift])
                {
                    shift += 1;
                }
                if shift == 0 {
                    continue;
                }
                // Moving the shared bases keeps every overlap valid: an entry that loses them
                // still ends with the node's (new) first `ov` bases.
                let head = qs[0][qs[0].len() - shift..].to_vec();
                for &((q, qf), _) in &preds {
                    let seq = &mut seqs[q];
                    if qf {
                        seq.truncate(seq.len() - shift);
                    } else {
                        seq.drain(..shift);
                    }
                }
                if vf {
                    let mut new = head;
                    new.extend_from_slice(&seqs[v]);
                    seqs[v] = new;
                } else {
                    seqs[v].extend_from_slice(&revcomp(&head));
                }
                moved += shift;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    moved
}

/// Rewrite overlapping links as 0M where one side has no other link.
///
/// The overlap of a link a -> b is sequence both a's end and b's start carry. When a's end has
/// no other link, those bases leave a (else, when b's start has no other link, they leave b):
/// every walk spells the same sequence and the break falls where the sequences diverge, which
/// after boundary normalization is the repeat's edge. Links whose two ends both branch keep
/// their overlap (trimming either side would break its other links).
fn blunt_links(seqs: &mut [Vec<u8>], links: &[ULink]) -> (Vec<ULink>, usize) {
    let mut degree: FxHashMap<(usize, bool), usize> = FxHashMap::default();
    // largest overlap of any link touching each segment: a trimmed segment must stay longer
    // than the overlaps its other links still need
    let mut max_ov = vec![0usize; seqs.len()];
    for &((a, b), ov) in links {
        *degree.entry(a).or_default() += 1;
        *degree.entry((b.0, !b.1)).or_default() += 1;
        max_ov[a.0] = max_ov[a.0].max(ov);
        max_ov[b.0] = max_ov[b.0].max(ov);
    }
    let mut out = Vec::with_capacity(links.len());
    let mut blunted = 0;
    for &((a, b), ov) in links {
        if ov == 0 || (a.0 == b.0 && a != b) {
            // a hairpin (a segment joined to its own reverse complement) keeps its overlap
            out.push(((a, b), ov));
            continue;
        }
        if degree[&a] == 1 && seqs[a.0].len() > ov + max_ov[a.0] {
            // leaving a at its right end when a is forward: its last `ov` bases go
            let seq = &mut seqs[a.0];
            if a.1 {
                seq.truncate(seq.len() - ov);
            } else {
                seq.drain(..ov);
            }
        } else if degree[&(b.0, !b.1)] == 1 && seqs[b.0].len() > ov + max_ov[b.0] {
            // entering b at its left end when b is forward: its first `ov` bases go
            let seq = &mut seqs[b.0];
            if b.1 {
                seq.drain(..ov);
            } else {
                seq.truncate(seq.len() - ov);
            }
        } else {
            out.push(((a, b), ov));
            continue;
        }
        blunted += 1;
        out.push(((a, b), 0));
    }
    (out, blunted)
}

/// Links between dead ends that overlap exactly by `min_ov..=max_ov` bases, when each is the
/// other's only longest match. A de Bruijn graph joins unitigs only through k - 1 bases; where
/// the k-mers across a junction are missing (short reads thin out over a hard stretch:
/// Nipponbare 153.9 kb, two unitigs sharing 84 bp at k 127) the two ends still share a shorter
/// stretch. The two ends of one unitig may join (a circle); one end to itself may not.
fn join_dead_ends(
    seqs: &[Vec<u8>],
    links: &[ULink],
    circular: &[bool],
    min_ov: usize,
    max_ov: usize,
) -> Vec<ULink> {
    let mut used: FxHashSet<(usize, bool)> = FxHashSet::default();
    for &((a, b), _) in links {
        used.insert(a); // right end of a.0 when a.1
        used.insert((b.0, !b.1)); // left end of b.0 when b.1
    }
    let oriented = |o: (usize, bool)| -> Vec<u8> {
        if o.1 {
            seqs[o.0].clone()
        } else {
            revcomp(&seqs[o.0])
        }
    };
    // leaving (u, f) uses u's right end when f; entering (v, f) uses v's left end when f
    let exits: Vec<(usize, bool)> = (0..seqs.len())
        .filter(|&u| !circular[u])
        .flat_map(|u| [(u, true), (u, false)])
        .filter(|&e| !used.contains(&e))
        .collect();
    let overlap = |x: &[u8], y: &[u8]| -> usize {
        let top = max_ov.min(x.len() - 1).min(y.len() - 1);
        (min_ov..=top)
            .rev()
            .find(|&l| x[x.len() - l..] == y[..l])
            .unwrap_or(0)
    };
    // best entry of every dead exit: (overlap, entry), None when tied
    type BestEntry = Option<(usize, (usize, bool))>;
    let mut best: FxHashMap<(usize, bool), BestEntry> = FxHashMap::default();
    for &x in &exits {
        let sx = oriented(x);
        let mut top: Option<(usize, (usize, bool))> = None;
        let mut tied = false;
        for &e in &exits {
            let y = (e.0, !e.1); // the entry through the same dead end
            if y.0 == x.0 && y.1 != x.1 {
                continue; // back into the end it leaves by
            }
            let l = overlap(&sx, &oriented(y));
            if l == 0 {
                continue;
            }
            match top {
                Some((m, _)) if l < m => {}
                Some((m, _)) if l == m => tied = true,
                _ => {
                    top = Some((l, y));
                    tied = false;
                }
            }
        }
        best.insert(x, if tied { None } else { top });
    }
    let mut out: BTreeSet<ULink> = BTreeSet::new();
    for (&x, &b) in &best {
        let Some((l, y)) = b else { continue };
        // mutual: the reverse walk from y's end comes back to x's end
        let back = best.get(&(y.0, !y.1)).copied().flatten();
        if back.map(|(_, z)| (z.0, !z.1)) != Some(x) {
            continue;
        }
        let rc = ((y.0, !y.1), (x.0, !x.1));
        let key = if (x, y) <= rc { (x, y) } else { rc };
        out.insert((key, l));
    }
    out.into_iter().collect()
}

/// Small k of the local assembly across a gap.
const GAP_K: usize = 31;
/// Bases at each dead end whose k-mers select the reads of the local assembly.
const GAP_WINDOW: usize = 150;
/// Longest stretch a gap fill may add (dense mode).
const GAP_MAX: usize = 600;
/// Longest stretch a gap fill may add between long-read dead ends: within one read's span.
const LONG_GAP_MAX: usize = 10_000;

fn canonical_small(kmer: &[u8]) -> Option<u64> {
    CanonicalKmers::new(kmer, kmer.len()).next().map(|(_, h)| h)
}

/// Dense mode: fill the gaps left between dead ends where no k-mer of the graph's k reached
/// the count threshold (short reads thin out over a hard stretch: Salvia at ~150x, k 77, lost
/// two ~50 bp stretches whose ends do not overlap). The reads touching a dead end's last
/// `GAP_WINDOW` bases are read again and assembled locally at k = 31: from each dead exit the
/// walk takes the one solid next base (or one at least 3x as common as the next best) until it
/// runs into the start of another dead end, whose first bases it must then spell. A gap found
/// from both of its ends is filled once, and an end claimed by two gaps by neither. The fill
/// extends the exit's unitig; the link to the entry overlaps by what the walk shares with it.
///
/// Long reads (sparse mode) leave such gaps too, at low depth: where every read differs from the
/// others somewhere in a stretch (HiFi homopolymer lengths), no k-mer of the graph's k = 1001
/// reaches the count threshold although the stretch is spanned. The Triglochin mitochondrion
/// (~20x) stayed a linear 219.9 kb unitig, its 810 bp closing stretch spanned by 21 reads whose
/// copies differ by a few bases (809-816 bp); Thuidium (~14x) likewise across 2,375 bp. There
/// the walk takes the most common next base (`dominance` 1: strictly more reads than the next
/// best), the reads' consensus, over up to `max_fill` bases.
fn fill_gaps(
    seqs: &mut [Vec<u8>],
    links: &mut Vec<ULink>,
    circular: &[bool],
    inputs: &[PathBuf],
    max_fill: usize,
    dominance: u32,
) -> Result<usize> {
    let mut used: FxHashSet<(usize, bool)> = FxHashSet::default();
    for &((a, b), _) in links.iter() {
        used.insert(a);
        used.insert((b.0, !b.1));
    }
    let exits: Vec<(usize, bool)> = (0..seqs.len())
        .filter(|&u| !circular[u] && seqs[u].len() >= GAP_K)
        .flat_map(|u| [(u, true), (u, false)])
        .filter(|&e| !used.contains(&e))
        .collect();
    if exits.len() < 2 {
        return Ok(0);
    }
    let oriented = |seqs: &[Vec<u8>], o: (usize, bool)| -> Vec<u8> {
        if o.1 {
            seqs[o.0].clone()
        } else {
            revcomp(&seqs[o.0])
        }
    };
    // tails of the dead exits select the reads; heads of the dead entries are the targets
    let mut anchors: FxHashSet<u64> = FxHashSet::default();
    type End = (usize, bool);
    let mut targets: FxHashMap<Vec<u8>, Option<(End, usize)>> = FxHashMap::default();
    for &x in &exits {
        let sx = oriented(seqs, x);
        let tail = &sx[sx.len().saturating_sub(GAP_WINDOW)..];
        anchors.extend(CanonicalKmers::new(tail, GAP_K).map(|(_, h)| h));
        let y = (x.0, !x.1); // the entry through the same dead end
        let sy = oriented(seqs, y);
        let head = &sy[..GAP_WINDOW.min(sy.len())];
        for p in 0..=head.len() - GAP_K {
            let e = targets
                .entry(head[p..p + GAP_K].to_vec())
                .or_insert(Some((y, p)));
            if *e != Some((y, p)) {
                *e = None; // the same 31-mer at two entries: not a target
            }
        }
    }
    // reads touching a dead end, and their k-mer counts
    let (rx, reader) = read_batches(inputs);
    let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
    for batch in rx {
        let hits: Vec<Vec<u64>> = batch
            .par_iter()
            .filter_map(|r| {
                let seq = r.seq.to_ascii_uppercase();
                let km: Vec<u64> = CanonicalKmers::new(&seq, GAP_K).map(|(_, h)| h).collect();
                km.iter().any(|h| anchors.contains(h)).then_some(km)
            })
            .collect();
        for km in hits {
            for h in km {
                *counts.entry(h).or_default() += 1;
            }
        }
    }
    reader.join().expect("reader thread panicked")?;
    let count = |kmer: &[u8]| {
        canonical_small(kmer)
            .and_then(|h| counts.get(&h).copied())
            .unwrap_or(0)
    };

    // (exit, entry, bases added, overlap with the entry)
    let mut fills: Vec<(End, End, Vec<u8>, usize)> = Vec::new();
    for &x in &exits {
        let sx = oriented(seqs, x);
        let tail = &sx[sx.len().saturating_sub(GAP_WINDOW)..];
        let mut tc: Vec<u32> = tail.windows(GAP_K).map(&count).collect();
        tc.sort_unstable();
        let solid = (tc.get(tc.len() / 2).copied().unwrap_or(0) / 5).max(2);
        let mut ext: Vec<u8> = sx[sx.len() - GAP_K..].to_vec();
        let mut added = Vec::new();
        let found = loop {
            if added.len() >= max_fill {
                break None;
            }
            let cur = &ext[ext.len() - (GAP_K - 1)..];
            let mut cands: Vec<(u32, u8)> = b"ACGT"
                .iter()
                .map(|&b| {
                    let mut km = cur.to_vec();
                    km.push(b);
                    (count(&km), b)
                })
                .filter(|&(c, _)| c >= solid)
                .collect();
            cands.sort_unstable_by(|a, b| b.cmp(a));
            let next = match cands.as_slice() {
                [] => break None,
                [only] => only.1,
                [first, second, ..] if first.0 >= dominance * second.0 && first.0 > second.0 => {
                    first.1
                }
                _ => break None,
            };
            ext.push(next);
            added.push(next);
            let last = &ext[ext.len() - GAP_K..];
            if let Some(&Some((y, p))) = targets.get(last) {
                if y == (x.0, !x.1) {
                    break None; // back into its own end
                }
                // the walk must spell the entry's first p + 31 bases
                let sy = oriented(seqs, y);
                let joined = [sx.as_slice(), added.as_slice()].concat();
                let ov = p + GAP_K;
                if joined.len() > ov && joined[joined.len() - ov..] == sy[..ov] {
                    break Some((y, ov));
                }
                break None;
            }
        };
        if let Some((y, ov)) = found {
            fills.push((x, y, added, ov));
        }
    }
    // one fill per gap (it is found from both ends) and per dead end
    let mut per_end: FxHashMap<(usize, bool), usize> = FxHashMap::default();
    for &(x, y, _, _) in &fills {
        *per_end.entry(x).or_default() += 1;
        *per_end.entry((y.0, !y.1)).or_default() += 1;
    }
    let mut done: FxHashSet<((usize, bool), (usize, bool))> = FxHashSet::default();
    let mut filled = 0;
    for (x, y, added, ov) in fills {
        let rc = ((y.0, !y.1), (x.0, !x.1));
        if done.contains(&rc) {
            continue;
        }
        // each end at most twice: once from its own walk, once from its partner's
        if per_end[&x] > 2 || per_end[&(y.0, !y.1)] > 2 {
            continue;
        }
        done.insert((x, y));
        if x.1 {
            seqs[x.0].extend_from_slice(&added);
        } else {
            let mut v = revcomp(&added);
            v.extend_from_slice(&seqs[x.0]);
            seqs[x.0] = v;
        }
        let key = if (x, y) <= rc { (x, y) } else { rc };
        links.push((key, ov));
        filled += 1;
    }
    Ok(filled)
}

#[derive(Default, Serialize)]
pub struct RescueStats {
    pub rounds: usize,
    /// Dead-end pairs reconnected.
    pub gaps: usize,
    /// Below-threshold k-mers taken back to do so.
    pub nodes: usize,
}

/// Graph over the kept k-mers: edges seen `min_count` times between kept nodes, plus `extra`.
/// Returns the node hashes (index = node id), adjacency and the number of edges kept.
fn build_graph(
    keep: &FxHashSet<u64>,
    edges: &FxHashMap<EdgeKey, u32>,
    extra: &[EdgeKey],
    min_count: u32,
) -> (Vec<u64>, AdjList, usize) {
    let mut ids: Vec<u64> = keep.iter().copied().collect();
    ids.sort_unstable();
    let index: FxHashMap<u64, usize> = ids.iter().enumerate().map(|(i, &h)| (h, i)).collect();
    let mut adj: Vec<Vec<(ONode, u32)>> = vec![Vec::new(); 2 * ids.len()];
    let slot = |o: ONode| 2 * o.0 + o.1 as usize;
    let mut kept: Vec<(ONode, ONode, u32)> = edges
        .iter()
        .filter(|(_, &c)| c >= min_count)
        .map(|(&e, _)| e)
        .chain(extra.iter().copied())
        .filter_map(|(a, b, d)| Some(((*index.get(&a.0)?, a.1), (*index.get(&b.0)?, b.1), d)))
        .collect();
    kept.sort_unstable();
    kept.dedup();
    let mut edges_kept = 0;
    for (a, b, d) in kept {
        adj[slot(a)].push((b, d));
        if (flip(b), flip(a)) != (a, b) {
            // not its own reverse complement (a hairpin), so the other strand is a second arc
            adj[slot(flip(b))].push((flip(a), d));
        }
        edges_kept += 1;
    }
    for v in &mut adj {
        v.sort_unstable();
        v.dedup();
    }
    (ids, adj, edges_kept)
}

/// Reconnect dead ends: for each unitig end with no link that some read leaves by, and a dead
/// entry that the read later reaches, take the run of sampled k-mers the read shows between
/// them. The most frequent run per dead exit, if at least two reads show it exactly, is
/// returned as (hash, canonical strand?, read position) steps from exit to entry. These are the
/// places where no k-mer reached the count threshold (e.g. a long homopolymer called
/// differently by most reads) but a group of reads still agrees.
type RunCounts = FxHashMap<Vec<Sample>, usize>;

fn rescue_gaps(
    ids: &[u64],
    unitigs: &[Unitig],
    links: &[ULink],
    read_samples: &[Vec<Sample>],
    k: usize,
) -> Vec<Vec<Sample>> {
    const MAX_SPAN: u32 = 30_000;
    let mut linked: FxHashSet<(usize, bool)> = FxHashSet::default();
    for &((a, b), _) in links {
        linked.insert(a); // leaves unitig a.0 at its right end when a.1
        linked.insert((b.0, !b.1)); // enters b.0 at its left end when b.1
    }
    // dead exits and entries as oriented k-mers (hash, strand) in walk direction
    let mut exits: FxHashSet<(u64, bool)> = FxHashSet::default();
    let mut entries: FxHashSet<(u64, bool)> = FxHashSet::default();
    for (u, t) in unitigs.iter().enumerate() {
        if t.circular {
            continue;
        }
        let (head, tail) = (t.nodes[0].0, t.nodes[t.nodes.len() - 1].0);
        if !linked.contains(&(u, true)) {
            exits.insert((ids[tail.0], tail.1));
            entries.insert((ids[tail.0], !tail.1));
        }
        if !linked.contains(&(u, false)) {
            exits.insert((ids[head.0], !head.1));
            entries.insert((ids[head.0], head.1));
        }
    }
    if exits.is_empty() {
        return Vec::new();
    }
    let live: FxHashSet<u64> = unitigs
        .iter()
        .flat_map(|t| t.nodes.iter().map(|&((i, _), _)| ids[i]))
        .collect();
    // dead exit -> k-mer run to a dead entry -> reads showing it
    let mut runs: FxHashMap<(u64, bool), RunCounts> = FxHashMap::default();
    for samples in read_samples {
        let len = samples.last().map_or(0, |s| s.2) + k as u32;
        let backward: Vec<Sample> = samples
            .iter()
            .rev()
            .map(|&(h, f, pos)| (h, !f, len - k as u32 - pos))
            .collect();
        for walk in [samples.as_slice(), backward.as_slice()] {
            for (i, &(h, f, pos)) in walk.iter().enumerate() {
                if !exits.contains(&(h, f)) {
                    continue;
                }
                for (j, &(h2, f2, pos2)) in walk.iter().enumerate().skip(i + 1) {
                    if pos2 - pos > MAX_SPAN {
                        break;
                    }
                    if entries.contains(&(h2, f2)) {
                        // the two ends may be one unitig's: closing it into a circle
                        let chain = &walk[i..=j];
                        if chain.windows(2).all(|w| ((w[1].2 - w[0].2) as usize) < k) {
                            // positions relative to the exit so identical runs compare equal
                            let rel: Vec<Sample> =
                                chain.iter().map(|&(a, b, q)| (a, b, q - pos)).collect();
                            *runs.entry((h, f)).or_default().entry(rel).or_default() += 1;
                        }
                        break;
                    }
                    if live.contains(&h2) {
                        // back in the graph somewhere other than a dead end: not a gap
                        break;
                    }
                }
            }
        }
    }
    let mut out = Vec::new();
    let mut used_entries: FxHashSet<(u64, bool)> = FxHashSet::default();
    let mut exit_keys: Vec<(u64, bool)> = runs.keys().copied().collect();
    exit_keys.sort_unstable();
    for key in exit_keys {
        let cands = &runs[&key];
        let (best, n) = cands
            .iter()
            .max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0)))
            .map(|(c, &n)| (c.clone(), n))
            .expect("non-empty");
        let last = best[best.len() - 1];
        // the reverse direction of the same gap is found from the other dead end: take it once
        let back = (last.0, !last.1);
        if n >= 2 && !used_entries.contains(&(last.0, last.1)) && !used_entries.contains(&back) {
            used_entries.insert((last.0, last.1));
            used_entries.insert((best[0].0, !best[0].1));
            out.push(best);
        }
    }
    out
}

#[derive(Default, Serialize)]
pub struct CleanStats {
    pub rounds: usize,
    pub tips_removed: usize,
    pub bubbles_popped: usize,
    pub isolated_removed: usize,
    /// Unitigs of components whose deepest unitig is below 30% of the genome's depth.
    pub shallow_component_removed: usize,
    /// Unitigs of components with no unitig reaching min(4k, 1 kb) and no circle (see
    /// `clean`).
    pub fragmented_component_removed: usize,
    /// Shallow arms spelling the same bases as the deeper way beside them (near copies: NUMTs,
    /// minor alleles), removed whatever their length (see `clean`).
    pub variant_arms_removed: usize,
    /// Bases of those arms beyond their shared overlap.
    pub variant_arm_bp: usize,
    /// Isolated circles shorter than k: tandem arrays whose length the graph cannot tell
    /// (e.g. a 62 bp nuclear satellite), removed (see `clean`).
    pub tandem_circles: Vec<TandemCircle>,
}

/// A removed tandem array: its repeat unit as the circle spells it.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct TandemCircle {
    pub period: usize,
    pub depth: f64,
    pub unit: String,
}

/// Word length and least share of an arm's words found in the deeper way beside it for the arm
/// to be a near copy of that way (see `clean`): one difference per 100 bp still leaves ~85%.
const VARIANT_W: usize = 15;
const VARIANT_SHARE: f64 = 0.8;

/// The first `want` bases (or all) of an oriented unitig, spelled from its k-mers.
fn unitig_prefix(t: &Unitig, kmers: &[&[u8]], k: usize, forward: bool, want: usize) -> Vec<u8> {
    let n = t.nodes.len();
    let mut seq: Vec<u8> = Vec::new();
    for i in 0..n {
        // walking a reversed unitig, node j follows node j + 1 at node j + 1's distance
        let (j, d) = if forward {
            (i, t.nodes[i].1)
        } else {
            (n - 1 - i, if i == 0 { 0 } else { t.nodes[n - i].1 })
        };
        let o = t.nodes[j].0;
        let km = if o.1 == forward {
            kmers[o.0].to_vec()
        } else {
            revcomp(kmers[o.0])
        };
        seq.extend_from_slice(if i == 0 { &km } else { &km[k - d as usize..] });
        if seq.len() >= want {
            break;
        }
    }
    seq
}

/// Share of `a`'s `VARIANT_W`-mers that occur in `b`.
fn shared_share(a: &[u8], b: &[u8]) -> f64 {
    fn words(s: &[u8]) -> impl Iterator<Item = u64> + '_ {
        s.windows(VARIANT_W).filter_map(|w| {
            w.iter().try_fold(0u64, |acc, &c| {
                let v = match c {
                    b'A' | b'a' => 0,
                    b'C' | b'c' => 1,
                    b'G' | b'g' => 2,
                    b'T' | b't' => 3,
                    _ => return None,
                };
                Some(acc << 2 | v)
            })
        })
    }
    let set: FxHashSet<u64> = words(b).collect();
    let (mut hit, mut all) = (0usize, 0usize);
    for w in words(a) {
        all += 1;
        hit += set.contains(&w) as usize;
    }
    if all == 0 {
        0.0
    } else {
        hit as f64 / all as f64
    }
}

/// One round of graph cleaning on unitigs; returns how many unitigs were removed.
///
/// Sequencing errors that recur in a minority of reads leave three shapes: short dead ends
/// (tips), short detours that rejoin the main path nearby (bubbles, simple or nested) and short
/// or shallow isolated pieces. Each is removed only when clearly shallower than its surroundings. Real
/// recombination joins distant parts of the genome and does not rejoin nearby, so it forms none
/// of these shapes. Near copies of the organelle sequence (NUMTs, minor alleles) leave arms of any
/// length that spell what the deeper way beside them spells; recombination leads into different
/// sequence, so the comparison tells them apart.
#[allow(clippy::too_many_arguments)]
fn clean(
    unitigs: &[Unitig],
    links: &[ULink],
    count: &[u32],
    k: usize,
    depth: f64,
    kmers: &[&[u8]],
    alive: &mut [bool],
    stats: &mut CleanStats,
) -> usize {
    let n = unitigs.len();
    let len: Vec<usize> = unitigs
        .iter()
        .map(|t| k + t.nodes[1..].iter().map(|(_, d)| *d as usize).sum::<usize>())
        .collect();
    let dp: Vec<f64> = unitigs
        .iter()
        .map(|t| t.nodes.iter().map(|(o, _)| count[o.0] as f64).sum::<f64>() / t.nodes.len() as f64)
        .collect();
    // neighbours at each end: [left, right]
    let mut nb: Vec<[Vec<usize>; 2]> = vec![[Vec::new(), Vec::new()]; n];
    for &(((a, af), (b, bf)), _) in links {
        nb[a][af as usize].push(b); // leaves a at its right end when a is forward
        nb[b][(!bf) as usize].push(a); // enters b at its left end when b is forward
    }
    for e in &mut nb {
        for side in e.iter_mut() {
            side.sort_unstable();
            side.dedup();
        }
    }
    // An error touches every k-mer over ~2k bp, and syncmer sampling can add up to k - s on each
    // side, so an error bubble or tip spans up to ~4k.
    // oriented unitig -> oriented unitigs it continues into
    let mut out: FxHashMap<(usize, bool), Vec<(usize, bool)>> = FxHashMap::default();
    for &(((a, af), (b, bf)), _) in links {
        out.entry((a, af)).or_default().push((b, bf));
        out.entry((b, !bf)).or_default().push((a, !af));
    }
    let short = 4 * k;
    // reference depth never below half the genome's depth
    let floor = 0.5 * depth;
    // The deepest of several ways on from a unitig end is never an error to remove: at a site
    // where systematic errors split the reads several ways (ONT homopolymers), the true branch
    // can hold barely half the reads, fall below half a repeat neighbour's depth and would
    // otherwise go in the same round as its error siblings, breaking the path. One such end is
    // enough (at its other end the true branch may compete with a repeat, which is deeper), but
    // an end with no rival does not count: an error branch is often the sole way into the short
    // piece where it rejoins.
    // Only one way is the deepest: at low depth the two lengths of a homopolymer that splits the
    // reads evenly often hold exactly equal counts (Tibet-0 at 150x: 11 and 11 against ~45 on
    // either side), and protecting both would leave the bubble standing. Ties go to more k-mer
    // observations, then to the lower index (the same way at both ends of a bubble).
    let mass: Vec<u64> = unitigs
        .iter()
        .map(|t| t.nodes.iter().map(|(o, _)| count[o.0] as u64).sum())
        .collect();
    let mut strongest_at = vec![false; n];
    for succ in out.values().filter(|s| s.len() > 1) {
        let top = succ.iter().map(|&(v, _)| v).max_by(|&x, &y| {
            dp[x]
                .total_cmp(&dp[y])
                .then(mass[x].cmp(&mass[y]))
                .then(y.cmp(&x))
        });
        if let Some(v) = top {
            strongest_at[v] = true;
        }
    }
    let mut remove = vec![false; n];
    for u in 0..n {
        let [l, r] = &nb[u];
        // Joined only to itself (a hairpin or self-loop, e.g. a (CG)n run that is its own
        // reverse complement): as isolated as a unitig with no links. A circle closed on itself
        // is compacted as `circular` and stays.
        let alone = l.iter().chain(r).all(|&v| v == u);
        match (l.is_empty() || alone, r.is_empty() || alone) {
            (true, true) => {
                let t = &unitigs[u];
                // A circle shorter than k is a tandem array with a unit shorter than k: no
                // molecule, and its length is not in the graph. Recruitment pulls in nuclear
                // satellites this way (Altai-5, Est-1: a 62 bp unit at half the organelle
                // depth, absent from the published mitogenome and plastome). The smallest
                // plant mitochondrial plasmids are over 1 kb, so real circles stay.
                let period = if t.circular {
                    t.nodes[1..].iter().map(|(_, d)| *d as usize).sum::<usize>() + t.close as usize
                } else {
                    usize::MAX
                };
                if period < k {
                    remove[u] = true;
                    let unit = unitig_prefix(t, kmers, k, true, period);
                    let circle = TandemCircle {
                        period,
                        depth: dp[u],
                        unit: String::from_utf8_lossy(&unit[..period.min(unit.len())]).into_owned(),
                    };
                    if !stats.tandem_circles.contains(&circle) {
                        stats.tandem_circles.push(circle);
                    }
                } else if !t.circular && (dp[u] < 0.25 * depth || len[u] < 2 * k) {
                    remove[u] = true;
                    stats.isolated_removed += 1;
                }
            }
            (true, false) | (false, true) => {
                let side = if l.is_empty() { r } else { l };
                let strongest = side.iter().map(|&v| dp[v]).fold(floor, f64::max);
                if len[u] < short && dp[u] < 0.5 * strongest && !strongest_at[u] {
                    remove[u] = true;
                    stats.tips_removed += 1;
                }
            }
            (false, false) => {
                let around = l.iter().chain(r).map(|&v| dp[v]).fold(floor, f64::max);
                // An arm shorter than 2k differs from its sibling in a few bases that reads
                // cannot place it by (ONT splits a homopolymer about half and half): unless it is
                // the deepest way on somewhere, it goes whatever its depth.
                if len[u] < short
                    && (dp[u] < 0.5 * around || len[u] < 2 * k)
                    && !strongest_at[u]
                    && has_detour(u, &out, &len, len[u] + 2 * k)
                {
                    remove[u] = true;
                    stats.bubbles_popped += 1;
                }
            }
        }
    }
    // A shallow arm that spells the bases of the deeper way beside it is a near copy of that
    // sequence: a NUMT recruited with the organelle reads (Col-0's chromosome 2 copy is ~99.9%
    // identical and, at ~0.1x the mitochondrial depth, sits at the min_count threshold) or a
    // minor allele. Its length is where the copy's differences happen to fall, far beyond an
    // error's 4k; a real alternative, such as a recombining repeat's other exit, leads into
    // different flanking sequence. It goes whatever its length, and is counted.
    let mut ov: FxHashMap<(OEnd, OEnd), usize> = FxHashMap::default();
    for &(((a, af), (b, bf)), o) in links {
        ov.insert(((a, af), (b, bf)), o);
        ov.insert(((b, !bf), (a, !af)), o);
    }
    let deepest_of = |succ: &[(usize, bool)], skip: usize| {
        succ.iter()
            .filter(|o| o.0 != skip)
            .copied()
            .max_by(|x, y| dp[x.0].total_cmp(&dp[y.0]))
    };
    let mut junctions: Vec<(&OEnd, &Vec<OEnd>)> = out.iter().filter(|(_, s)| s.len() > 1).collect();
    junctions.sort_unstable_by_key(|(from, _)| **from);
    for (&from, succ) in junctions {
        let Some(best) = deepest_of(succ, usize::MAX) else {
            continue;
        };
        for &arm in succ {
            let u = arm.0;
            if u == best.0
                || remove[u]
                || strongest_at[u]
                || unitigs[u].circular
                || dp[u] >= 0.5 * dp[best.0]
            {
                continue;
            }
            let a = unitig_prefix(&unitigs[u], kmers, k, arm.1, usize::MAX);
            let a = &a[ov[&(from, arm)].min(a.len())..];
            if a.len() < 2 * VARIANT_W {
                continue;
            }
            // the deeper way, followed through its deepest continuations for as long as the arm
            let want = a.len() + a.len() / 10 + 100;
            let mut b: Vec<u8> = Vec::new();
            let (mut prev, mut cur) = (from, best);
            for _ in 0..16 {
                let need = want - b.len() + ov[&(prev, cur)];
                let s = unitig_prefix(&unitigs[cur.0], kmers, k, cur.1, need);
                b.extend_from_slice(&s[ov[&(prev, cur)].min(s.len())..]);
                if b.len() >= want {
                    break;
                }
                match out.get(&cur).and_then(|s| deepest_of(s, u)) {
                    Some(nx) => (prev, cur) = (cur, nx),
                    None => break,
                }
            }
            b.truncate(want);
            if shared_share(a, &b) >= VARIANT_SHARE {
                remove[u] = true;
                stats.variant_arms_removed += 1;
                stats.variant_arm_bp += a.len();
            }
        }
    }
    let mut comp: Vec<usize> = (0..n).collect();
    fn root(c: &mut [usize], mut x: usize) -> usize {
        while c[x] != x {
            c[x] = c[c[x]];
            x = c[x];
        }
        x
    }
    for &(((a, _), (b, _)), _) in links {
        let (ra, rb) = (root(&mut comp, a), root(&mut comp, b));
        comp[ra] = rb;
    }
    let mut deepest: FxHashMap<usize, f64> = FxHashMap::default();
    // longest unitig of each component, or usize::MAX when it holds a circle
    let mut longest: FxHashMap<usize, usize> = FxHashMap::default();
    for (u, &d) in dp.iter().enumerate() {
        let r = root(&mut comp, u);
        let e = deepest.entry(r).or_insert(0.0);
        *e = e.max(d);
        let l = longest.entry(r).or_insert(0);
        *l = (*l).max(if unitigs[u].circular {
            usize::MAX
        } else {
            len[u]
        });
    }
    // A component with no unitig of some length is a tangle of a high-copy nuclear repeat or
    // low-complexity sequence: only unique sequence compacts into long unitigs, and an
    // organelle molecule, even a small plasmid, has some. (Nipponbare short reads, k 127: 25
    // such components, longest unitig <= 251 bp, none of it mitochondrial; every organelle
    // component of the benchmarks has one of >= 60 kb.)
    let fragmented = (4 * k).min(1000);
    for u in 0..n {
        let r = root(&mut comp, u);
        if remove[u] || unitigs[u].circular {
            continue;
        }
        if deepest[&r] < 0.3 * depth {
            remove[u] = true;
            stats.shallow_component_removed += 1;
        } else if longest[&r] < fragmented {
            remove[u] = true;
            stats.fragmented_component_removed += 1;
        }
    }
    let mut removed = 0;
    for t in unitigs
        .iter()
        .zip(&remove)
        .filter(|(_, &r)| r)
        .map(|(t, _)| t)
    {
        for (o, _) in &t.nodes {
            alive[o.0] = false;
        }
        removed += 1;
    }
    removed
}

/// k values tried by `run_auto_k`, largest first.
pub const K_LADDER: [usize; 8] = [1001, 701, 501, 351, 251, 201, 151, 101];
/// Node depth a k must reach to be taken by `run_auto_k`.
pub const AUTO_K_NODE_DEPTH: f64 = 22.0;

/// Assemble at the largest k of `K_LADDER` whose node depth reaches `AUTO_K_NODE_DEPTH`, or at
/// the smallest if none does.
///
/// A node is a k-mer seen error-free, and a HiFi read gives one over 1001 bp only about a third
/// of the time, so the node depth at k = 1001 is about a third of the read depth: enough at
/// Nipponbare's (node depth 226), not at 30x, where the graph falls apart (Zou cohort, 12
/// accessions thinned: one circle in 2/12 at k = 1001, 12/12 at k = 251). A shorter k sees
/// more error-free k-mers but resolves fewer repeats on its own (Nipponbare at 501: 10-14 of
/// 16 junctions, at 1001: 16), and too short merges short repeats (k = 101 at 30x: 10/12), so
/// k drops only as far as the depth demands. Across 15-30x the accessions came out as one
/// circle once the node depth reached about 22. Each k tried is a fresh pass over the reads,
/// cheap at the depths that need it.
///
/// When even k = 101 stays under the target (about 20x and less), the reads are self-corrected
/// (`ovasm correct` with its defaults) and the ladder is climbed again on the corrected reads,
/// written next to the graph (`<out>.corrected.fa`, named in the report). Correction makes a
/// read's k-mers error-free over longer stretches, which is what the node depth lacks there; at
/// higher depth it can only lose detail, so it is not used (12 Zou accessions thinned: one
/// circle within 10 bp in 9/12 at 15x and 12/12 at 20x with correction, 8 and 11 without; at
/// 30x 12/12 without, 11 with).
pub fn run_auto_k(
    inputs: Vec<PathBuf>,
    s: usize,
    min_count: Option<u32>,
    target: f64,
    out_gfa: &Path,
    out_json: &Path,
) -> Result<AssembleReport> {
    let ladder = |reads: &[PathBuf], tried: &mut Vec<(usize, f64)>| -> Result<AssembleReport> {
        let mut last = None;
        for (i, &k) in K_LADDER.iter().enumerate() {
            let r = run(
                reads.to_vec(),
                AssembleParams { k, s, min_count },
                out_gfa,
                out_json,
            )?;
            tried.push((k, r.node_depth));
            let done = r.node_depth >= target || i + 1 == K_LADDER.len();
            last = Some(r);
            if done {
                break;
            }
        }
        Ok(last.expect("the ladder is not empty"))
    };
    let mut tried = Vec::new();
    let mut r = ladder(&inputs, &mut tried)?;
    if r.node_depth < target {
        let corrected = out_gfa.with_extension("corrected.fa");
        crate::correct::run(
            inputs,
            Vec::new(),
            &[17, 25],
            &[15, 21, 25],
            crate::correct::CorrectParams {
                k: 17,
                min_count: None,
                length_tolerance: 0.2,
                max_gap: 1000,
                max_steps: 20_000,
                keep_ends: false,
                relative: 0.0,
                peak_fraction: 0.3,
            },
            &corrected,
            &out_gfa.with_extension("correct.json"),
        )?;
        r = ladder(std::slice::from_ref(&corrected), &mut tried)?;
        r.corrected_reads = Some(corrected.display().to_string());
    }
    r.k_tried = tried;
    fs::write(out_json, serde_json::to_string_pretty(&r)?)?;
    Ok(r)
}

pub fn run(
    inputs: Vec<PathBuf>,
    p: AssembleParams,
    out_gfa: &Path,
    out_json: &Path,
) -> Result<AssembleReport> {
    let dense = p.s == 0;
    if p.k % 2 == 0 || (!dense && (p.s >= p.k || p.s > crate::kmer::MAX_K || p.s < 5)) {
        bail!(
            "k must be odd and larger than s, and s in 5..={} (or 0: every k-mer)",
            crate::kmer::MAX_K
        );
    }
    let t0 = Instant::now();
    let (rx, reader) = read_batches(&inputs);
    // canonical hash -> (count, canonical sequence)
    let mut nodes: FxHashMap<u64, (u32, Vec<u8>)> = FxHashMap::default();
    // canonical ((a, fwd), (b, fwd), distance) -> count
    let mut edges: FxHashMap<EdgeKey, u32> = FxHashMap::default();
    let (mut reads, mut bases, mut sampled) = (0u64, 0u64, 0u64);
    // every read's sampled k-mers in read order, for gap rescue
    let mut read_samples: Vec<Vec<Sample>> = Vec::new();
    // Dense mode (short reads): a random error makes k new k-mers and edges seen once each,
    // tens of millions at organelle depth. Their first sighting only goes into a Bloom filter,
    // so counts start from the second. Gap rescue (for systematic long-read errors) is off:
    // keeping every read's k-mers would take gigabytes.
    let mut seen_nodes = dense.then(|| crate::correct::Bloom::new(30));
    let mut seen_edges = dense.then(|| crate::correct::Bloom::new(30));
    for mut batch in rx {
        for r in &mut batch {
            r.seq.make_ascii_uppercase();
        }
        let per_read: Vec<Vec<Sample>> =
            batch.par_iter().map(|r| sample_read(&r.seq, &p)).collect();
        for (r, samples) in batch.iter().zip(per_read) {
            reads += 1;
            bases += r.seq.len() as u64;
            sampled += samples.len() as u64;
            for &(h, fwd, pos) in &samples {
                if let Some(b) = seen_nodes.as_mut() {
                    if !b.insert(h) {
                        continue;
                    }
                }
                let e = nodes.entry(h).or_insert_with(|| {
                    let km = &r.seq[pos as usize..pos as usize + p.k];
                    (0, if fwd { km.to_vec() } else { revcomp(km) })
                });
                e.0 += 1;
            }
            for w in samples.windows(2) {
                let d = w[1].2 - w[0].2;
                if d as usize >= p.k {
                    continue; // a gap (non-ACGT run) between the two k-mers
                }
                let (a, b) = ((w[0].0, w[0].1), (w[1].0, w[1].1));
                let rc = ((b.0, !b.1), (a.0, !a.1));
                let key = if (a, b) <= rc { (a, b) } else { rc };
                if let Some(bl) = seen_edges.as_mut() {
                    let mut hs = FxHasher::default();
                    key.hash(&mut hs);
                    d.hash(&mut hs);
                    if !bl.insert(hs.finish()) {
                        continue;
                    }
                }
                *edges.entry((key.0, key.1, d)).or_default() += 1;
            }
            if !dense {
                read_samples.push(samples);
            }
        }
    }
    reader.join().expect("reader thread panicked")?;

    // Depth and threshold from the node-count distribution.
    let mut counts: Vec<u32> = nodes.values().map(|v| v.0).filter(|&c| c >= 3).collect();
    counts.sort_unstable();
    let node_depth = {
        let total: u64 = counts.iter().map(|&c| c as u64).sum();
        let mut acc = 0u64;
        counts
            .iter()
            .find(|&&c| {
                acc += c as u64;
                acc * 2 >= total
            })
            .copied()
            .unwrap_or(0) as f64
    };
    let min_count = p
        .min_count
        .unwrap_or_else(|| ((node_depth / 10.0) as u32).max(3));

    // Keep nodes and edges seen at least min_count times, plus whatever gap rescue adds.
    let mut keep: FxHashSet<u64> = nodes
        .iter()
        .filter(|(_, v)| v.0 >= min_count)
        .map(|(&h, _)| h)
        .collect();
    let mut extra: Vec<EdgeKey> = Vec::new();
    let mut cleaning = CleanStats::default();
    let mut rescue = RescueStats::default();
    let (ids, edges_kept, unitigs, links) = loop {
        let (ids, adj, edges_kept) = build_graph(&keep, &edges, &extra, min_count);
        let count: Vec<u32> = ids.iter().map(|h| nodes[h].0).collect();
        let kmers: Vec<&[u8]> = ids.iter().map(|h| nodes[h].1.as_slice()).collect();
        let mut alive = vec![true; ids.len()];
        let (unitigs, links) = loop {
            let unitigs = compact(&adj, &alive);
            let links = unitig_links(&adj, &alive, &unitigs, p.k);
            let removed = clean(
                &unitigs,
                &links,
                &count,
                p.k,
                node_depth,
                &kmers,
                &mut alive,
                &mut cleaning,
            );
            cleaning.rounds += 1;
            if removed == 0 || cleaning.rounds >= 50 {
                break (unitigs, links);
            }
        };
        if rescue.rounds >= 2 {
            break (ids, edges_kept, unitigs, links);
        }
        rescue.rounds += 1;
        let found = rescue_gaps(&ids, &unitigs, &links, &read_samples, p.k);
        if found.is_empty() {
            break (ids, edges_kept, unitigs, links);
        }
        for chain in &found {
            rescue.gaps += 1;
            for w in chain.windows(2) {
                let (a, b, d) = ((w[0].0, w[0].1), (w[1].0, w[1].1), w[1].2 - w[0].2);
                let rc = ((b.0, !b.1), (a.0, !a.1));
                let key = if (a, b) <= rc { (a, b) } else { rc };
                extra.push((key.0, key.1, d));
            }
            for &(h, _, _) in chain {
                if keep.insert(h) {
                    rescue.nodes += 1;
                }
            }
        }
    };

    let oseq = |o: ONode| -> Vec<u8> {
        let s = &nodes[&ids[o.0]].1;
        if o.1 {
            s.clone()
        } else {
            revcomp(s)
        }
    };
    let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(unitigs.len());
    for (u, t) in unitigs.iter().enumerate() {
        let mut seq = oseq(t.nodes[0].0);
        for &(n, d) in &t.nodes[1..] {
            let o = oseq(n);
            let ov = p.k - d as usize;
            if seq[seq.len() - ov..] != o[..ov] {
                bail!("unitig {u}: consecutive k-mers disagree on their overlap (hash collision?)");
            }
            seq.extend_from_slice(&o[ov..]);
        }
        seqs.push(seq);
    }
    // Dense mode keeps no per-read samples for gap rescue; dead ends facing each other are
    // joined by their exact sequence overlap instead.
    let mut links = links;
    let mut overlap_joins = 0;
    let circular: Vec<bool> = unitigs.iter().map(|t| t.circular).collect();
    let gaps_filled = if dense {
        let joins = join_dead_ends(&seqs, &links, &circular, (p.k / 3).max(31), p.k - 2);
        overlap_joins = joins.len();
        links.extend(joins);
        fill_gaps(&mut seqs, &mut links, &circular, &inputs, GAP_MAX, 3)?
    } else {
        fill_gaps(&mut seqs, &mut links, &circular, &inputs, LONG_GAP_MAX, 1)?
    };
    let boundary_bp_moved = normalize_boundaries(&mut seqs, &links);
    let (links, links_blunted) = blunt_links(&mut seqs, &links);
    let mut gfa = String::from("H\tVN:Z:1.0\n");
    let mut lengths = Vec::new();
    for (u, (t, seq)) in unitigs.iter().zip(&seqs).enumerate() {
        let mut counts: Vec<u32> = t.nodes.iter().map(|(n, _)| nodes[&ids[n.0]].0).collect();
        let depth: f64 = counts.iter().map(|&c| c as f64).sum::<f64>() / counts.len() as f64;
        // dm: the median node count. A stretch read far more often than the rest of a unitig
        // (a plastid-derived insert in the mitochondrion, recruited plastid reads piling onto
        // it) lifts the mean, not the median: Nipponbare's short-read unitig of 37 kb holding
        // 7 kb of plastid sequence has mean 1.72x the genome's depth and would pass as a
        // two-copy repeat, while a real repeat is deep along its whole length.
        counts.sort_unstable();
        let median = counts[counts.len() / 2];
        writeln!(
            gfa,
            "S\tu{u}\t{}\tLN:i:{}\tdp:f:{depth:.1}\tdm:f:{median}.0\tNC:i:{}",
            std::str::from_utf8(seq)?,
            seq.len(),
            t.nodes.len()
        )?;
        lengths.push(seq.len());
    }
    let sign = |f: bool| if f { '+' } else { '-' };
    for &(((a, af), (b, bf)), ov) in &links {
        writeln!(gfa, "L\tu{a}\t{}\tu{b}\t{}\t{ov}M", sign(af), sign(bf))?;
    }
    fs::create_dir_all(out_gfa.parent().unwrap_or(Path::new(".")))?;
    fs::write(out_gfa, &gfa)?;

    lengths.sort_unstable_by(|a, b| b.cmp(a));
    let total: usize = lengths.iter().sum();
    let mut acc = 0;
    let n50 = lengths
        .iter()
        .find(|&&l| {
            acc += l;
            acc * 2 >= total
        })
        .copied()
        .unwrap_or(0);
    let report = AssembleReport {
        k: p.k,
        s: p.s,
        reads,
        bases,
        sampled_kmers: sampled,
        distinct_nodes: nodes.len(),
        min_count,
        min_count_auto: p.min_count.is_none(),
        node_depth,
        nodes_kept: ids.len(),
        edges_kept,
        cleaning,
        rescue,
        overlap_joins,
        gaps_filled,
        boundary_bp_moved,
        links_blunted,
        links_with_overlap: links.iter().filter(|l| l.1 > 0).count(),
        unitigs: unitigs.len(),
        circular_unitigs: unitigs.iter().filter(|t| t.circular).count(),
        links: links.len(),
        total_length: total,
        n50,
        longest: lengths.first().copied().unwrap_or(0),
        elapsed_seconds: t0.elapsed().as_secs_f64(),
        peak_rss_mb: peak_rss_mb(),
        k_tried: Vec::new(),
        corrected_reads: None,
    };
    fs::create_dir_all(out_json.parent().unwrap_or(Path::new(".")))?;
    fs::write(out_json, serde_json::to_string_pretty(&report)?)
        .with_context(|| format!("cannot write {}", out_json.display()))?;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(n: usize, seed: u64) -> Vec<u8> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(s >> 62) as usize]
            })
            .collect()
    }

    #[test]
    fn syncmers_are_strand_symmetric_and_never_too_far_apart() {
        let (k, s) = (51, 11);
        let x = seq(5000, 3);
        let fw = closed_syncmers(&x, k, s);
        let rc = closed_syncmers(&revcomp(&x), k, s);
        let mirrored: Vec<usize> = rc.iter().rev().map(|&p| x.len() - k - p).collect();
        assert_eq!(fw, mirrored);
        assert!(
            fw.windows(2).all(|w| w[1] - w[0] <= k - s),
            "window guarantee"
        );
        let density = fw.len() as f64 / (x.len() - k + 1) as f64;
        assert!(
            (density - 2.0 / (k - s + 1) as f64).abs() < 0.02,
            "{density}"
        );
    }

    #[test]
    fn a_non_acgt_base_splits_runs() {
        let mut x = seq(400, 4);
        x[200] = b'N';
        let got = closed_syncmers(&x, 31, 9);
        assert!(got.iter().all(|&p| p + 31 <= 200 || p > 200));
    }

    fn fastq(dir: &Path, reads: &[Vec<u8>]) -> PathBuf {
        let p = dir.join("reads.fa");
        let text: String = reads
            .iter()
            .enumerate()
            .map(|(i, r)| format!(">r{i}\n{}\n", std::str::from_utf8(r).unwrap()))
            .collect();
        fs::write(&p, text).unwrap();
        p
    }

    fn tmpdir(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("ovasm-assemble-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    /// Tile reads over a circular genome, both strands, at the given depth.
    fn tile(genome: &[u8], len: usize, step: usize) -> Vec<Vec<u8>> {
        let circ = [genome, &genome[..len]].concat();
        (0..genome.len())
            .step_by(step)
            .enumerate()
            .map(|(i, st)| {
                let r = circ[st..st + len].to_vec();
                if i % 2 == 0 {
                    r
                } else {
                    revcomp(&r)
                }
            })
            .collect()
    }

    fn params() -> AssembleParams {
        AssembleParams {
            k: 101,
            s: 15,
            min_count: Some(3),
        }
    }

    #[test]
    fn a_long_read_gap_is_filled_by_the_reads_majority() {
        // Two dead ends 300 bp apart; five reads span the gap, three with A and two with C at
        // one gap position (a split call). The dense-mode rule (3x the next base) stops there;
        // the long-read rule takes the majority and reaches the entry.
        let (x, gap, y) = (seq(400, 21), seq(300, 22), seq(400, 23));
        let reads: Vec<Vec<u8>> = (0..5)
            .map(|i| {
                let mut g = gap.clone();
                g[150] = if i < 3 { b'A' } else { b'C' };
                [&x[200..], &g[..], &y[..200]].concat()
            })
            .collect();
        let dir = tmpdir("majority-fill");
        let input = fastq(&dir, &reads);
        let mut expected = gap.clone();
        expected[150] = b'A';
        for (dominance, filled) in [(3, false), (1, true)] {
            let mut seqs = vec![x.clone(), y.clone()];
            let mut links: Vec<ULink> = Vec::new();
            let n = fill_gaps(
                &mut seqs,
                &mut links,
                &[false, false],
                std::slice::from_ref(&input),
                LONG_GAP_MAX,
                dominance,
            )
            .unwrap();
            assert_eq!(n > 0, filled, "dominance {dominance}");
            if filled {
                // the exit unitig now spells the gap (and the link overlaps the entry's start)
                let joined = &seqs[0];
                assert!(joined
                    .windows(expected.len())
                    .any(|w| w == expected.as_slice()));
                assert_eq!(links.len(), 1);
            }
        }
        fs::remove_dir_all(&dir).ok();
    }

    fn gfa_segments(gfa: &str) -> Vec<String> {
        gfa.lines()
            .filter(|l| l.starts_with("S\t"))
            .map(|l| l.split('\t').nth(2).unwrap().to_string())
            .collect()
    }

    #[test]
    fn a_circular_genome_assembles_into_one_closed_unitig() {
        let g = seq(20_000, 11);
        let dir = tmpdir("circle");
        let reads = tile(&g, 3000, 100);
        let rep = run(
            vec![fastq(&dir, &reads)],
            params(),
            &dir.join("a.gfa"),
            &dir.join("a.json"),
        )
        .unwrap();
        assert_eq!(rep.unitigs, 1);
        let gfa = fs::read_to_string(dir.join("a.gfa")).unwrap();
        let u = gfa_segments(&gfa).remove(0);
        let links: Vec<&str> = gfa.lines().filter(|l| l.starts_with("L\t")).collect();
        assert_eq!(links.len(), 1, "closed by a self link");
        let ov: usize = links[0]
            .split('\t')
            .nth(5)
            .unwrap()
            .trim_end_matches('M')
            .parse()
            .unwrap();
        let circle = &u[..u.len() - ov];
        assert_eq!(circle.len(), g.len());
        let doubled = [g.clone(), g.clone()].concat();
        let rc = revcomp(circle.as_bytes());
        assert!(
            doubled
                .windows(g.len())
                .any(|w| w == circle.as_bytes() || w == rc.as_slice()),
            "the unitig is a rotation of the genome"
        );
    }

    #[test]
    fn recurring_errors_close_together_are_cleaned_back_to_one_circle() {
        // Two systematic errors 300 bp apart, each in a different fifth of the reads: above the
        // count threshold, their k-mers cross into each other and form a nested bubble.
        let g = seq(20_000, 31);
        let mut reads = tile(&g, 3000, 50);
        for (i, rd) in reads.iter_mut().enumerate() {
            let fwd = i % 2 == 0;
            let start = (i * 50) % g.len();
            for (err, every) in [(5000usize, 5usize), (5300, 7)] {
                let within = (err + g.len() - start) % g.len();
                if i % every == 0 && within < rd.len() {
                    let p = if fwd { within } else { rd.len() - 1 - within };
                    rd[p] = if rd[p] == b'A' { b'C' } else { b'A' };
                }
            }
        }
        let dir = tmpdir("nested");
        let rep = run(
            vec![fastq(&dir, &reads)],
            AssembleParams {
                k: 101,
                s: 15,
                min_count: None,
            },
            &dir.join("a.gfa"),
            &dir.join("a.json"),
        )
        .unwrap();
        assert!(
            rep.cleaning.bubbles_popped + rep.cleaning.tips_removed > 0,
            "errors formed structure"
        );
        assert_eq!(
            rep.unitigs,
            1,
            "{:?}",
            fs::read_to_string(dir.join("a.gfa"))
                .unwrap()
                .lines()
                .filter(|l| !l.starts_with("S"))
                .collect::<Vec<_>>()
        );
        assert_eq!(rep.circular_unitigs, 1);
    }

    #[test]
    fn a_detour_is_seen_where_its_arm_rejoins_a_step_later() {
        // A -> {T, E}; T -> Y; E -> F -> Y. E's exit F is itself part of the error arm, so the
        // arms rejoin only at Y, one step past E's exits.
        type Out = FxHashMap<(usize, bool), Vec<(usize, bool)>>;
        let link = |out: &mut Out, a: usize, b: usize| {
            out.entry((a, true)).or_default().push((b, true));
            out.entry((b, false)).or_default().push((a, false));
        };
        let (a, t, e, f, y) = (0, 1, 2, 3, 4);
        let len = vec![5000, 300, 300, 300, 5000];
        let mut out = Out::default();
        for (x, z) in [(a, t), (a, e), (t, y), (e, f), (f, y)] {
            link(&mut out, x, z);
        }
        assert!(has_detour(e, &out, &len, 300 + 2 * 101));
        // E -> F dead-ends: a branch, not a bubble
        let mut out = Out::default();
        for (x, z) in [(a, t), (a, e), (t, y), (e, f)] {
            link(&mut out, x, z);
        }
        assert!(!has_detour(e, &out, &len, 300 + 2 * 101));
    }

    #[test]
    fn dead_ends_sharing_a_shorter_overlap_are_joined_mutually() {
        // X ends with the first 60 bp of Y (k - 1 would be 100). Z also starts with 40 bp of
        // X's end: a shorter match that loses to Y.
        let x = seq(1000, 71);
        let y = [x[940..].to_vec(), seq(1000, 72)].concat();
        let z = [x[960..].to_vec(), seq(1000, 73)].concat();
        let seqs = vec![x, y, z];
        let joins = join_dead_ends(&seqs, &[], &[false; 3], 33, 99);
        assert_eq!(joins, vec![(((0, true), (1, true)), 60)]);
        // too short an overlap is not a join
        assert!(join_dead_ends(&seqs[..2], &[], &[false; 2], 61, 99).is_empty());
        // a linear unitig whose ends overlap closes into a circle
        let c = seq(1000, 74);
        let circle = [c.clone(), c[..50].to_vec()].concat();
        assert_eq!(
            join_dead_ends(&[circle], &[], &[false], 33, 99),
            vec![(((0, false), (0, false)), 50)] // the canonical side of the self-link
        );
    }

    #[test]
    fn a_short_run_joined_only_to_itself_is_removed_like_an_isolated_one() {
        // A circle plus reads of a (CG)n run: its k-mers are their own reverse complements'
        // neighbours, so the run's unitig links only to itself.
        let g = seq(20_000, 61);
        let mut reads = tile(&g, 3000, 100);
        let cg = b"CG".repeat(150);
        reads.extend((0..20).map(|_| cg.clone()));
        let dir = tmpdir("selfonly");
        let rep = run(
            vec![fastq(&dir, &reads)],
            params(),
            &dir.join("a.gfa"),
            &dir.join("a.json"),
        )
        .unwrap();
        assert_eq!((rep.unitigs, rep.circular_unitigs), (1, 1));
        let gfa = fs::read_to_string(dir.join("a.gfa")).unwrap();
        assert!(!gfa.contains("CGCGCGCGCGCGCGCGCGCGCGCG"), "the run is gone");
    }

    #[test]
    fn a_true_branch_holding_under_half_the_reads_survives_its_error_siblings() {
        // At 10,000 only 45% of the reads carry the true base; the rest split over the three
        // others (ONT-style systematic error). The true branch is below half its neighbours'
        // depth like every error branch, but it is the deepest way on, so it stays.
        let g = seq(20_000, 51);
        let mut reads = tile(&g, 3000, 50);
        for (i, rd) in reads.iter_mut().enumerate() {
            let start = (i * 50) % g.len();
            let within = (10_000 + g.len() - start) % g.len();
            if i % 20 >= 9 && within < rd.len() {
                let p = if i % 2 == 0 {
                    within
                } else {
                    rd.len() - 1 - within
                };
                let others: Vec<u8> = b"ACGT".iter().copied().filter(|&b| b != rd[p]).collect();
                rd[p] = others[i % 3];
            }
        }
        let dir = tmpdir("split");
        let rep = run(
            vec![fastq(&dir, &reads)],
            params(),
            &dir.join("a.gfa"),
            &dir.join("a.json"),
        )
        .unwrap();
        assert!(
            rep.cleaning.bubbles_popped > 0,
            "error branches formed bubbles"
        );
        assert_eq!((rep.unitigs, rep.circular_unitigs), (1, 1));
        let u = gfa_segments(&fs::read_to_string(dir.join("a.gfa")).unwrap()).remove(0);
        let (fw, rc) = (u.as_bytes().to_vec(), revcomp(u.as_bytes()));
        let site = &g[9_900..10_100];
        assert!(
            fw.windows(200).any(|w| w == site) || rc.windows(200).any(|w| w == site),
            "the true base was kept"
        );
    }

    #[test]
    fn a_gap_below_the_count_threshold_is_rescued_by_reads_that_agree() {
        // Every read over 10,000..10,028 carries its own error there except two clean ones, so
        // the true k-mers across it are seen twice (< min_count 3) and the circle breaks open.
        let g = seq(20_000, 41);
        let mut reads = tile(&g, 3000, 100);
        for (i, rd) in reads.iter_mut().enumerate() {
            let start = (i * 100) % g.len();
            let within = (10_000 + i % 28 + g.len() - start) % g.len();
            if start != 8000 && start != 8500 && within < rd.len() {
                let p = if i % 2 == 0 {
                    within
                } else {
                    rd.len() - 1 - within
                };
                rd[p] = if rd[p] == b'A' { b'C' } else { b'A' };
            }
        }
        let dir = tmpdir("rescue");
        let rep = run(
            vec![fastq(&dir, &reads)],
            params(),
            &dir.join("a.gfa"),
            &dir.join("a.json"),
        )
        .unwrap();
        assert_eq!(rep.rescue.gaps, 1);
        assert_eq!((rep.unitigs, rep.circular_unitigs), (1, 1));
        let gfa = fs::read_to_string(dir.join("a.gfa")).unwrap();
        let u = gfa_segments(&gfa).remove(0);
        let doubled = [g.clone(), g.clone()].concat();
        let rc = revcomp(u.as_bytes());
        assert!(
            doubled
                .windows(u.len().min(g.len()))
                .any(|w| u.as_bytes().starts_with(w) || rc.starts_with(w)),
            "the rescued circle spells the genome"
        );
    }

    #[test]
    fn a_repeat_collapses_into_a_branching_unitig_and_errors_drop_out() {
        // Circle: A R B R' ... with the same 600 bp repeat twice (k = 101 < 600).
        let (a, b, r) = (seq(8000, 21), seq(8000, 22), seq(600, 23));
        let g = [a.clone(), r.clone(), b.clone(), r.clone()].concat();
        let mut reads = tile(&g, 3000, 60);
        // Sprinkle one substitution into every 7th read: its k-mers become singletons.
        for (i, rd) in reads.iter_mut().enumerate().filter(|(i, _)| i % 7 == 0) {
            let p = (i * 37) % rd.len();
            rd[p] = if rd[p] == b'A' { b'C' } else { b'A' };
        }
        let dir = tmpdir("repeat");
        let rep = run(
            vec![fastq(&dir, &reads)],
            params(),
            &dir.join("a.gfa"),
            &dir.join("a.json"),
        )
        .unwrap();
        let gfa = fs::read_to_string(dir.join("a.gfa")).unwrap();
        let segs = gfa_segments(&gfa);
        assert_eq!(
            rep.unitigs, 3,
            "A-side, B-side and the collapsed repeat: {segs:?}"
        );
        assert_eq!(rep.links, 4);
        // The repeat unitig holds the sampled k-mers inside R: a substring of R (or its
        // reverse complement); the links' overlaps carry R's ends.
        let rstr = std::str::from_utf8(&r).unwrap();
        let rrc = String::from_utf8(revcomp(&r)).unwrap();
        // After boundary normalization the repeat unitig holds all of R (a few bases more only
        // where A and B happen to agree with each other next to R).
        let rep_seg: Vec<&String> = segs
            .iter()
            .filter(|s| s.contains(rstr) || s.contains(&rrc))
            .collect();
        assert_eq!(rep_seg.len(), 1, "{segs:?}");
        assert!(rep_seg[0].len() <= 600 + 8, "{}", rep_seg[0].len());
        assert!(rep.boundary_bp_moved > 0);
    }

    #[test]
    fn unitig_prefix_spells_both_orientations() {
        // ACGTTGCAAT as 5-mers at 0, 3 and 5; the middle one stored reverse complemented
        let kmers: Vec<&[u8]> = vec![b"ACGTT", b"TGCAA", b"GCAAT"];
        let t = Unitig {
            nodes: vec![((0, true), 0), ((1, false), 3), ((2, true), 2)],
            circular: false,
            close: 0,
        };
        assert_eq!(
            unitig_prefix(&t, &kmers, 5, true, usize::MAX),
            b"ACGTTGCAAT"
        );
        assert_eq!(
            unitig_prefix(&t, &kmers, 5, false, usize::MAX),
            b"ATTGCAACGT"
        );
        // stops once it has what was asked for
        assert_eq!(unitig_prefix(&t, &kmers, 5, true, 4), b"ACGTT");
        assert!(b"ATTGCAACGT".starts_with(&unitig_prefix(&t, &kmers, 5, false, 6)));
    }

    #[test]
    fn a_near_copy_shares_its_words_and_other_sequence_does_not() {
        let base: Vec<u8> = (0..2000u64)
            .map(|i| b"ACGT"[(mix64(i) % 4) as usize])
            .collect();
        let mut copy = base.clone();
        for i in (100..2000).step_by(400) {
            copy[i] = if copy[i] == b'A' { b'C' } else { b'A' }; // 5 differences in 2 kb
        }
        assert!(shared_share(&copy, &base) > 0.95);
        let other: Vec<u8> = (5000..7000u64)
            .map(|i| b"ACGT"[(mix64(i) % 4) as usize])
            .collect();
        assert!(shared_share(&other, &base) < 0.05);
        assert_eq!(shared_share(b"ACGT", &base), 0.0);
    }

    #[test]
    fn a_tandem_array_shorter_than_k_is_reported_not_kept_as_a_molecule() {
        // 40 copies of a 62 bp unit and nothing else: one isolated circle of period 62
        let unit = b"TCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGATACACACAGACTGTGAGTA";
        let array: Vec<u8> = unit.iter().copied().cycle().take(62 * 40).collect();
        let reads: Vec<Vec<u8>> = (0..30)
            .map(|i| array[(i * 37) % 600..(i * 37) % 600 + 1800].to_vec())
            .collect();
        let dir = tmpdir("tandem");
        let rep = run(
            vec![fastq(&dir, &reads)],
            AssembleParams {
                k: 101,
                s: 15,
                min_count: Some(3),
            },
            &dir.join("a.gfa"),
            &dir.join("a.json"),
        )
        .unwrap();
        assert_eq!(rep.unitigs, 0);
        assert_eq!(rep.cleaning.tandem_circles.len(), 1);
        let c = &rep.cleaning.tandem_circles[0];
        assert_eq!(c.period, 62);
        // the reported unit is a rotation of the true unit, on either strand
        let doubled = [unit.as_slice(), unit.as_slice()].concat();
        let rc: Vec<u8> = crate::evidence::revcomp(&doubled);
        let u = c.unit.as_bytes();
        assert!(
            doubled.windows(62).any(|w| w == u) || rc.windows(62).any(|w| w == u),
            "{}",
            c.unit
        );
    }

    #[test]
    fn automatic_k_keeps_1001_when_deep_and_steps_down_when_not() {
        let g = seq(12_000, 7);
        // deep: a read every 50 bp of 3 kb (60x) -> node depth well over the target at 1001
        let dir = tmpdir("autok-deep");
        let rep = run_auto_k(
            vec![fastq(&dir, &tile(&g, 3000, 50))],
            31,
            None,
            AUTO_K_NODE_DEPTH,
            &dir.join("a.gfa"),
            &dir.join("a.json"),
        )
        .unwrap();
        assert_eq!(rep.k, 1001);
        assert_eq!(rep.k_tried.len(), 1);
        // shallow: a read every 600 bp (5x) -> no k reaches the target; every k is tried,
        // largest first, and the smallest is kept
        let dir = tmpdir("autok-shallow");
        let rep = run_auto_k(
            vec![fastq(&dir, &tile(&g, 3000, 600))],
            31,
            None,
            AUTO_K_NODE_DEPTH,
            &dir.join("a.gfa"),
            &dir.join("a.json"),
        )
        .unwrap();
        // raw reads first, every k, then the same ladder on the self-corrected reads
        let ks: Vec<usize> = rep.k_tried.iter().map(|x| x.0).collect();
        assert_eq!(ks, [K_LADDER.to_vec(), K_LADDER.to_vec()].concat());
        assert_eq!(rep.k, *K_LADDER.last().unwrap());
        let corrected = rep
            .corrected_reads
            .clone()
            .expect("shallow reads are corrected");
        assert!(Path::new(&corrected).is_file());
        // the report on disk carries both ladders too
        let json: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(dir.join("a.json")).unwrap()).unwrap();
        assert_eq!(
            json["k_tried"].as_array().unwrap().len(),
            2 * K_LADDER.len()
        );
        assert_eq!(json["corrected_reads"].as_str(), Some(corrected.as_str()));
    }
}
