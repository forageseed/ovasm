//! Stage 2: read evidence on an assembly graph (`ovasm evidence`).
//!
//! Every GFA link becomes a junction sequence (end of the left segment + start of the right
//! segment, `flank` bp each side). A read supports a link when its exact k-mer hits on that
//! junction are collinear (one diagonal) and anchor at least `min_anchor` bp on both sides of
//! the joint. From the supported links of one read, consecutive joints spaced by a segment's
//! length form traversals `X -> S -> Y`; at a repeat with two entries and two exits these give
//! the two pairings (the "master" and the recombinant configuration) and the minority share.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::Serialize;

use crate::kmer::CanonicalKmers;
use crate::recruit::{peak_rss_mb, read_batches};

/// Oriented segment: (segment index, forward?).
pub type Oriented = (usize, bool);

pub struct Graph {
    pub names: Vec<String>,
    pub seqs: Vec<Vec<u8>>,
    /// (from, to) with GFA orientation.
    pub links: Vec<(Oriented, Oriented)>,
    /// Overlap (bp) of each link: the end of `from` equals the start of `to` over this length.
    pub overlaps: Vec<usize>,
}

pub fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

pub(crate) fn flip(o: Oriented) -> Oriented {
    (o.0, !o.1)
}

pub fn parse_gfa(path: &Path) -> Result<Graph> {
    let text =
        fs::read_to_string(path).with_context(|| format!("cannot read {}", path.display()))?;
    let mut g = Graph {
        names: Vec::new(),
        seqs: Vec::new(),
        links: Vec::new(),
        overlaps: Vec::new(),
    };
    let mut index: FxHashMap<String, usize> = FxHashMap::default();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.first() == Some(&"S") {
            if f.len() < 3 || f[2] == "*" {
                bail!(
                    "segment {:?} has no sequence; evidence needs sequences in the GFA",
                    f.get(1)
                );
            }
            index.insert(f[1].to_string(), g.names.len());
            g.names.push(f[1].to_string());
            g.seqs.push(f[2].as_bytes().to_ascii_uppercase());
        }
    }
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.first() != Some(&"L") {
            continue;
        }
        if f.len() < 6 {
            bail!("malformed link line: {line}");
        }
        let overlap = match f[5] {
            "*" => 0,
            cigar => cigar
                .strip_suffix('M')
                .and_then(|n| n.parse::<usize>().ok())
                .with_context(|| format!("link {line:?}: overlap must be <n>M"))?,
        };
        let seg = |name: &str| {
            index
                .get(name)
                .copied()
                .with_context(|| format!("link to unknown segment {name}"))
        };
        let (a, b) = ((seg(f[1])?, f[2] == "+"), (seg(f[3])?, f[4] == "+"));
        if overlap >= g.seqs[a.0].len() || overlap >= g.seqs[b.0].len() {
            bail!("link {line:?}: overlap is not shorter than both segments");
        }
        g.links.push((a, b));
        g.overlaps.push(overlap);
    }
    Ok(g)
}

impl Graph {
    pub(crate) fn oriented_seq(&self, o: Oriented) -> Vec<u8> {
        if o.1 {
            self.seqs[o.0].clone()
        } else {
            revcomp(&self.seqs[o.0])
        }
    }

    pub fn label(&self, o: Oriented) -> String {
        format!("{}{}", self.names[o.0], if o.1 { '+' } else { '-' })
    }

    /// Canonical key of a link: a link and its reverse complement are the same adjacency.
    fn link_key(a: Oriented, b: Oriented) -> (Oriented, Oriented) {
        let rc = (flip(b), flip(a));
        if (a, b) <= rc {
            (a, b)
        } else {
            rc
        }
    }
}

#[derive(Clone, Copy)]
pub struct EvidenceParams {
    pub k: usize,
    pub flank: usize,
    pub min_anchor: usize,
    /// Hits whose diagonals differ by at most this many bp are chained.
    pub diag_tol: i64,
    pub min_side_hits: usize,
    /// Allowed deviation of joint spacing from the middle segment length (fraction, floor 50 bp).
    pub spacing_tol: f64,
    pub bootstrap: usize,
    pub seed: u64,
    /// Longest segment still treated as a repeat inside phased paths.
    pub max_repeat_len: usize,
}

struct Junction {
    link: usize,
    /// Joint = end of the left segment in junction coordinates (the overlap lies just before it).
    left_len: usize,
    overlap: usize,
    /// Required anchor on each side: `min_anchor`, capped at 90% of a short segment
    /// (a read covering nearly all of a short segment is as strong as it gets).
    anchor_left: i64,
    anchor_right: i64,
}

struct Index {
    k: usize,
    junctions: Vec<Junction>,
    /// canonical k-mer -> (junction, position in junction, k-mer is forward-canonical there)
    map: FxHashMap<u64, Vec<(u32, u32, bool)>>,
}

fn forward_is_canonical(kmer_seq: &[u8], canon: u64) -> bool {
    let mut v = 0u64;
    for &b in kmer_seq {
        v = (v << 2)
            | match b.to_ascii_uppercase() {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                _ => 3,
            };
    }
    v == canon
}

fn build_index(g: &Graph, p: &EvidenceParams) -> Index {
    let mut idx = Index {
        k: p.k,
        junctions: Vec::new(),
        map: FxHashMap::default(),
    };
    for (li, &(a, b)) in g.links.iter().enumerate() {
        let sa = g.oriented_seq(a);
        let sb = g.oriented_seq(b);
        // keep room beyond the anchor so one error at the flank edge cannot veto support
        let flank = p.flank.max(p.min_anchor + 500);
        // Left keeps the overlap at its end; right starts after it. Support must reach the
        // unique sequence on both sides of the overlap.
        let ov = g.overlaps[li];
        let left = &sa[sa.len().saturating_sub(flank + ov)..];
        let right = &sb[ov..sb.len().min(ov + flank)];
        let mut j = left.to_vec();
        j.extend_from_slice(right);
        let ji = idx.junctions.len() as u32;
        for (pos, kmer) in CanonicalKmers::new(&j, p.k) {
            let fwd = forward_is_canonical(&j[pos..pos + p.k], kmer);
            idx.map.entry(kmer).or_default().push((ji, pos as u32, fwd));
        }
        idx.junctions.push(Junction {
            link: li,
            left_len: left.len(),
            overlap: ov,
            anchor_left: p.min_anchor.min((left.len() - ov) * 9 / 10) as i64,
            anchor_right: p.min_anchor.min(right.len() * 9 / 10) as i64,
        });
    }
    idx
}

/// A supported link in one read, in read coordinates and read direction: the read leaves
/// `from` at `from_end` and enters `to` at `to_start` (`to_start = from_end - overlap`).
#[derive(Clone, Copy, Debug)]
struct Crossing {
    link: usize,
    along: bool,
    to_start: i64,
    from_end: i64,
    /// Collinear k-mer hits supporting this crossing (both sides).
    score: u32,
}

fn crossings(read: &[u8], idx: &Index, p: &EvidenceParams) -> Vec<Crossing> {
    let k = idx.k as i64;
    // (junction, along) -> [(diagonal, junction position)]
    let mut hits: FxHashMap<(u32, bool), Vec<(i64, u32)>> = FxHashMap::default();
    for (rpos, kmer) in CanonicalKmers::new(read, idx.k) {
        let Some(entries) = idx.map.get(&kmer) else {
            continue;
        };
        let rfwd = forward_is_canonical(&read[rpos..rpos + idx.k], kmer);
        for &(ji, jpos, jfwd) in entries {
            let along = rfwd == jfwd;
            let diag = if along {
                rpos as i64 - jpos as i64
            } else {
                rpos as i64 + jpos as i64 + k
            };
            hits.entry((ji, along)).or_default().push((diag, jpos));
        }
    }
    let mut out = Vec::new();
    for ((ji, along), mut h) in hits {
        h.sort_unstable();
        let junction = &idx.junctions[ji as usize];
        let joint = junction.left_len as i64;
        let unique_left_end = joint - junction.overlap as i64;
        let (mut best_left, mut best_right, mut best_diag) = (0usize, 0usize, 0i64);
        let mut start = 0;
        while start < h.len() {
            let mut end = start + 1;
            while end < h.len() && h[end].0 - h[end - 1].0 <= p.diag_tol {
                end += 1;
            }
            let chain = &h[start..end];
            let min_j = chain.iter().map(|x| x.1 as i64).min().unwrap();
            let max_j = chain.iter().map(|x| x.1 as i64 + k).max().unwrap();
            let left = chain
                .iter()
                .filter(|x| (x.1 as i64 + k) <= unique_left_end)
                .count();
            let right = chain.iter().filter(|x| x.1 as i64 >= joint).count();
            if unique_left_end - min_j >= junction.anchor_left
                && max_j - joint >= junction.anchor_right
                && left >= p.min_side_hits
                && right >= p.min_side_hits
                && left.min(right) > best_left.min(best_right)
            {
                best_left = left;
                best_right = right;
                best_diag = chain[chain.len() / 2].0;
            }
            start = end;
        }
        if best_left > 0 {
            let ov = junction.overlap as i64;
            // Along: [.. left unique | overlap | right ..] -> `to` starts at the overlap.
            // Against: the read shows the reverse complement, so the overlap follows the joint.
            let (to_start, from_end) = if along {
                (best_diag + joint - ov, best_diag + joint)
            } else {
                (best_diag - joint, best_diag - joint + ov)
            };
            out.push(Crossing {
                link: junction.link,
                along,
                to_start,
                from_end,
                score: (best_left + best_right) as u32,
            });
        }
    }
    out
}

/// Decide between crossings that compete for the same place in a read: alternatives out of
/// (or into) the same segment end, e.g. `u7+ -> u9+` vs `u7+ -> u10+` when u9 and u10 begin
/// with the same sequence. `scores[i]` is the number of collinear k-mer hits of alternative i;
/// the shared side contributes equally, so differences come from the side that diverges.
/// Returns the index of the alternative the read supports, or `None` when the read cannot
/// tell them apart (it is then reported as ambiguous and counted for none of them).
fn pick_winner(scores: &[u32]) -> Option<usize> {
    // Best must beat the runner-up by both an absolute margin (noise on noisy reads) and a
    // relative one (the shared side inflates both scores on long reads).
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(scores[i]));
    let (best, second) = (scores[order[0]], scores.get(order[1]).copied().unwrap_or(0));
    (best >= second + 10.max(second / 20)).then_some(order[0])
}

/// Oriented (from, to) of a crossing in read direction.
fn crossing_ends(g: &Graph, c: &Crossing) -> (Oriented, Oriented) {
    let (a, b) = g.links[c.link];
    if c.along {
        (a, b)
    } else {
        (flip(b), flip(a))
    }
}

/// Whether the read shows `m`'s `to` segment in full: another crossing leaves that segment about
/// its length after `m` enters it (or, with `backward`, enters `m`'s `from` about its length
/// before `m` leaves it).
fn walked_through(g: &Graph, cs: &[Crossing], m: &Crossing, backward: bool, tol: f64) -> bool {
    let (mf, mt) = crossing_ends(g, m);
    let seg = if backward { mf } else { mt };
    let len = g.seqs[seg.0].len() as f64;
    cs.iter().any(|c| {
        let (cf, ct) = crossing_ends(g, c);
        let span = if backward {
            (ct == seg).then(|| (m.from_end - c.to_start) as f64)
        } else {
            (cf == seg).then(|| (c.from_end - m.to_start) as f64)
        };
        span.is_some_and(|s| (s - len).abs() <= (len * tol).max(50.0))
    })
}

/// Competitive assignment within one read. Returns the crossings kept and the links of groups
/// the read could not resolve.
///
/// Rivals are first compared by their collinear hits. When that is a tie, e.g. because one exit
/// is entirely a prefix of the other so the link's flanks cannot separate them, the read's next
/// (or previous) crossing decides: the rival whose segment the read demonstrably walks through
/// wins, if exactly one does.
fn assign(g: &Graph, cs: Vec<Crossing>, tol: f64) -> (Vec<Crossing>, Vec<usize>) {
    const SAME_PLACE: i64 = 100;
    let ends = |c: &Crossing| crossing_ends(g, c);
    let n = cs.len();
    let mut group: Vec<usize> = (0..n).collect();
    for i in 0..n {
        for j in i + 1..n {
            let ((fi, ti), (fj, tj)) = (ends(&cs[i]), ends(&cs[j]));
            let rival = cs[i].link != cs[j].link
                && ((fi == fj && (cs[i].from_end - cs[j].from_end).abs() <= SAME_PLACE)
                    || (ti == tj && (cs[i].to_start - cs[j].to_start).abs() <= SAME_PLACE));
            if rival {
                let (gi, gj) = (group[i], group[j]);
                group.iter_mut().filter(|x| **x == gj).for_each(|x| *x = gi);
            }
        }
    }
    let (mut kept, mut ambiguous) = (Vec::new(), Vec::new());
    let mut ids: Vec<usize> = group.clone();
    ids.sort_unstable();
    ids.dedup();
    for id in ids {
        let members: Vec<&Crossing> = (0..n).filter(|&i| group[i] == id).map(|i| &cs[i]).collect();
        if members.len() == 1 {
            kept.push(*members[0]);
            continue;
        }
        let scores: Vec<u32> = members.iter().map(|c| c.score).collect();
        let by_walk = || {
            // rivals sharing their `from` differ in where the read goes next, else where it came from
            let backward = ends(members[0]).0 != ends(members[1]).0;
            let through: Vec<usize> = (0..members.len())
                .filter(|&i| walked_through(g, &cs, members[i], backward, tol))
                .collect();
            (through.len() == 1).then(|| through[0])
        };
        match pick_winner(&scores).or_else(by_walk) {
            Some(w) => kept.push(*members[w]),
            None => ambiguous.extend(members.iter().map(|c| c.link)),
        }
    }
    (kept, ambiguous)
}

/// Directed step in read order: (from, to, read start of `to`, read end of `from`).
type Step = (Oriented, Oriented, i64, i64);

fn read_steps(g: &Graph, cs: &[Crossing]) -> Vec<Step> {
    let mut steps: Vec<Step> = cs
        .iter()
        .map(|c| {
            let (a, b) = g.links[c.link];
            if c.along {
                (a, b, c.to_start, c.from_end)
            } else {
                (flip(b), flip(a), c.to_start, c.from_end)
            }
        })
        .collect();
    steps.sort_by_key(|s| s.2);
    steps
}

/// Whether step `s2` continues step `s1` through their shared segment: the read enters the
/// segment at `s1` and leaves it at `s2`, spanning about its length.
fn continues(g: &Graph, s1: &Step, s2: &Step, tol: f64) -> bool {
    if s2.0 != s1.1 {
        return false;
    }
    let len = g.seqs[s1.1 .0].len() as f64;
    let span = (s2.3 - s1.2) as f64;
    (span - len).abs() <= (len * tol).max(50.0)
}

/// Walks in one read: maximal chains of steps, each continuing the previous one.
fn walks(g: &Graph, steps: &[Step], tol: f64) -> Vec<Vec<Oriented>> {
    let n = steps.len();
    let mut used = vec![false; n];
    let mut out = Vec::new();
    for i in 0..n {
        if used[i] {
            continue;
        }
        used[i] = true;
        let mut walk = vec![steps[i].0, steps[i].1];
        let mut cur = i;
        while let Some(j) =
            (cur + 1..n).find(|&j| !used[j] && continues(g, &steps[cur], &steps[j], tol))
        {
            used[j] = true;
            walk.push(steps[j].1);
            cur = j;
        }
        out.push(walk);
    }
    out
}

/// A walk and its reverse complement are the same path; keep the smaller.
fn canonical_walk(w: &[Oriented]) -> Vec<Oriented> {
    let rc: Vec<Oriented> = w.iter().rev().map(|&o| flip(o)).collect();
    if w <= rc.as_slice() {
        w.to_vec()
    } else {
        rc
    }
}

/// Traversals X -> S -> Y: three consecutive segments of a walk.
fn traversals(g: &Graph, steps: &[Step], tol: f64) -> Vec<(Oriented, Oriented, Oriented)> {
    let mut out: Vec<(Oriented, Oriented, Oriented)> = walks(g, steps, tol)
        .iter()
        .flat_map(|w| w.windows(3).map(canonical_walk).collect::<Vec<_>>())
        .map(|w| (w[0], w[1], w[2]))
        .collect();
    out.sort();
    out.dedup();
    out
}

/// A collapsed repeat is at least this many times deeper than the graph's median depth. On the
/// Nipponbare gold graph repeats sit at 1.56-2.35x and unique segments at most 1.19x.
const REPEAT_FOLD: f64 = 1.35;
/// A segment end's second-strongest way must hold at least this share of the strongest's reads
/// for the end to count as a crossroad end (see `crossroads`).
const CROSSROAD_SHARE: f64 = 1.0 / 3.0;

/// Mean depth of every segment from its GFA tags (dp, DP, KC/LN or RC/LN), if all have one.
fn segment_depths(text: &str, g: &Graph) -> Option<Vec<f64>> {
    let tags = crate::unify::segment_tags(text);
    g.names
        .iter()
        .zip(&g.seqs)
        .map(|(n, s)| {
            tags.get(n)
                .and_then(|t| crate::unify::depth(t, s.len()))
                .map(|d| d.0)
        })
        .collect()
}

/// Depths for calling repeats: ovasm assemble's median node depth (`dm`) where present, else
/// the backend's depth. A stretch read far more often than the rest of a segment (plastid
/// reads piling onto a plastid-derived insert) lifts the mean, not the median, while a real
/// repeat is deep along its whole length. Only repeat calls use it; the reported depths and
/// every other rule keep the mean, which near-threshold decisions were tuned on.
fn repeat_depths(text: &str, g: &Graph) -> Option<Vec<f64>> {
    let tags = crate::unify::segment_tags(text);
    g.names
        .iter()
        .zip(&g.seqs)
        .map(|(n, s)| {
            let t = tags.get(n)?;
            t.get("dm")
                .and_then(|v| v.parse::<f64>().ok())
                .or_else(|| crate::unify::depth(t, s.len()).map(|d| d.0))
        })
        .collect()
}

/// Length-weighted median depth: the depth of the base in the middle of the graph.
fn median_depth(g: &Graph, depth: &[f64]) -> f64 {
    let mut v: Vec<(f64, usize)> = depth
        .iter()
        .copied()
        .zip(g.seqs.iter().map(Vec::len))
        .collect();
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

/// Repeat segments: a branching end (two or more distinct read-supported links) and, when the
/// GFA records depth, at least `REPEAT_FOLD` times the median depth. Without depth the inference
/// is limited to segments short enough for reads to span. A depth-backed repeat remains a
/// repeat even when no read spans it; linearize reports its topological pairings as ambiguous.
/// Everything else anchors phased paths. Only supported links count: a spurious link from
/// the assembler must not turn a unique segment into a repeat. Depth matters because a graph's
/// branch points move with its k: a repeat close to k in length is not collapsed into a node of
/// its own but smeared onto the ends of unique segments, which then branch without being repeats.
/// Segments where two well-read ways meet at both ends. On one molecule a single-copy segment
/// has one way in and one way out; two ways at each end, both read often, mean two passes, so
/// the segment is a repeat whatever its depth says. At low depth a short repeat's depth is
/// noise (Zou accessions thinned to 15x, k = 101: 142-555 bp repeats at 0.84-1.37x the median,
/// each end with ways read 5-13 times), and as an anchor it can be passed once only, so the
/// molecule cannot close (Geg-14 at 30x: its 18.9 kb repeat at 1.37x, ends read 19/23 and
/// 18/22, came out linear and one copy short). A minority configuration leaves one way much
/// weaker than the other and stays out (`CROSSROAD_SHARE`).
fn crossroads(g: &Graph, link_reads: &[u64]) -> Vec<bool> {
    // (segment, right end?) -> reads of each distinct way (other segment end) leaving it
    let mut ends: FxHashMap<(usize, bool), FxHashMap<Oriented, u64>> = FxHashMap::default();
    for (&(a, b), &n) in g.links.iter().zip(link_reads) {
        if n == 0 || a.0 == b.0 {
            continue;
        }
        // leaves a's right end if a is forward; enters b's left end if b is forward
        for (end, other) in [((a.0, a.1), b), ((b.0, !b.1), flip(a))] {
            let e = ends.entry(end).or_default().entry(other).or_default();
            *e = (*e).max(n);
        }
    }
    let crossing = |end: (usize, bool)| {
        ends.get(&end).is_some_and(|ways| {
            let mut r: Vec<u64> = ways.values().copied().collect();
            r.sort_unstable_by(|x, y| y.cmp(x));
            r.len() >= 2 && r[1] as f64 >= CROSSROAD_SHARE * r[0] as f64
        })
    };
    (0..g.seqs.len())
        .map(|s| crossing((s, true)) && crossing((s, false)))
        .collect()
}

fn repeat_segments(
    g: &Graph,
    max_len: usize,
    supported: &[bool],
    depth: Option<&[f64]>,
) -> Vec<bool> {
    // (segment, right end?) -> distinct links touching that end
    let mut degree: FxHashMap<(usize, bool), usize> = FxHashMap::default();
    let mut keys: Vec<(Oriented, Oriented)> = g
        .links
        .iter()
        .zip(supported)
        .filter(|(_, &s)| s)
        .map(|(&(a, b), _)| Graph::link_key(a, b))
        .collect();
    keys.sort();
    keys.dedup();
    for (a, b) in keys {
        *degree.entry((a.0, a.1)).or_default() += 1; // leaves a's right end if a is forward
        *degree.entry((b.0, !b.1)).or_default() += 1; // enters b's left end if b is forward
    }
    let median = depth.map_or(0.0, |d| median_depth(g, d));
    let deep = |s: usize| match depth {
        Some(d) => d[s] >= REPEAT_FOLD * median,
        None => true,
    };
    let mut repeat: Vec<bool> = (0..g.seqs.len())
        .map(|s| {
            let branching = degree.get(&(s, true)).copied().unwrap_or(0) >= 2
                || degree.get(&(s, false)).copied().unwrap_or(0) >= 2;
            branching && (depth.is_some() || g.seqs[s].len() <= max_len) && deep(s)
        })
        .collect();
    // A deep segment inside a chain of repeats branches at neither end (its neighbours do):
    // one whose neighbours at both ends are all repeats is part of the chain (Nipponbare short
    // reads: a 4.4 kb two-copy segment between two repeats). Graph links count here, since
    // short reads cannot test every link; the depth keeps a deep unique segment out (an MTPT
    // at 3.75x between unique neighbours stays an anchor).
    if depth.is_some() {
        let mut side: FxHashMap<(usize, bool), Vec<usize>> = FxHashMap::default();
        for &(a, b) in &g.links {
            side.entry((a.0, a.1)).or_default().push(b.0);
            side.entry((b.0, !b.1)).or_default().push(a.0);
        }
        loop {
            let add: Vec<usize> = (0..g.seqs.len())
                .filter(|&s| !repeat[s] && deep(s))
                .filter(|&s| {
                    [true, false].iter().all(|&end| {
                        side.get(&(s, end))
                            .is_some_and(|v| v.iter().all(|&n| n != s && repeat[n]))
                    })
                })
                .collect();
            if add.is_empty() {
                break;
            }
            for s in add {
                repeat[s] = true;
            }
        }
    }
    repeat
}

/// Depth (fraction of the median) from which an ambiguously joined segment is a repeat.
const AMBIGUOUS_FOLD: f64 = 1.1;

/// Repeats the depth rule missed whose junction reads say so: a segment branching at both ends
/// (counting links that only ambiguous reads cross) with a link whose every crossing read fits
/// another link as well, i.e. whose junction occurs twice, at >= 1.1x the median depth. Salvia
/// short reads: a 628 bp inverted repeat (two copies) at 1.18x, two of its four links crossed
/// by 31 reads each and all of them ambiguous. A unique segment's junction reads are placed
/// (Nipponbare ONT: a 19.5 kb segment branching at both ends at 1.19x is not).
fn ambiguous_repeats(
    g: &Graph,
    repeat: &[bool],
    link_reads: &[u64],
    link_ambiguous: &[u64],
    depth: Option<&[f64]>,
    max_len: usize,
) -> Vec<usize> {
    let Some(d) = depth else {
        return Vec::new();
    };
    let median = median_depth(g, d);
    let mut ends: FxHashMap<(usize, bool), std::collections::BTreeSet<(Oriented, Oriented)>> =
        FxHashMap::default();
    let mut touches = vec![false; g.seqs.len()];
    for (li, &(a, b)) in g.links.iter().enumerate() {
        let (n, amb) = (link_reads[li], link_ambiguous[li]);
        if n == 0 && amb == 0 {
            continue;
        }
        let key = Graph::link_key(a, b);
        ends.entry((a.0, a.1)).or_default().insert(key);
        ends.entry((b.0, !b.1)).or_default().insert(key);
        if n == 0 {
            touches[a.0] = true;
            touches[b.0] = true;
        }
    }
    let branching = |s: usize, end: bool| ends.get(&(s, end)).map_or(0, |v| v.len()) >= 2;
    (0..g.seqs.len())
        .filter(|&s| {
            !repeat[s]
                && touches[s]
                && d[s] >= AMBIGUOUS_FOLD * median
                && g.seqs[s].len() <= max_len
                && branching(s, true)
                && branching(s, false)
        })
        .collect()
}

/// A bridge carrying at least this share of the reads leaving its anchor end is where that end
/// goes (not one of several configurations).
const HIDDEN_DOMINANT: f64 = 0.8;

/// Candidates for segments passed more than once that depth did not flag (it compresses on
/// noisy reads: a two-copy repeat measured 1.0-1.3x the median on ONT). Bridges are phased paths
/// and direct links between non-repeats. If an end of an unflagged segment is where two
/// different segments' ends each send nearly all their reads, both walks may cross it. A unique
/// segment next to a recombining repeat instead splits its reads between configurations, so no
/// single bridge dominates. This also fires on a plain branch whose two sides dead-end into it,
/// hence only candidates: `crossed_in_two_contexts` decides.
fn hidden_repeats(
    g: &Graph,
    repeat: &[bool],
    phased: &BTreeMap<Vec<Oriented>, u64>,
    link_reads: &[u64],
    max_len: usize,
) -> Vec<usize> {
    let mut br: Vec<(Oriented, Oriented, u64)> = Vec::new();
    for (w, &n) in phased {
        let (a, b) = (w[0], w[w.len() - 1]);
        br.push((a, b, n));
        br.push((flip(b), flip(a), n));
    }
    let mut seen = std::collections::BTreeSet::new();
    for (li, &(a, b)) in g.links.iter().enumerate() {
        let n = link_reads[li];
        if repeat[a.0] || repeat[b.0] || n == 0 || !seen.insert(Graph::link_key(a, b)) {
            continue;
        }
        br.push((a, b, n));
        if (flip(b), flip(a)) != (a, b) {
            br.push((flip(b), flip(a), n));
        }
    }
    let mut total: BTreeMap<Oriented, u64> = BTreeMap::new();
    for &(a, _, n) in &br {
        *total.entry(a).or_default() += n;
    }
    // entered end -> segments whose end sends nearly all its reads there
    let mut into: BTreeMap<Oriented, std::collections::BTreeSet<usize>> = BTreeMap::new();
    for &(a, b, n) in &br {
        if a.0 != b.0 && n >= 3 && n as f64 >= HIDDEN_DOMINANT * total[&a] as f64 {
            into.entry(b).or_default().insert(a.0);
        }
    }
    let mut out: Vec<usize> = into
        .iter()
        .filter(|(b, from)| from.len() >= 2 && !repeat[b.0] && g.seqs[b.0].len() <= max_len)
        .map(|(b, _)| b.0)
        .collect();
    out.sort_unstable();
    out.dedup();
    out
}

/// Depth (fraction of the median) from which a segment is part of the main genome.
const MAJOR_DEPTH: f64 = 0.5;

/// Whether reads cross `s` (inside a phased path, counted as a repeat) in two contexts with
/// different segments on both sides, each seen in at least 3 reads and flanked by `major`
/// segments: two copies of it in the genome. A unique segment is crossed in one context, or
/// also in a low-frequency configuration whose flanks are shallow (a conformation that reuses
/// it, not a second copy: at Nipponbare 55 kb a 1.3 kb segment is crossed by 157 reads between
/// flanks at 0.8x and by 115 between flanks at 0.3-0.4x).
fn crossed_in_two_contexts(
    phased: &BTreeMap<Vec<Oriented>, u64>,
    s: usize,
    major: &[bool],
) -> bool {
    let mut contexts: Vec<(Oriented, Oriented)> = Vec::new();
    for (w, &n) in phased {
        if n < 3 || !major[w[0].0] || !major[w[w.len() - 1].0] {
            continue;
        }
        for i in 1..w.len().saturating_sub(1) {
            if w[i].0 != s {
                continue;
            }
            // oriented so that s reads forward
            let (a, b) = if w[i].1 {
                (w[0], w[w.len() - 1])
            } else {
                (flip(w[w.len() - 1]), flip(w[0]))
            };
            contexts.push((a, b));
        }
    }
    contexts.iter().enumerate().any(|(i, x)| {
        contexts[i + 1..]
            .iter()
            .any(|y| x.0 .0 != y.0 .0 && x.1 .0 != y.1 .0)
    })
}

/// Locally phased paths in a walk: unique segment, one or more repeats, unique segment.
fn phased_paths(walk: &[Oriented], repeat: &[bool]) -> Vec<Vec<Oriented>> {
    let anchors: Vec<usize> = (0..walk.len()).filter(|&i| !repeat[walk[i].0]).collect();
    anchors
        .windows(2)
        .filter(|w| w[1] - w[0] >= 2)
        .map(|w| canonical_walk(&walk[w[0]..=w[1]]))
        .collect()
}

#[derive(Serialize)]
pub struct LinkSupport {
    pub from: String,
    pub to: String,
    /// Reads assigned to this link (competing alternatives resolved per read).
    pub reads: u64,
    /// Reads that matched this link and a rival equally well (counted for neither).
    pub ambiguous_reads: u64,
}

#[derive(Serialize)]
pub struct Alternative {
    pub to: String,
    pub reads: u64,
    pub ambiguous_reads: u64,
}

#[derive(Serialize)]
pub struct BranchSupport {
    /// Oriented segment whose outgoing end branches (reverse-complemented ends are merged).
    pub from: String,
    pub alternatives: Vec<Alternative>,
    pub minority_fraction: Option<f64>,
    pub minority_ci95: Option<(f64, f64)>,
}

#[derive(Serialize)]
pub struct Traversal {
    pub path: String,
    pub reads: u64,
}

#[derive(Serialize)]
pub struct RepeatPairing {
    /// Middle segment (repeat) in the orientation the pairing is reported.
    pub repeat: String,
    pub entries: Vec<String>,
    pub exits: Vec<String>,
    /// Reads for (entry0->exit0 + entry1->exit1) and (entry0->exit1 + entry1->exit0).
    pub pairing_reads: (u64, u64),
    pub minority_fraction: Option<f64>,
    pub minority_ci95: Option<(f64, f64)>,
}

#[derive(Serialize)]
pub struct EvidenceReport {
    pub graph: String,
    pub reads_scanned: u64,
    pub reads_with_support: u64,
    pub links: Vec<LinkSupport>,
    pub links_supported: usize,
    pub branches: Vec<BranchSupport>,
    pub traversals: Vec<Traversal>,
    pub repeats: Vec<RepeatPairing>,
    /// Segments treated as repeats. Depth-backed repeats have no length cap; without depth,
    /// topology-only inference is limited to `max_repeat_len`.
    pub repeat_segments: Vec<String>,
    /// Of those, the ones found from phasing rather than depth (see `hidden_repeats`).
    pub hidden_repeats: Vec<String>,
    /// Of those, the ones found from junction reads that fit two copies (see
    /// `ambiguous_repeats`).
    pub ambiguous_repeats: Vec<String>,
    /// Per-segment depth from the GFA tags, when every segment has one.
    pub segment_depths: Option<BTreeMap<String, f64>>,
    /// Length-weighted median of `segment_depths`.
    pub median_depth: Option<f64>,
    /// Unique -> repeat(s) -> unique paths observed in single reads, most supported first.
    pub phased_paths: Vec<Traversal>,
    pub params: serde_json::Value,
    pub elapsed_seconds: f64,
    pub peak_rss_mb: Option<f64>,
}

fn splitmix(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Point estimate and bootstrap 95% interval of the share of `minority` among `minority + majority`.
pub fn minority_share(
    minority: u64,
    majority: u64,
    reps: usize,
    seed: u64,
) -> (Option<f64>, Option<(f64, f64)>) {
    let n = minority + majority;
    if n == 0 {
        return (None, None);
    }
    let est = minority as f64 / n as f64;
    if reps == 0 {
        return (Some(est), None);
    }
    let mut state = seed ^ n;
    let mut draws: Vec<f64> = (0..reps)
        .map(|_| {
            let mut m = 0u64;
            for _ in 0..n {
                if splitmix(&mut state) % n < minority {
                    m += 1;
                }
            }
            m as f64 / n as f64
        })
        .collect();
    draws.sort_by(f64::total_cmp);
    let lo = draws[(reps as f64 * 0.025) as usize];
    let hi = draws[((reps as f64 * 0.975) as usize).min(reps - 1)];
    (Some(est), Some((lo, hi)))
}

pub fn run(
    gfa: &Path,
    inputs: Vec<PathBuf>,
    p: EvidenceParams,
    out: &Path,
) -> Result<EvidenceReport> {
    let t0 = Instant::now();
    let g = parse_gfa(gfa)?;
    let gfa_text = fs::read_to_string(gfa)?;
    let depths = segment_depths(&gfa_text, &g);
    let rdepths = repeat_depths(&gfa_text, &g);
    let idx = build_index(&g, &p);
    let (rx, reader) = read_batches(&inputs);

    let mut link_reads = vec![0u64; g.links.len()];
    let mut link_ambiguous = vec![0u64; g.links.len()];
    let mut traversal_reads: BTreeMap<(Oriented, Oriented, Oriented), u64> = BTreeMap::new();
    let (mut scanned, mut supported) = (0u64, 0u64);
    type PerRead = (
        Vec<Crossing>,
        Vec<usize>,
        Vec<(Oriented, Oriented, Oriented)>,
        Vec<Vec<Oriented>>,
    );
    // Walks of each informative read; phased paths need the repeats, and the repeats need the
    // link support of the whole scan.
    let mut read_walks: Vec<Vec<Vec<Oriented>>> = Vec::new();
    for batch in rx {
        scanned += batch.len() as u64;
        let per_read: Vec<PerRead> = batch
            .par_iter()
            .map(|r| {
                let (cs, ambiguous) = assign(&g, crossings(&r.seq, &idx, &p), p.spacing_tol);
                let steps = read_steps(&g, &cs);
                let tr = traversals(&g, &steps, p.spacing_tol);
                let ws: Vec<Vec<Oriented>> = walks(&g, &steps, p.spacing_tol)
                    .into_iter()
                    .filter(|w| w.len() >= 3)
                    .collect();
                (cs, ambiguous, tr, ws)
            })
            .collect();
        for (cs, ambiguous, tr, ws) in per_read {
            if !ws.is_empty() {
                read_walks.push(ws);
            }
            if !cs.is_empty() {
                supported += 1;
            }
            for (links, counts) in [
                (
                    cs.iter().map(|c| c.link).collect::<Vec<_>>(),
                    &mut link_reads,
                ),
                (ambiguous, &mut link_ambiguous),
            ] {
                let mut seen = links;
                seen.sort_unstable();
                seen.dedup();
                for l in seen {
                    counts[l] += 1;
                }
            }
            for t in tr {
                *traversal_reads.entry(t).or_default() += 1;
            }
        }
    }
    reader.join().expect("reader thread panicked")?;

    let supported_links: Vec<bool> = link_reads.iter().map(|&n| n > 0).collect();
    let mut repeat = repeat_segments(&g, p.max_repeat_len, &supported_links, rdepths.as_deref());
    for (r, c) in repeat.iter_mut().zip(crossroads(&g, &link_reads)) {
        *r |= c;
    }
    let ambiguous = ambiguous_repeats(
        &g,
        &repeat,
        &link_reads,
        &link_ambiguous,
        depths.as_deref(),
        p.max_repeat_len,
    );
    for &s in &ambiguous {
        repeat[s] = true;
    }
    let phase = |repeat: &[bool]| {
        let mut out: BTreeMap<Vec<Oriented>, u64> = BTreeMap::new();
        for ws in &read_walks {
            let mut ph: Vec<Vec<Oriented>> =
                ws.iter().flat_map(|w| phased_paths(w, repeat)).collect();
            ph.sort();
            ph.dedup();
            for path in ph {
                *out.entry(path).or_default() += 1;
            }
        }
        out
    };
    let mut phased_reads = phase(&repeat);
    // Segments at genome-level depth: a context between two of them is the main genome, not a
    // low-frequency configuration (the same cut as linearize's optional anchors).
    let major: Vec<bool> = match depths.as_deref() {
        Some(d) => {
            let md = median_depth(&g, d);
            d.iter().map(|&x| x >= MAJOR_DEPTH * md).collect()
        }
        None => vec![true; g.seqs.len()],
    };
    // Repeats the depth missed, shown by the reads: candidates are tried as repeats and kept
    // when reads cross them in two different contexts; rephase until none is left.
    let mut hidden: Vec<usize> = Vec::new();
    for _ in 0..4 {
        let cands = hidden_repeats(&g, &repeat, &phased_reads, &link_reads, p.max_repeat_len);
        if cands.is_empty() {
            break;
        }
        let mut trial = repeat.clone();
        for &s in &cands {
            trial[s] = true;
        }
        let tried = phase(&trial);
        let confirmed: Vec<usize> = cands
            .into_iter()
            .filter(|&s| crossed_in_two_contexts(&tried, s, &major))
            .collect();
        if confirmed.is_empty() {
            break;
        }
        for s in confirmed {
            repeat[s] = true;
            hidden.push(s);
        }
        phased_reads = phase(&repeat);
    }

    // Branches: group links by their canonical outgoing oriented end.
    let mut by_end: BTreeMap<Oriented, Vec<(Oriented, u64, u64)>> = BTreeMap::new();
    let mut seen_keys = std::collections::BTreeSet::new();
    for (li, &(a, b)) in g.links.iter().enumerate() {
        if !seen_keys.insert(Graph::link_key(a, b)) {
            continue;
        }
        let (n, amb) = (link_reads[li], link_ambiguous[li]);
        by_end.entry(a).or_default().push((b, n, amb));
        if (flip(b), flip(a)) != (a, b) {
            by_end.entry(flip(b)).or_default().push((flip(a), n, amb));
        }
    }
    let mut branches = Vec::new();
    let mut seed = p.seed;
    for (from, mut alts) in by_end {
        if alts.len() < 2 {
            continue;
        }
        alts.sort_by_key(|x| std::cmp::Reverse(x.1));
        let total: u64 = alts.iter().map(|a| a.1).sum();
        let (frac, ci) = minority_share(total - alts[0].1, alts[0].1, p.bootstrap, seed);
        seed = seed.wrapping_add(1);
        branches.push(BranchSupport {
            from: g.label(from),
            alternatives: alts
                .iter()
                .map(|(to, n, amb)| Alternative {
                    to: g.label(*to),
                    reads: *n,
                    ambiguous_reads: *amb,
                })
                .collect(),
            minority_fraction: frac,
            minority_ci95: ci,
        });
    }

    // Repeat pairings: a middle segment with exactly two entries and two exits among traversals
    // (in one orientation; the reverse orientation is the same repeat).
    let mut repeats = Vec::new();
    let mut by_mid: BTreeMap<Oriented, Vec<(Oriented, Oriented, u64)>> = BTreeMap::new();
    for (&(x, s, y), &n) in &traversal_reads {
        let key = if s.1 {
            (x, s, y)
        } else {
            (flip(y), flip(s), flip(x))
        };
        by_mid.entry(key.1).or_default().push((key.0, key.2, n));
    }
    for (mid, trs) in &by_mid {
        let mut entries: Vec<Oriented> = trs.iter().map(|t| t.0).collect();
        let mut exits: Vec<Oriented> = trs.iter().map(|t| t.1).collect();
        entries.sort();
        entries.dedup();
        exits.sort();
        exits.dedup();
        if entries.len() != 2 || exits.len() != 2 {
            continue;
        }
        let count = |x: Oriented, y: Oriented| {
            trs.iter()
                .filter(|t| t.0 == x && t.1 == y)
                .map(|t| t.2)
                .sum::<u64>()
        };
        let a = count(entries[0], exits[0]) + count(entries[1], exits[1]);
        let b = count(entries[0], exits[1]) + count(entries[1], exits[0]);
        let (frac, ci) = minority_share(a.min(b), a.max(b), p.bootstrap, seed);
        seed = seed.wrapping_add(1);
        repeats.push(RepeatPairing {
            repeat: g.label(*mid),
            entries: entries.iter().map(|&o| g.label(o)).collect(),
            exits: exits.iter().map(|&o| g.label(o)).collect(),
            pairing_reads: (a, b),
            minority_fraction: frac,
            minority_ci95: ci,
        });
    }

    let links: Vec<LinkSupport> = g
        .links
        .iter()
        .zip(link_reads.iter().zip(&link_ambiguous))
        .map(|(&(a, b), (&n, &amb))| LinkSupport {
            from: g.label(a),
            to: g.label(b),
            reads: n,
            ambiguous_reads: amb,
        })
        .collect();
    let report = EvidenceReport {
        graph: gfa.display().to_string(),
        reads_scanned: scanned,
        reads_with_support: supported,
        links_supported: link_reads.iter().filter(|&&n| n > 0).count(),
        links,
        branches,
        traversals: traversal_reads
            .iter()
            .map(|(&(x, s, y), &n)| Traversal {
                path: format!("{} {} {}", g.label(x), g.label(s), g.label(y)),
                reads: n,
            })
            .collect(),
        repeats,
        segment_depths: depths
            .as_ref()
            .map(|d| g.names.iter().cloned().zip(d.iter().copied()).collect()),
        median_depth: depths.as_deref().map(|d| median_depth(&g, d)),
        repeat_segments: (0..g.seqs.len())
            .filter(|&s| repeat[s])
            .map(|s| g.names[s].clone())
            .collect(),
        hidden_repeats: hidden.iter().map(|&s| g.names[s].clone()).collect(),
        ambiguous_repeats: ambiguous.iter().map(|&s| g.names[s].clone()).collect(),
        phased_paths: {
            let mut v: Vec<Traversal> = phased_reads
                .iter()
                .map(|(w, &n)| Traversal {
                    path: w.iter().map(|&o| g.label(o)).collect::<Vec<_>>().join(" "),
                    reads: n,
                })
                .collect();
            v.sort_by(|a, b| b.reads.cmp(&a.reads).then_with(|| a.path.cmp(&b.path)));
            v
        },
        params: serde_json::json!({
            "k": p.k, "flank": p.flank, "min_anchor": p.min_anchor, "diag_tol": p.diag_tol,
            "min_side_hits": p.min_side_hits, "spacing_tol": p.spacing_tol, "bootstrap": p.bootstrap,
            "max_repeat_len": p.max_repeat_len,
            "repeat_rule": if depths.is_some() { "branching + depth, then phasing" } else { "branching, then phasing" },
        }),
        elapsed_seconds: t0.elapsed().as_secs_f64(),
        peak_rss_mb: peak_rss_mb(),
    };
    fs::create_dir_all(out.parent().unwrap_or(Path::new(".")))?;
    fs::write(out, serde_json::to_string_pretty(&report)?)?;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rng_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut s = seed;
        (0..n)
            .map(|_| b"ACGT"[(splitmix(&mut s) % 4) as usize])
            .collect()
    }

    fn params() -> EvidenceParams {
        EvidenceParams {
            k: 21,
            flank: 500,
            min_anchor: 200,
            diag_tol: 20,
            min_side_hits: 5,
            spacing_tol: 0.05,
            bootstrap: 200,
            seed: 1,
            max_repeat_len: 20_000,
        }
    }

    /// X -> R1 -> R2 -> Y with side links Z -> R1 and R2 -> W, so R1 and R2 both branch.
    fn repeat_chain() -> (Graph, [Vec<u8>; 6]) {
        let [x, r1, r2, y, z, w] = [41, 42, 43, 44, 45, 46].map(|i| rng_seq(2000, i));
        let g = Graph {
            names: ["X", "R1", "R2", "Y", "Z", "W"].map(String::from).to_vec(),
            seqs: vec![
                x.clone(),
                r1.clone(),
                r2.clone(),
                y.clone(),
                z.clone(),
                w.clone(),
            ],
            links: vec![
                ((0, true), (1, true)),
                ((4, true), (1, true)),
                ((1, true), (2, true)),
                ((2, true), (3, true)),
                ((2, true), (5, true)),
            ],
            overlaps: vec![0; 5],
        };
        (g, [x, r1, r2, y, z, w])
    }

    #[test]
    fn repeats_are_branching_segments_within_the_length_limit() {
        let (g, _) = repeat_chain();
        assert_eq!(
            repeat_segments(&g, 20_000, &[true; 5], None),
            vec![false, true, true, false, false, false]
        );
        assert_eq!(repeat_segments(&g, 1000, &[true; 5], None), vec![false; 6]);
        // Without reads on Z -> R1, R1 no longer branches.
        assert_eq!(
            repeat_segments(&g, 20_000, &[true, false, true, true, true], None),
            vec![false, false, true, false, false, false]
        );
    }

    #[test]
    fn a_branching_segment_at_single_copy_depth_is_not_a_repeat() {
        let (g, _) = repeat_chain();
        // R1 at 2x depth stays a repeat; R2 branches too but sits at unique-segment depth.
        let depth = [100.0, 200.0, 105.0, 100.0, 95.0, 100.0];
        assert_eq!(
            repeat_segments(&g, 20_000, &[true; 5], Some(&depth)),
            vec![false, true, false, false, false, false]
        );
        let text = "S\tX\tACGT\tdp:f:10\nS\tR1\tACGT\tKC:i:80\nS\tR2\tACGT\n";
        let tiny = Graph {
            names: ["X", "R1", "R2"].map(String::from).to_vec(),
            seqs: vec![b"ACGT".to_vec(); 3],
            links: vec![],
            overlaps: vec![],
        };
        assert_eq!(
            segment_depths(text, &tiny),
            None,
            "one segment has no depth tag"
        );
        let text = "S\tX\tACGT\tdp:f:10\nS\tR1\tACGT\tKC:i:80\nS\tR2\tACGT\tDP:f:7\n";
        assert_eq!(segment_depths(text, &tiny), Some(vec![10.0, 20.0, 7.0]));
    }

    #[test]
    fn depth_backed_repeats_do_not_depend_on_the_read_span_limit() {
        let (g, _) = repeat_chain();
        let depth = [100.0, 200.0, 200.0, 100.0, 100.0, 100.0];
        // Both repeat segments exceed the topology-only limit. Their copy-depth evidence
        // still excludes them from the single-copy anchors used for representative paths.
        assert_eq!(
            repeat_segments(&g, 1000, &[true; 5], Some(&depth)),
            vec![false, true, true, false, false, false]
        );
    }

    #[test]
    fn a_read_across_two_repeats_gives_one_phased_path() {
        let (g, [x, r1, r2, y, ..]) = repeat_chain();
        let idx = build_index(&g, &params());
        let repeat = repeat_segments(&g, 20_000, &[true; 5], None);
        let read = cat(&[&x[1200..], &r1, &r2, &y[..900]]);
        for read in [read.clone(), revcomp(&read)] {
            let (cs, _) = assign(&g, crossings(&read, &idx, &params()), 0.05);
            let ws = walks(&g, &read_steps(&g, &cs), 0.05);
            assert_eq!(ws.len(), 1);
            let paths = phased_paths(&ws[0], &repeat);
            let fw = vec![(0, true), (1, true), (2, true), (3, true)];
            assert_eq!(paths, vec![canonical_walk(&fw)]);
        }
        // Stopping inside R2 phases nothing: the read never reaches a unique exit.
        let partial = cat(&[&x[1200..], &r1, &r2[..1500]]);
        let (cs, _) = assign(&g, crossings(&partial, &idx, &params()), 0.05);
        let ws = walks(&g, &read_steps(&g, &cs), 0.05);
        assert!(ws.iter().all(|w| phased_paths(w, &repeat).is_empty()));
    }

    /// Repeat R between A..B and C..D: graph links A->R, C->R, R->B, R->D.
    fn repeat_graph() -> (Graph, [Vec<u8>; 5]) {
        let [a, b, c, d, r] = [1, 2, 3, 4, 5].map(|i| rng_seq(2000, i));
        let g = Graph {
            names: ["A", "B", "C", "D", "R"].map(String::from).to_vec(),
            seqs: vec![a.clone(), b.clone(), c.clone(), d.clone(), r.clone()],
            links: vec![
                ((0, true), (4, true)),
                ((2, true), (4, true)),
                ((4, true), (1, true)),
                ((4, true), (3, true)),
            ],
            overlaps: vec![0; 4],
        };
        (g, [a, b, c, d, r])
    }

    /// S -> U1 and S -> U2, where U1 and U2 begin with `shared` bp of identical sequence.
    fn twin_exits(shared: usize) -> (Graph, Vec<u8>, Vec<u8>) {
        let s = rng_seq(3000, 31);
        let common = rng_seq(shared, 32);
        let u1 = cat(&[&common, &rng_seq(3000, 33)]);
        let u2 = cat(&[&common, &rng_seq(3000, 34)]);
        let g = Graph {
            names: ["S", "U1", "U2"].map(String::from).to_vec(),
            seqs: vec![s.clone(), u1.clone(), u2],
            links: vec![((0, true), (1, true)), ((0, true), (2, true))],
            overlaps: vec![0, 0],
        };
        (g, s, u1)
    }

    #[test]
    fn read_goes_to_the_exit_it_continues_into() {
        // shared 300 bp < flank: the read's hits beyond the shared part decide.
        let (g, s, u1) = twin_exits(300);
        let idx = build_index(&g, &params());
        let read = cat(&[&s[1500..], &u1[..1500]]);
        let (kept, ambiguous) = assign(&g, crossings(&read, &idx, &params()), 0.05);
        assert_eq!(kept.iter().map(|c| c.link).collect::<Vec<_>>(), vec![0]);
        assert!(ambiguous.is_empty());
    }

    #[test]
    fn a_short_exit_inside_a_longer_one_is_decided_by_the_next_crossing() {
        // S -> U1 and S -> U2, where U1 (1200 bp) is a prefix of U2; U1 continues into V.
        // The S-joint flanks cannot tell them apart, but the read walks through U1 into V.
        let s = rng_seq(3000, 51);
        let u1 = rng_seq(1200, 52);
        let u2 = cat(&[&u1, &rng_seq(2000, 53)]);
        let v = rng_seq(3000, 54);
        let g = Graph {
            names: ["S", "U1", "U2", "V"].map(String::from).to_vec(),
            seqs: vec![s.clone(), u1.clone(), u2, v.clone()],
            links: vec![
                ((0, true), (1, true)),
                ((0, true), (2, true)),
                ((1, true), (3, true)),
            ],
            overlaps: vec![0, 0, 0],
        };
        let idx = build_index(&g, &params());
        let read = cat(&[&s[1500..], &u1, &v[..1500]]);
        for read in [read.clone(), revcomp(&read)] {
            let (kept, ambiguous) = assign(&g, crossings(&read, &idx, &params()), 0.05);
            let mut links: Vec<usize> = kept.iter().map(|c| c.link).collect();
            links.sort();
            assert_eq!(links, vec![0, 2], "S->U1 and U1->V");
            assert!(ambiguous.is_empty());
        }
        // Without the continuation into V the read stays ambiguous between U1 and U2.
        let stops = cat(&[&s[1500..], &u1]);
        let (kept, ambiguous) = assign(&g, crossings(&stops, &idx, &params()), 0.05);
        assert!(kept.iter().all(|c| c.link != 1), "never credited to U2");
        let _ = ambiguous;
    }

    #[test]
    fn exits_identical_across_the_flank_are_ambiguous() {
        // shared 2000 bp > flank: nothing in the indexed window separates U1 from U2.
        let (g, s, u1) = twin_exits(2000);
        let idx = build_index(&g, &params());
        let read = cat(&[&s[1500..], &u1[..2500]]);
        let (kept, mut ambiguous) = assign(&g, crossings(&read, &idx, &params()), 0.05);
        assert!(kept.is_empty());
        ambiguous.sort();
        assert_eq!(ambiguous, vec![0, 1]);
    }

    #[test]
    fn overlapping_links_need_reads_through_the_overlap() {
        // X -300M- S -300M- Y, with S reversed in the GFA to exercise orientation.
        let ov = 300;
        let x = rng_seq(3000, 21);
        let body = rng_seq(1500, 22);
        let y_tail = rng_seq(3000, 23);
        let s: Vec<u8> = cat(&[&x[x.len() - ov..], &body, &y_tail[..ov]]);
        let y: Vec<u8> = y_tail.clone();
        let g = Graph {
            names: ["X", "S", "Y"].map(String::from).to_vec(),
            seqs: vec![x.clone(), revcomp(&s), y.clone()],
            links: vec![((0, true), (1, false)), ((1, false), (2, true))],
            overlaps: vec![ov, ov],
        };
        let idx = build_index(&g, &params());
        let read = cat(&[&x[1500..], &body, &y[..1200]]);
        for read in [read.clone(), revcomp(&read)] {
            let cs = crossings(&read, &idx, &params());
            assert_eq!(cs.len(), 2, "both overlapping links supported");
            let tr = traversals(&g, &read_steps(&g, &cs), 0.05);
            assert_eq!(tr.len(), 1, "S spans len(S) between the two joints");
        }
        // A read that stops inside the overlap of S -> Y proves nothing about Y.
        let short = cat(&[&x[1500..], &body, &y[..ov]]);
        let cs = crossings(&short, &idx, &params());
        assert_eq!(cs.iter().map(|c| c.link).collect::<Vec<_>>(), vec![0]);
    }

    fn cat(parts: &[&[u8]]) -> Vec<u8> {
        parts.concat()
    }

    #[test]
    fn spanning_read_supports_both_joints_and_the_traversal() {
        let (g, [a, b, _, _, r]) = repeat_graph();
        let idx = build_index(&g, &params());
        let read = cat(&[&a[1400..], &r, &b[..600]]);
        for read in [read.clone(), revcomp(&read)] {
            let cs = crossings(&read, &idx, &params());
            let mut links: Vec<usize> = cs.iter().map(|c| c.link).collect();
            links.sort();
            assert_eq!(links, vec![0, 2]);
            let tr = traversals(&g, &read_steps(&g, &cs), 0.05);
            assert_eq!(tr.len(), 1);
            let (x, s, y) = tr[0];
            assert!(
                ((x, s, y) == ((0, true), (4, true), (1, true)))
                    || ((x, s, y) == ((1, false), (4, false), (0, false)))
            );
        }
    }

    #[test]
    fn short_anchor_or_wrong_exit_is_not_support() {
        let (g, [a, b, _, _, r]) = repeat_graph();
        let idx = build_index(&g, &params());
        // only 100 bp into B: below min_anchor
        let cs = crossings(&cat(&[&a[1400..], &r[..800]]), &idx, &params());
        assert_eq!(cs.iter().map(|c| c.link).collect::<Vec<_>>(), vec![0]);
        let cs = crossings(&cat(&[&r[1200..], &b[..100]]), &idx, &params());
        assert!(cs.is_empty());
        // R followed by unrelated sequence supports neither exit
        let cs = crossings(&cat(&[&r[1000..], &rng_seq(800, 99)]), &idx, &params());
        assert!(cs.is_empty());
    }

    #[test]
    fn a_hidden_repeat_is_crossed_in_two_major_contexts() {
        // Segments 0..4 are anchors A B C D, 4 is a known repeat R, 5 is S. Reads cross S after
        // R in two contexts: A R S B and C R S D.
        let o = |s: usize| (s, true);
        let mut phased: BTreeMap<Vec<Oriented>, u64> = BTreeMap::new();
        phased.insert(vec![o(0), o(4), o(5), o(1)], 20);
        phased.insert(vec![o(2), o(4), o(5), o(3)], 15);
        assert!(crossed_in_two_contexts(&phased, 5, &[true; 6]));
        // the same seen from the other strand is still two contexts
        let mut rc: BTreeMap<Vec<Oriented>, u64> = BTreeMap::new();
        rc.insert(vec![o(0), o(4), o(5), o(1)], 20);
        rc.insert(vec![(3, false), (5, false), (4, false), (2, false)], 15);
        assert!(crossed_in_two_contexts(&rc, 5, &[true; 6]));
        // the second context between shallow flanks is a low-frequency configuration
        let mut shallow = [true; 6];
        shallow[2] = false;
        assert!(!crossed_in_two_contexts(&phased, 5, &shallow));
        // too few reads, or one context only
        phased.insert(vec![o(2), o(4), o(5), o(3)], 2);
        assert!(!crossed_in_two_contexts(&phased, 5, &[true; 6]));
    }

    #[test]
    fn minority_share_brackets_the_estimate() {
        let (est, ci) = minority_share(30, 70, 500, 7);
        let (lo, hi) = ci.unwrap();
        assert!((est.unwrap() - 0.3).abs() < 1e-9);
        assert!(lo < 0.3 && hi > 0.3 && lo > 0.15 && hi < 0.45);
    }

    #[test]
    fn soft_masked_reads_count_like_uppercase() {
        let (g, [a, b, _, _, r]) = repeat_graph();
        let idx = build_index(&g, &params());
        let read = cat(&[&a[1400..], &r, &b[..600]]).to_ascii_lowercase();
        let rc = revcomp(&read).to_ascii_lowercase();
        for read in [read, rc] {
            assert_eq!(crossings(&read, &idx, &params()).len(), 2);
        }
    }

    #[test]
    fn repeat_calls_take_the_median_depth_and_reports_keep_the_mean() {
        let text = "S\ta\tACGT\tdp:f:172.0\tdm:f:100.0\nS\tb\tACGT\tdp:f:100.0\n";
        let g = Graph {
            names: vec!["a".into(), "b".into()],
            seqs: vec![b"ACGT".to_vec(), b"ACGT".to_vec()],
            links: vec![],
            overlaps: vec![],
        };
        assert_eq!(repeat_depths(text, &g), Some(vec![100.0, 100.0]));
        assert_eq!(segment_depths(text, &g), Some(vec![172.0, 100.0]));
    }

    #[test]
    fn two_well_read_ways_at_both_ends_make_a_repeat_whatever_the_depth() {
        // X between A/B on its left and C/D on its right
        let g = Graph {
            names: ["A", "B", "X", "C", "D"].map(String::from).to_vec(),
            seqs: vec![b"ACGT".to_vec(); 5],
            links: vec![
                ((0, true), (2, true)),
                ((1, true), (2, true)),
                ((2, true), (3, true)),
                ((2, true), (4, true)),
            ],
            overlaps: vec![0; 4],
        };
        assert_eq!(
            crossroads(&g, &[9, 8, 7, 11]),
            vec![false, false, true, false, false]
        );
        // a minority way at one end (2 of 20 reads) is a minority configuration, not a copy
        assert_eq!(crossroads(&g, &[18, 2, 7, 11]), vec![false; 5]);
        // unread links do not count
        assert_eq!(crossroads(&g, &[9, 0, 7, 11]), vec![false; 5]);
    }
}
