//! Read crossings of junctions cut from one backbone, from a single index of the backbone.
//!
//! The candidate junctions of `junctions::restore` are joins of slices of the backbone: every
//! segment of the hypothesis graph is a forward slice of a backbone segment, and the two segments
//! of a link share the repeat (at least the seed word, longer than k). Every junction k-mer is
//! then a backbone k-mer, and a read's hits on a junction are its backbone hits inside the
//! junction's window, moved to junction coordinates. `crossings` on the per-junction index
//! (`build_index`) gets the same hits by listing each backbone k-mer once per junction window
//! that holds it: with thousands of short inverted repeats that is hundreds of times (Clitopilus
//! prunulus mitochondrion, ERR13660011, 125 kb: 2,650 repeat pairs, 10,600 junctions, 43 million
//! junction k-mers; 290-480 s of a 650 s run on 16 threads). Here a read is looked up once and
//! its backbone hits chained once, per strand; each chain is cut into pieces whose diagonals
//! span at most `diag_tol` (any subset of a piece is one chain), and a junction is answered from
//! range counts on the pieces in its window, joined into chains in diagonal order. The crossings
//! are those of `crossings` on the per-junction index, ordered by (link, orientation).

use rustc_hash::{FxHashMap, FxHashSet};

use super::{
    forward_is_canonical, junction_pieces, revcomp, BestChain, Crossing, EvidenceParams, Graph,
    Junction, Oriented,
};
use crate::kmer::CanonicalKmers;

/// A segment of the hypothesis graph as the forward slice of a backbone segment from `start`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct SliceSource {
    pub segment: usize,
    pub start: usize,
}

/// Junction k-mer starts q in [q0, q1] that are the backbone k-mers of `segment` starting at
/// P = q + shift (`forward`) or P = shift - q (the junction shows the reverse strand there).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Run {
    segment: usize,
    q0: i64,
    q1: i64,
    forward: bool,
    shift: i64,
}

impl Run {
    fn q(&self, p: i64) -> i64 {
        if self.forward {
            p - self.shift
        } else {
            self.shift - p
        }
    }

    fn p(&self, q: i64) -> i64 {
        if self.forward {
            q + self.shift
        } else {
            self.shift - q
        }
    }

    /// Backbone k-mer starts of the junction k-mer starts [a, b] in this run.
    fn p_range(&self, a: i64, b: i64) -> Option<(i64, i64)> {
        let (a, b) = (a.max(self.q0), b.min(self.q1));
        (a <= b).then(|| {
            let (x, y) = (self.p(a), self.p(b));
            (x.min(y), x.max(y))
        })
    }
}

/// A read k-mer at `rpos` that is the backbone k-mer starting at `p`; `same`: both show the
/// strand of the canonical k-mer, or both its reverse.
#[derive(Clone, Copy)]
struct Hit {
    p: i64,
    rpos: i64,
    same: bool,
}

/// A junction of one forward slice: backbone k-mer starts [`first`, `last`], P = q + `shift`.
#[derive(Clone, Copy)]
struct Window {
    first: i64,
    last: i64,
    shift: i64,
    junction: u32,
}

pub(crate) struct Layout {
    k: usize,
    junctions: Vec<Junction>,
    runs: Vec<Vec<Run>>,
    /// backbone canonical k-mer -> (segment, start, the k-mer's forward strand is canonical)
    index: FxHashMap<u64, Vec<(u32, u32, bool)>>,
    /// junctions of one forward slice, per backbone segment, by `first`
    windows: FxHashMap<usize, Vec<Window>>,
    longest: FxHashMap<usize, i64>,
    /// the other junctions
    split: Vec<u32>,
}

/// Bases [from, to) of oriented segment `o` of `g`.
fn oriented_slice(g: &Graph, o: Oriented, from: usize, to: usize) -> Vec<u8> {
    let s = &g.seqs[o.0];
    if o.1 {
        s[from..to].to_vec()
    } else {
        revcomp(&s[s.len() - to..s.len() - from])
    }
}

/// The same bases to `CanonicalKmers` and `forward_is_canonical`: equal ignoring case, or both
/// not ACGT (a k-mer over either is skipped). The overlap of a repeat's inverted join is found
/// from canonical words, so its two copies can differ in case.
fn same_kmer_bases(x: &[u8], y: &[u8]) -> bool {
    let code = |b: u8| match b.to_ascii_uppercase() {
        c @ (b'A' | b'C' | b'G' | b'T') => c,
        _ => 0,
    };
    x.len() == y.len() && x.iter().zip(y).all(|(&a, &b)| code(a) == code(b))
}

/// Junction k-mer starts [q0, q1] whose bases are oriented segment `o` from index q + `offset`.
fn run(g: &Graph, sources: &[SliceSource], o: Oriented, offset: i64, q0: i64, q1: i64, k: usize) -> Run {
    let src = sources[o.0];
    let (start, len) = (src.start as i64, g.seqs[o.0].len() as i64);
    // forward: k-mer [i, i + k) of the segment starts at backbone start + i; reverse: the
    // segment's reverse strand from i is the backbone k-mer ending at start + len - i
    let (forward, shift) = if o.1 {
        (true, start + offset)
    } else {
        (false, start + len - offset - k as i64)
    };
    Run { segment: src.segment, q0, q1, forward, shift }
}

impl Layout {
    /// None where the junctions are not cut from `backbone` as `sources` says (the caller then
    /// uses the per-junction index).
    pub(crate) fn new(
        g: &Graph,
        backbone: &Graph,
        sources: &[SliceSource],
        p: &EvidenceParams,
    ) -> Option<Self> {
        Self::build(g, backbone, sources, p).ok()
    }

    /// `new`, or why not.
    fn build(g: &Graph, backbone: &Graph, sources: &[SliceSource], p: &EvidenceParams) -> Result<Self, String> {
        let k = p.k;
        // odd k: no k-mer is its own reverse complement, so a reversed slice flips the strand
        if k % 2 == 0 || sources.len() != g.seqs.len() {
            return Err(format!("k {k}, {} sources for {} segments", sources.len(), g.seqs.len()));
        }
        for (i, (seq, s)) in g.seqs.iter().zip(sources).enumerate() {
            let slice = backbone.seqs.get(s.segment).and_then(|b| b.get(s.start..s.start + seq.len()));
            if slice != Some(&seq[..]) {
                return Err(format!("segment {i} is not backbone {}:{}", s.segment, s.start));
            }
        }
        let ki = k as i64;
        let (mut junctions, mut runs) = (Vec::new(), Vec::new());
        for (junction, piece) in junction_pieces(g, p) {
            // a side shorter than `min_anchor` goes on into the segments beyond it (`onward` in
            // `build_index`): several junctions per link, which only the index path builds
            if (piece.left_len - piece.overlap).min(piece.right_len) < p.min_anchor {
                return Err(format!("link {}: a side shorter than min_anchor goes onward", junction.link));
            }
            let (left_len, ov) = (piece.left_len as i64, piece.overlap as i64);
            let last = left_len + piece.right_len as i64 - ki;
            // k-mers inside the left part, and the rest: right of the joint or across it,
            // where the left part ends in the overlap, i.e. the right segment's beginning
            let a = run(g, sources, piece.a, piece.left_skip as i64, 0, left_len - ki, k);
            let b_first = (left_len - ki + 1).max(0);
            let b = run(g, sources, piece.b, ov - left_len, b_first, last, k);
            let mut list = Vec::new();
            if a.q0 <= a.q1 {
                list.push(a);
            }
            if b.q0 <= b.q1 {
                let across = left_len - b_first;
                if across > ov {
                    return Err(format!("link {}: overlap {ov} shorter than k - 1", junction.link));
                }
                if across > 0 {
                    let la = g.seqs[piece.a.0].len();
                    let tail = oriented_slice(g, piece.a, la - across as usize, la);
                    let head = oriented_slice(g, piece.b, (ov - across) as usize, ov as usize);
                    if !same_kmer_bases(&tail, &head) {
                        return Err(format!(
                            "link {}: {} then {} differ in the overlap",
                            junction.link,
                            String::from_utf8_lossy(&tail),
                            String::from_utf8_lossy(&head)
                        ));
                    }
                }
                match list.last_mut() {
                    Some(prev)
                        if prev.forward
                            && b.forward
                            && prev.segment == b.segment
                            && prev.shift == b.shift =>
                    {
                        prev.q1 = b.q1
                    }
                    _ => list.push(b),
                }
            }
            junctions.push(junction);
            runs.push(list);
        }
        let used: FxHashSet<usize> = sources.iter().map(|s| s.segment).collect();
        let mut index: FxHashMap<u64, Vec<(u32, u32, bool)>> = FxHashMap::default();
        for &segment in &used {
            let seq = &backbone.seqs[segment];
            for (pos, kmer) in CanonicalKmers::new(seq, k) {
                let fwd = forward_is_canonical(&seq[pos..pos + k], kmer);
                index.entry(kmer).or_default().push((segment as u32, pos as u32, fwd));
            }
        }
        let (mut windows, mut split): (FxHashMap<usize, Vec<Window>>, Vec<u32>) = Default::default();
        for (ji, list) in runs.iter().enumerate() {
            match list.as_slice() {
                [] => {}
                [r] if r.forward => windows.entry(r.segment).or_default().push(Window {
                    first: r.p(r.q0),
                    last: r.p(r.q1),
                    shift: r.shift,
                    junction: ji as u32,
                }),
                _ => split.push(ji as u32),
            }
        }
        let mut longest = FxHashMap::default();
        for (&segment, w) in windows.iter_mut() {
            w.sort_unstable_by_key(|w| (w.first, w.junction));
            longest.insert(segment, w.iter().map(|w| w.last - w.first).max().unwrap_or(0));
        }
        Ok(Layout { k, junctions, runs, index, windows, longest, split })
    }

    /// The crossings of `super::crossings` on the per-junction index of the same graph, by
    /// (link, orientation).
    pub(crate) fn crossings(&self, read: &[u8], p: &EvidenceParams) -> Vec<Crossing> {
        let k = self.k as i64;
        let mut hits: FxHashMap<usize, Vec<Hit>> = FxHashMap::default();
        for (rpos, kmer) in CanonicalKmers::new(read, self.k) {
            let Some(entries) = self.index.get(&kmer) else {
                continue;
            };
            let rfwd = forward_is_canonical(&read[rpos..rpos + self.k], kmer);
            for &(segment, pos, bfwd) in entries {
                hits.entry(segment as usize).or_default().push(Hit {
                    p: pos as i64,
                    rpos: rpos as i64,
                    same: rfwd == bfwd,
                });
            }
        }
        // per segment and strand (`Hit::same`), the read's chains
        let mut chains: FxHashMap<(usize, bool), Chains> = FxHashMap::default();
        for (&segment, hs) in &hits {
            for same in [true, false] {
                let c = Chains::new(hs, same, p.diag_tol);
                if !c.e.is_empty() {
                    chains.insert((segment, same), c);
                }
            }
        }
        let mut best: FxHashMap<(u32, bool), BestChain> = FxHashMap::default();
        let mut scratch = Scratch::default();
        for (&(segment, same), c) in &chains {
            if let Some(windows) = self.windows.get(&segment) {
                self.forward_windows(c, same, windows, self.longest[&segment], k, p, &mut best, &mut scratch);
            }
        }
        let mut out = Vec::new();
        for &ji in &self.split {
            for along in [true, false] {
                out.extend(self.split_junction(ji, along, &chains, k, p, &mut scratch));
            }
        }
        for ((ji, along), b) in &best {
            out.extend(b.crossing(&self.junctions[*ji as usize], *along));
        }
        out.sort_unstable_by_key(|c| (c.link, c.along));
        out
    }

    /// One orientation (`same` strand as the backbone: along) of the junctions of one forward
    /// slice. Their hits are the read's backbone hits in the window with diagonal
    /// key + const (along: key = rpos - P, + shift; against: key = rpos + P, - shift + k) and
    /// q = P - shift. Chains of a junction's hits never join hits of two chains of the whole read
    /// (consecutive keys there differ by more than `diag_tol`), so a chain of the read with fewer
    /// than 2 `min_side_hits` hits cannot anchor in any window, and each chain of the read is
    /// cut to the window on its own (`Chains::items`, `offer_groups`).
    #[allow(clippy::too_many_arguments)]
    fn forward_windows(
        &self,
        c: &Chains,
        same: bool,
        windows: &[Window],
        longest: i64,
        k: i64,
        p: &EvidenceParams,
        best: &mut FxHashMap<(u32, bool), BestChain>,
        s: &mut Scratch,
    ) {
        let along = same;
        for chain in &c.chains {
            if chain.range.len() < 2 * p.min_side_hits {
                continue;
            }
            let from = windows.partition_point(|w| w.first < chain.pmin - longest);
            for w in &windows[from..] {
                if w.first > chain.pmax {
                    break;
                }
                if w.last < chain.pmin {
                    continue;
                }
                let junction = &self.junctions[w.junction as usize];
                let joint = junction.left_len as i64;
                let unique_left_end = joint - junction.overlap as i64;
                let constant = if along { w.shift } else { k - w.shift };
                s.items.clear();
                s.diagonals.clear();
                // left: q + k <= unique_left_end; right: q >= joint
                c.items(
                    chain,
                    (w.first, w.last),
                    (w.shift, w.shift + unique_left_end - k),
                    (w.shift + joint, w.last),
                    constant,
                    |pp| pp - w.shift,
                    s,
                );
                if s.items.is_empty() {
                    continue;
                }
                // the pieces of one chain come in key order
                let entry = best.entry((w.junction, along)).or_default();
                offer_groups(s, junction, k, p, entry);
            }
        }
    }

    /// One orientation of a junction of two slices (or of one reversed slice). Each slice maps
    /// the read's hits on one strand to junction diagonals key + const, the pieces of the read's
    /// chains cut to the slice; a chain on the junction is a run of pieces from either slice.
    fn split_junction(
        &self,
        ji: u32,
        along: bool,
        chains: &FxHashMap<(usize, bool), Chains>,
        k: i64,
        p: &EvidenceParams,
        s: &mut Scratch,
    ) -> Option<Crossing> {
        let junction = &self.junctions[ji as usize];
        let runs = &self.runs[ji as usize];
        let joint = junction.left_len as i64;
        let unique_left_end = joint - junction.overlap as i64;
        let last = runs.iter().map(|r| r.q1).max()?;
        s.items.clear();
        s.diagonals.clear();
        for r in runs {
            // a forward slice shows the read's `along` strand as the backbone's, a reversed one
            // the other
            let Some(c) = chains.get(&(r.segment, along == r.forward)) else { continue };
            let Some(window) = r.p_range(r.q0, r.q1) else { continue };
            let constant = match (r.forward, along) {
                (true, true) => r.shift,
                (true, false) => k - r.shift,
                (false, true) => -r.shift,
                (false, false) => r.shift + k,
            };
            let none = (1, 0);
            let left = r.p_range(0, unique_left_end - k).unwrap_or(none);
            let right = r.p_range(joint, last).unwrap_or(none);
            c.over(window.0, window.1, |chain| c.items(chain, window, left, right, constant, |pp| r.q(pp), s));
        }
        // an anchored chain has min_side_hits on each side
        let (l, rt) = s.items.iter().fold((0, 0), |(a, b), x| (a + x.left as usize, b + x.right as usize));
        if l < p.min_side_hits || rt < p.min_side_hits {
            return None;
        }
        s.items.sort_unstable_by_key(|x| (x.lo, x.hi));
        let mut best = BestChain::default();
        offer_groups(s, junction, k, p, &mut best);
        best.crossing(junction, along)
    }
}

/// Offer the chains of one junction and orientation to `best`, in diagonal order: runs of the
/// items (sorted by lowest diagonal) without a gap over `diag_tol`. Each item is one chain by
/// itself, so in the sorted hits of the junction these runs are exactly `diagonal_chains`.
fn offer_groups(s: &Scratch, junction: &Junction, k: i64, p: &EvidenceParams, best: &mut BestChain) {
    let joint = junction.left_len as i64;
    let unique_left_end = joint - junction.overlap as i64;
    let items = &s.items;
    let mut start = 0;
    while start < items.len() {
        let mut end = start + 1;
        let mut high = items[start].hi;
        while end < items.len() && items[end].lo - high <= p.diag_tol {
            high = high.max(items[end].hi);
            end += 1;
        }
        let group = &items[start..end];
        start = end;
        let left: usize = group.iter().map(|x| x.left as usize).sum();
        let right: usize = group.iter().map(|x| x.right as usize).sum();
        if left < p.min_side_hits || right < p.min_side_hits {
            continue;
        }
        let min_j = group.iter().map(|x| x.min_q).min().expect("nonempty group");
        let max_j = group.iter().map(|x| x.max_q).max().expect("nonempty group") + k;
        // `anchored_sides`
        if unique_left_end - min_j < junction.anchor_left || max_j - joint < junction.anchor_right {
            continue;
        }
        let n: usize = group.iter().map(|x| x.n as usize).sum();
        best.offer(junction, k, (left, right, min_j, max_j), || {
            let mut d: Vec<(i64, usize)> = group
                .iter()
                .flat_map(|x| s.diagonals[x.diagonals.start as usize..x.diagonals.end as usize].iter().copied())
                .collect();
            d.sort_unstable();
            // the diagonal of hit n / 2 of the chain in (diagonal, q) order
            let mut rank = n / 2;
            for (diag, m) in d {
                if rank < m {
                    return diag;
                }
                rank -= m;
            }
            unreachable!("the middle hit is one of the chain's {n}")
        });
    }
}

/// Chains spanning at most this many backbone bases are found by their first start.
const SHORT_CHAIN: i64 = 4096;

/// Hits of the read whose keys span at most `diag_tol`: any subset of them is one chain.
/// `stretches` index `Chains::stretches` (a key's hits each, sorted by P), in key order.
struct Piece {
    stretches: std::ops::Range<u32>,
    pmin: i64,
    pmax: i64,
}

/// A chain of the read's hits on one backbone segment and strand: `Chains::e[range]`, cut into
/// `Chains::pieces[pieces]` in key order.
struct ChainInfo {
    range: std::ops::Range<usize>,
    pieces: std::ops::Range<u32>,
    pmin: i64,
    pmax: i64,
}

/// The read's hits on one backbone segment and strand (`Hit::same`) as (key, P): key = rpos - P
/// on the backbone's strand, rpos + P on the other (a junction's diagonal is key + a constant),
/// sorted, and their chains of keys at most `diag_tol` apart, in key order.
struct Chains {
    e: Vec<(i64, i64)>,
    stretches: Vec<std::ops::Range<usize>>,
    pieces: Vec<Piece>,
    chains: Vec<ChainInfo>,
    short: Vec<(i64, u32)>,
    long: Vec<u32>,
}

impl Chains {
    fn new(hs: &[Hit], same: bool, tol: i64) -> Self {
        let mut e: Vec<(i64, i64)> = hs
            .iter()
            .filter(|h| h.same == same)
            .map(|h| (if same { h.rpos - h.p } else { h.rpos + h.p }, h.p))
            .collect();
        e.sort_unstable();
        let (mut stretches, mut pieces, mut chains) = (Vec::new(), Vec::new(), Vec::new());
        let mut start = 0;
        while start < e.len() {
            let mut end = start + 1;
            while end < e.len() && e[end].0 - e[end - 1].0 <= tol {
                end += 1;
            }
            let first_piece = pieces.len() as u32;
            let mut i = start;
            while i < end {
                let base = e[i].0;
                let first_stretch = stretches.len() as u32;
                let (mut pmin, mut pmax) = (i64::MAX, i64::MIN);
                while i < end && e[i].0 - base <= tol {
                    let j = i + e[i..end].partition_point(|x| x.0 == e[i].0);
                    stretches.push(i..j);
                    pmin = pmin.min(e[i].1);
                    pmax = pmax.max(e[j - 1].1);
                    i = j;
                }
                pieces.push(Piece { stretches: first_stretch..stretches.len() as u32, pmin, pmax });
            }
            let cut = &pieces[first_piece as usize..];
            let pmin = cut.iter().map(|x| x.pmin).min().expect("a chain has a piece");
            let pmax = cut.iter().map(|x| x.pmax).max().expect("a chain has a piece");
            chains.push(ChainInfo { range: start..end, pieces: first_piece..pieces.len() as u32, pmin, pmax });
            start = end;
        }
        let (mut short, mut long) = (Vec::new(), Vec::new());
        for (i, c) in chains.iter().enumerate() {
            if c.pmax - c.pmin <= SHORT_CHAIN {
                short.push((c.pmin, i as u32));
            } else {
                long.push(i as u32);
            }
        }
        short.sort_unstable();
        Chains { e, stretches, pieces, chains, short, long }
    }

    /// Each chain with starts on both sides of or inside [lo, hi] (in no particular order).
    fn over(&self, lo: i64, hi: i64, mut each: impl FnMut(&ChainInfo)) {
        let from = self.short.partition_point(|x| x.0 < lo - SHORT_CHAIN);
        for x in self.short[from..].iter().take_while(|x| x.0 <= hi) {
            let c = &self.chains[x.1 as usize];
            if c.pmax >= lo {
                each(c);
            }
        }
        for &i in &self.long {
            let c = &self.chains[i as usize];
            if c.pmin <= hi && c.pmax >= lo {
                each(c);
            }
        }
    }

    /// The pieces of `chain` cut to backbone starts [window.0, window.1], as items (in key
    /// order): hits in `left` and `right` (ranges inside the window; an empty one as (1, 0)),
    /// diagonals key + `constant`, junction positions `q`.
    #[allow(clippy::too_many_arguments)]
    fn items(
        &self,
        chain: &ChainInfo,
        window: (i64, i64),
        left: (i64, i64),
        right: (i64, i64),
        constant: i64,
        q: impl Fn(i64) -> i64,
        s: &mut Scratch,
    ) {
        for piece in &self.pieces[chain.pieces.start as usize..chain.pieces.end as usize] {
            if piece.pmax < window.0 || piece.pmin > window.1 {
                continue;
            }
            let at = s.diagonals.len() as u32;
            let (mut n, mut nl, mut nr) = (0u32, 0u32, 0u32);
            let (mut min_p, mut max_p) = (i64::MAX, i64::MIN);
            for st in &self.stretches[piece.stretches.start as usize..piece.stretches.end as usize] {
                let e = &self.e[st.clone()];
                let lo = e.partition_point(|x| x.1 < window.0);
                let hi = e.partition_point(|x| x.1 <= window.1);
                if lo >= hi {
                    continue;
                }
                n += (hi - lo) as u32;
                min_p = min_p.min(e[lo].1);
                max_p = max_p.max(e[hi - 1].1);
                nl += count_starts(e, left) as u32;
                nr += count_starts(e, right) as u32;
                s.diagonals.push((e[0].0 + constant, hi - lo));
            }
            if n == 0 {
                continue;
            }
            let to = s.diagonals.len() as u32;
            let (qa, qb) = (q(min_p), q(max_p));
            s.items.push(Item {
                lo: s.diagonals[at as usize].0,
                hi: s.diagonals[to as usize - 1].0,
                n,
                left: nl,
                right: nr,
                min_q: qa.min(qb),
                max_q: qa.max(qb),
                diagonals: at..to,
            });
        }
    }
}

/// Hits of a key's stretch (sorted by P) with a start in [a, b]; none for a > b.
fn count_starts(e: &[(i64, i64)], (a, b): (i64, i64)) -> usize {
    e.partition_point(|x| x.1 <= b).saturating_sub(e.partition_point(|x| x.1 < a))
}

/// A piece of a read chain cut to a junction window: diagonals [lo, hi] (each with its hits at
/// `Scratch::diagonals[diagonals]`, in order), hits, hits left and right of the joint, and the
/// lowest and highest junction position.
struct Item {
    lo: i64,
    hi: i64,
    n: u32,
    left: u32,
    right: u32,
    min_q: i64,
    max_q: i64,
    diagonals: std::ops::Range<u32>,
}

/// Buffers reused across the junctions of one read.
#[derive(Default)]
struct Scratch {
    items: Vec<Item>,
    diagonals: Vec<(i64, usize)>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::{build_index, crossings};

    fn rng(seed: &mut u64) -> u64 {
        *seed ^= *seed << 13;
        *seed ^= *seed >> 7;
        *seed ^= *seed << 17;
        *seed
    }

    fn random_seq(n: usize, seed: &mut u64) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[(rng(seed) % 4) as usize]).collect()
    }

    fn key(c: &Crossing) -> (usize, bool, i64, i64, u32, u32) {
        (c.link, c.along, c.to_start, c.from_end, c.score, c.share.to_bits())
    }

    /// A read from `template` at a random place and strand, with substitutions and small indels.
    fn read_from(template: &[u8], len: usize, errors: usize, seed: &mut u64) -> Vec<u8> {
        let len = len.min(template.len());
        let at = (rng(seed) as usize) % (template.len() - len + 1);
        let mut r = template[at..at + len].to_vec();
        for _ in 0..errors {
            let i = (rng(seed) as usize) % r.len();
            match rng(seed) % 3 {
                0 => r[i] = b"ACGT"[(rng(seed) % 4) as usize],
                1 => {
                    r.remove(i);
                }
                _ => r.insert(i, b"ACGT"[(rng(seed) % 4) as usize]),
            }
        }
        if rng(seed) % 2 == 0 {
            r = revcomp(&r);
        }
        r
    }

    /// Backbones with many short inverted repeats (and a tandem copy and an N run): the
    /// backbone layout gives exactly the crossings of the per-junction index, for reads of the
    /// backbone and of the inverted (alternative) form.
    #[test]
    fn backbone_layout_matches_the_junction_index() {
        let p = crate::junctions::hifi_evidence();
        let mut checked = 0;
        for trial in 0..4u64 {
            let mut seed = 0x9e37_79b9_7f4a_7c15 ^ (trial + 1);
            let len = 25_000 + (rng(&mut seed) % 10_000) as usize;
            let mut seq = random_seq(len, &mut seed);
            for _ in 0..(20 + trial as usize * 15) {
                let n = 31 + (rng(&mut seed) % 40) as usize;
                let a = 600 + (rng(&mut seed) as usize) % (len - 2 * 600 - 2 * n - 1200);
                let b = a + n + 1100 + (rng(&mut seed) as usize) % (len - a - n - 1100 - 600 - n).max(1);
                if b + n + 600 >= len {
                    continue;
                }
                let rc = revcomp(&seq[a..a + n]);
                seq[b..b + n].copy_from_slice(&rc);
            }
            // a tandem copy and an N run (assembly graphs are upper case: `inverted_repeats`
            // orders a word's strands by bytes, so in lower case it can pair direct copies,
            // whose joins are no backbone slices and are left to the per-junction index)
            let copy = seq[5_000..5_800].to_vec();
            seq[9_000..9_800].copy_from_slice(&copy);
            for x in &mut seq[15_000..15_010] {
                *x = b'N';
            }
            let backbone = Graph {
                names: vec!["u0".into()],
                seqs: vec![seq.clone()],
                links: Vec::new(),
                overlaps: Vec::new(),
            };
            let repeats = crate::junctions::inverted_repeats(&backbone, 1001, 31, p.min_anchor, false);
            assert!(repeats.len() >= 5, "trial {trial}: {} repeats", repeats.len());
            let (h, sources) = crate::junctions::hypothesis_graph(&backbone, &repeats, p.flank);
            let layout = Layout::build(&h, &backbone, &sources, &p)
                .unwrap_or_else(|e| panic!("trial {trial}: {e}; {:?}", &repeats[..repeats.len().min(3)]));
            assert!(!layout.split.is_empty() && !layout.windows.is_empty());
            let idx = build_index(&h, &p);
            // reads of the backbone and of each repeat's inverted form
            let mut templates = vec![seq.clone()];
            for r in repeats.iter().take(8) {
                let mut alt = seq[..r.left + r.length].to_vec();
                alt.extend_from_slice(&revcomp(&seq[r.left..r.right]));
                alt.extend_from_slice(&seq[r.right + r.length..]);
                templates.push(alt);
            }
            for (t, template) in templates.iter().enumerate() {
                for i in 0..25 {
                    let read = read_from(template, 3_000 + (rng(&mut seed) % 9_000) as usize, i % 7, &mut seed);
                    let mut want: Vec<_> = crossings(&read, &idx, &p).iter().map(key).collect();
                    want.sort_unstable();
                    let got: Vec<_> = layout.crossings(&read, &p).iter().map(key).collect();
                    assert_eq!(got, want, "trial {trial} template {t} read {i}");
                    checked += want.len();
                }
            }
        }
        assert!(checked > 300, "only {checked} crossings compared");
    }

    #[test]
    fn layout_refuses_segments_that_are_not_backbone_slices() {
        let p = crate::junctions::hifi_evidence();
        let mut seed = 7;
        let seq = random_seq(5_000, &mut seed);
        let backbone = Graph { names: vec!["u0".into()], seqs: vec![seq.clone()], links: vec![], overlaps: vec![] };
        let mut other = seq[100..900].to_vec();
        other[400] = if other[400] == b'A' { b'C' } else { b'A' };
        let g = Graph { names: vec!["x".into()], seqs: vec![other], links: vec![], overlaps: vec![] };
        assert!(Layout::new(&g, &backbone, &[SliceSource { segment: 0, start: 100 }], &p).is_none());
        let g = Graph { names: vec!["x".into()], seqs: vec![seq[100..900].to_vec()], links: vec![], overlaps: vec![] };
        assert!(Layout::new(&g, &backbone, &[SliceSource { segment: 0, start: 100 }], &p).is_some());
        let even = EvidenceParams { k: 20, ..p };
        assert!(Layout::new(&g, &backbone, &[SliceSource { segment: 0, start: 100 }], &even).is_none());
    }
}
