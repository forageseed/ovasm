//! Closing a linear short-read molecule into its circle with read pairs, as GetOrganelle's
//! assemblies close on paired reads.
//!
//! Some sequence is all but missing from short reads: ZM4's plastid has an AT-rich inverted
//! repeat (a hairpin) where reads fall to a tenth of the plastid's depth, and the de novo graph
//! stops on both sides of it at every k, leaving a linear 125,122 bp molecule. Its fragments are
//! still sequenced from both ends: 21 read pairs have one mate running off the molecule's end and
//! the other onto its start, 175 bp apart by their insert size, and the reads anchored at either
//! end, with their mates, spell the 178 bp between (base by base, by majority) until they reach
//! the start. The circle of 125,300 bp matches the held-out HiFi plastid (125,509 bp) at 99.9%.
//!
//! When the reads do not reach the start, the pairs still say the end and the start are joined
//! and how far apart (insert size from the pairs within each end); the reference graph's walks
//! may then hold the sequence between, which `reference_fill` takes only when its length agrees
//! with that distance. Those bases are the reference's, not the reads', and are labelled so.
//! Mates share their read name's first word.

use std::path::PathBuf;

use anyhow::Result;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;

use crate::evidence::revcomp;
use crate::kmer::CanonicalKmers;
use crate::recruit::for_each_rec;

/// Word length of placing reads and of the extension.
const K: usize = 31;
/// Each end of the molecule this long holds the reads placed at it (more than most inserts:
/// ZM4's are 225-445 bp, 5-95%).
const WINDOW: usize = 600;
/// Fewest read pairs joining the end to the start for the circle to be tried.
pub const MIN_JOINING_PAIRS: usize = 5;
/// Fewest local reads behind each base of the extension, and the share of them it needs.
const MIN_SUPPORT: usize = 3;
const DOMINANCE: f64 = 0.8;
/// Longest sequence spelled between the end and the start.
const MAX_FILL: usize = 2000;
/// Bases of the extension a read must hold before the next one, longest first (`close`).
const CONTEXTS: [usize; 4] = [90, 70, 50, K];

#[derive(Debug, Clone, Serialize, Default)]
pub struct Closure {
    /// Read pairs with one mate forward in the last `WINDOW` bp and the other reverse in the
    /// first: fragments across the join.
    pub joining_pairs: usize,
    /// Reads placed at either end, and with their mates the reads the extension spells from.
    pub anchored_names: usize,
    pub local_reads: usize,
    /// Median insert of the pairs with both mates, facing, within one end.
    pub insert_median: Option<f64>,
    /// Median distance between the end and the start by the joining pairs and that insert.
    pub gap_estimate: Option<f64>,
    /// Bases spelled between the end and the start; None when the reads did not reach it.
    pub fill: Option<usize>,
    /// Why the extension stopped, when it did not reach the start.
    pub stopped: Option<String>,
    #[serde(skip)]
    pub sequence: Vec<u8>,
}

fn name_of(id: &[u8]) -> Vec<u8> {
    let name = id.split(|&b| b == b' ' || b == b'\t').next().unwrap_or(id);
    name.strip_suffix(b"/1").or_else(|| name.strip_suffix(b"/2")).unwrap_or(name).to_vec()
}

/// Close `molecule` (linear) with the read pairs of `reads`; `Closure::fill` is Some when the
/// reads spell the join, and `sequence` is then what goes between the end and the start.
pub fn close(molecule: &[u8], reads: &[PathBuf]) -> Result<Closure> {
    let mut c = Closure::default();
    if molecule.len() < 2 * WINDOW {
        c.stopped = Some("molecule shorter than both windows".into());
        return Ok(c);
    }
    let molecule = molecule.to_ascii_uppercase();
    let len = molecule.len() as i64;
    // k-mers of the two windows, either strand: (at the end?, forward?, molecule position), None
    // for a word with two places. ZM4's plastid ends at a hairpin: the end's words are the
    // start's reverse complements, and placing reads by them put the end's forward reads at the
    // start reversed (1 joining pair instead of 21).
    let mut placed: FxHashMap<Vec<u8>, Option<(bool, bool, i64)>> = FxHashMap::default();
    for (offset, at_end) in [(0, false), (molecule.len() - WINDOW, true)] {
        for (j, x) in molecule[offset..offset + WINDOW].windows(K).enumerate() {
            let at = (offset + j) as i64;
            for (word, fw) in [(x.to_vec(), true), (revcomp(x), false)] {
                placed
                    .entry(word)
                    .and_modify(|p| {
                        if *p != Some((at_end, fw, at)) {
                            *p = None;
                        }
                    })
                    .or_insert(Some((at_end, fw, at)));
            }
        }
    }
    // each name: its reads placed, as (at the end?, forward?, leftmost, rightmost) on the molecule
    type Placed = (bool, bool, i64, i64);
    let mut where_: FxHashMap<Vec<u8>, Vec<Placed>> = FxHashMap::default();
    for path in reads {
        for_each_rec(path, |r| {
            let seq = r.seq.to_ascii_uppercase();
            let n = seq.len() as i64;
            let hit = seq
                .windows(K)
                .enumerate()
                .find_map(|(i, x)| placed.get(x).copied().flatten().map(|p| (i as i64, p)));
            if let Some((i, (at_end, fw, at))) = hit {
                let left = if fw { at - i } else { at - (n - i - K as i64) };
                where_
                    .entry(name_of(&r.id))
                    .or_default()
                    .push((at_end, fw, left, left + n));
            }
        })?;
    }
    // Ambiguous duplicate/colliding IDs must not contribute to local extension either:
    // otherwise they could amplify bases after other pairs establish an apparent join.
    where_.retain(|_, placements| placements.len() <= 2);
    c.anchored_names = where_.len();
    let median = |mut v: Vec<i64>| -> Option<f64> {
        v.sort_unstable();
        (!v.is_empty()).then(|| v[v.len() / 2] as f64)
    };
    let mut inserts = Vec::new();
    let mut spans = Vec::new();
    for v in where_.values() {
        // A read name identifies at most one pair, not a set of independent fragments.
        // More than two placements mean duplicate/colliding names: abstain instead of
        // turning their Cartesian product into spurious support for a join.
        if v.len() != 2 { continue; }
        for a in v.iter().filter(|a| a.1) {
            for b in v.iter().filter(|b| !b.1) {
                if a.0 == b.0 && b.3 > a.2 && b.3 - a.2 < 2000 {
                    inserts.push(b.3 - a.2);
                }
                if a.0 && !b.0 {
                    // forward off the end, reverse onto the start: on the molecule are its last
                    // len - left bases and its first `right`
                    spans.push((len - a.2) + b.3);
                }
            }
        }
    }
    c.joining_pairs = spans.len();
    c.insert_median = median(inserts);
    if let (Some(ins), Some(span)) = (c.insert_median, median(spans)) {
        c.gap_estimate = Some(ins - span as f64);
    }
    if c.joining_pairs < MIN_JOINING_PAIRS {
        c.stopped = Some(format!("{} read pairs join the ends", c.joining_pairs));
        return Ok(c);
    }
    // the anchored reads and their mates, both strands
    let mut local: Vec<Vec<u8>> = Vec::new();
    for path in reads {
        for_each_rec(path, |r| {
            if where_.contains_key(&name_of(&r.id)) {
                let seq = r.seq.to_ascii_uppercase();
                local.push(revcomp(&seq));
                local.push(seq);
            }
        })?;
    }
    c.local_reads = local.len() / 2;
    let mut index: FxHashMap<&[u8], Vec<(usize, usize)>> = FxHashMap::default();
    for (i, r) in local.iter().enumerate() {
        for (p, x) in r.windows(K).enumerate() {
            index.entry(x).or_default().push((i, p));
        }
    }
    let target = &molecule[..K];
    // the extension starts with the molecule's last `CONTEXTS[0]` bases (it is longer than both windows)
    let init = CONTEXTS[0];
    let mut cur = molecule[molecule.len() - init..].to_vec();
    let mut seen: FxHashSet<Vec<u8>> = FxHashSet::default();
    while cur.len() - init < MAX_FILL + K {
        let suffix = &cur[cur.len() - K..];
        if !seen.insert(suffix.to_vec()) {
            c.stopped = Some(format!("the extension loops after {} bp", cur.len() - init));
            return Ok(c);
        }
        // the next base after the longest context enough reads hold: across an inverted repeat
        // longer than K, the last K bases are followed by the repeat's two other sides
        let hits = index.get(suffix).map(Vec::as_slice).unwrap_or(&[]);
        let mut next = [0usize; 4];
        for &context in &CONTEXTS {
            if context > cur.len() {
                continue;
            }
            let before = &cur[cur.len() - context..cur.len() - K];
            next = [0; 4];
            for &(i, p) in hits {
                let r = &local[i];
                if p < before.len() || &r[p - before.len()..p] != before {
                    continue;
                }
                if let Some(b) = r.get(p + K).and_then(|b| b"ACGT".iter().position(|x| x == b)) {
                    next[b] += 1;
                }
            }
            if next.iter().sum::<usize>() >= MIN_SUPPORT {
                break;
            }
        }
        let total: usize = next.iter().sum();
        let (b, &n) = next.iter().enumerate().max_by_key(|x| x.1).unwrap();
        if total == 0 || n < MIN_SUPPORT || (n as f64) < DOMINANCE * total as f64 {
            c.stopped = Some(format!(
                "after {} bp the local reads continue A {} C {} G {} T {}",
                cur.len() - init,
                next[0],
                next[1],
                next[2],
                next[3]
            ));
            return Ok(c);
        }
        cur.push(b"ACGT"[b]);
        if cur.ends_with(target) {
            if cur.len() < init + K {
                // the start is reached inside the molecule's own last bases: its ends overlap,
                // which spelling forward does not close
                c.stopped = Some("the start's first bases are the end's last: the ends overlap".into());
                return Ok(c);
            }
            c.sequence = cur[init..cur.len() - K].to_vec();
            c.fill = Some(c.sequence.len());
            return Ok(c);
        }
    }
    c.stopped = Some(format!("no join within {MAX_FILL} bp"));
    Ok(c)
}

/// Fewest homopolymer-compressed units of the molecule's end that its start must repeat for an
/// overlap closure (`overlap_join`): 12 units are at least 12 bases with no two alike in a row,
/// found by chance in a 600 bp window about once in 300 molecules, and the read pairs' distance
/// must agree as well.
const MIN_OVERLAP_UNITS: usize = 12;
/// The bases an overlap closure removes may differ from the read pairs' (negative) end-to-start
/// distance by at most this (rice: 45 against 53, 35 against 49, 75 against 80).
const OVERLAP_GAP_TOLERANCE: f64 = 100.0;

/// An overlap closure: the circle is `molecule[head..len - tail]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct OverlapJoin {
    pub head: usize,
    pub tail: usize,
    /// Homopolymer-compressed units the end and the start share.
    pub units: usize,
}

/// Homopolymer runs of `s`: (base, first position).
fn runs(s: &[u8]) -> Vec<(u8, usize)> {
    let mut out: Vec<(u8, usize)> = Vec::new();
    for (i, &b) in s.iter().enumerate() {
        if out.last().map_or(true, |r| r.0 != b) {
            out.push((b, i));
        }
    }
    out
}

/// Close a linear molecule whose start repeats its end: the short-read graph stops on both sides
/// of a homopolymer it cannot count, and the stretch beside it, assembled from both ends, comes
/// out twice with the homopolymer at two lengths (rice plastids SAMEA115396711, 717 and 721 at
/// one poly(A) tract: the read pairs put the start 49-80 bp before the end, and spelling forward
/// stops at once where the reads differ on the next base). The longest run of the end's last
/// homopolymer-compressed units that occurs once in the start's window (and not at its last
/// unit) is the overlap, of at least `MIN_OVERLAP_UNITS` units. Each homopolymer is taken where
/// both its sides are in the molecule: the overlap's last unit from the start, the others from
/// the end. The bases removed must agree with `gap` (the pairs' end-to-start distance, negative
/// for an overlap) within `OVERLAP_GAP_TOLERANCE`.
pub fn overlap_join(molecule: &[u8], gap: f64) -> Option<OverlapJoin> {
    let m = molecule.to_ascii_uppercase();
    let n = m.len();
    if n < 2 * WINDOW || gap >= 0.0 {
        return None;
    }
    let tail_off = n - WINDOW;
    let (hr, tr) = (runs(&m[..WINDOW]), runs(&m[tail_off..]));
    let hb: Vec<u8> = hr.iter().map(|r| r.0).collect();
    let tb: Vec<u8> = tr.iter().map(|r| r.0).collect();
    // the end's first overlap unit must have a base before it in the window
    for units in (MIN_OVERLAP_UNITS..tb.len()).rev() {
        let suffix = &tb[tb.len() - units..];
        // the start's last overlap unit must have a base after it
        let mut at = (0..hb.len().saturating_sub(units)).filter(|&j| &hb[j..j + units] == suffix);
        let (Some(j), None) = (at.next(), at.next()) else {
            continue;
        };
        let head = hr[j + units - 1].1;
        let t_last = tail_off + tr[tr.len() - 1].1;
        if t_last <= head {
            return None;
        }
        let removed = n - (t_last - head);
        if (removed as f64 + gap).abs() > OVERLAP_GAP_TOLERANCE {
            return None;
        }
        return Some(OverlapJoin { head, tail: n - t_last, units });
    }
    None
}

/// Most bases a molecule's end (or start) may be trimmed back to find an exact anchor in a
/// reference walk; the walk must then spell those bases as the molecule does.
const MAX_TRIM: usize = 200;
/// Longest reference sequence taken between the end and the start.
pub const MAX_REFERENCE_FILL: usize = 5000;
/// The reference's sequence between the anchors must be this close to the pairs' distance:
/// within this many bases, or this share of the distance, whichever is more.
const FILL_TOLERANCE_BP: f64 = 100.0;
const FILL_TOLERANCE_SHARE: f64 = 0.3;

#[derive(Debug, Clone, Serialize)]
pub struct ReferenceFill {
    /// The walk (index) and strand the fill was taken from.
    pub walk: usize,
    pub reverse: bool,
    pub length: usize,
    /// Distinct fills of an agreeing length the walks held (more than one: another was as good).
    pub candidates: usize,
    /// Share of the molecule's sampled k-mers the walk holds (`similarity`): the fill is taken
    /// from the agreeing walk most like the individual.
    pub similarity: f64,
    /// The most similar walk with a different agreeing fill, and its similarity.
    pub runner_up: Option<(usize, f64)>,
    #[serde(skip)]
    pub sequence: Vec<u8>,
}

/// The sequence a reference walk (`walks`: sequence, circular) holds between the molecule's end
/// and its start, when its length agrees with `gap` (the pairs' distance). Anchors are exact
/// `K`-mers of the molecule's end and start, moved inwards up to `MAX_TRIM` bases where the
/// individual differs from every walk at the very end. None when no walk agrees.
pub fn reference_fill(molecule: &[u8], walks: &[(Vec<u8>, bool)], gap: f64) -> Option<ReferenceFill> {
    let m = molecule.to_ascii_uppercase();
    if m.len() < 2 * (K + MAX_TRIM) {
        return None;
    }
    let n = m.len();
    let mut anchors: FxHashMap<&[u8], Vec<(bool, usize)>> = FxHashMap::default();
    for d in 0..=MAX_TRIM {
        anchors.entry(&m[n - K - d..n - d]).or_default().push((true, d));
        anchors.entry(&m[d..d + K]).or_default().push((false, d));
    }
    let tolerance = FILL_TOLERANCE_BP.max(FILL_TOLERANCE_SHARE * gap.abs());
    let mut found: Vec<(f64, Vec<u8>, usize, bool)> = Vec::new();
    for (w, (seq, circular)) in walks.iter().enumerate() {
        for reverse in [false, true] {
            let mut s = if reverse { revcomp(seq) } else { seq.to_ascii_uppercase() };
            if *circular {
                let wrap: Vec<u8> = s.iter().cycle().take(MAX_REFERENCE_FILL + 2 * (K + MAX_TRIM)).copied().collect();
                s.extend(wrap);
            }
            let (mut ends, mut starts) = (Vec::new(), Vec::new());
            for (p, x) in s.windows(K).enumerate() {
                for &(is_end, d) in anchors.get(x).map(Vec::as_slice).unwrap_or(&[]) {
                    if is_end {
                        ends.push((d, p));
                    } else {
                        starts.push((d, p));
                    }
                }
            }
            ends.sort_unstable();
            starts.sort_unstable_by_key(|&(e, q)| (q, e));
            for &(d, p) in &ends {
                let from = p + K + d;
                if s.get(p + K..from) != Some(&m[n - d..]) {
                    continue;
                }
                // the nearest start after it whose trimmed bases the walk spells too
                let next = starts.iter().find(|&&(e, q)| {
                    q >= from + e && q - e - from <= MAX_REFERENCE_FILL && s[q - e..q] == m[..e]
                });
                if let Some(&(e, q)) = next {
                    let fill = s[from..q - e].to_vec();
                    let off = (fill.len() as f64 - gap).abs();
                    if off <= tolerance {
                        found.push((off, fill, w, reverse));
                    }
                    break;
                }
            }
        }
    }
    if found.is_empty() {
        return None;
    }
    // of the walks whose fill agrees, the one most like the individual, then the closest length
    let like = similarity(&m, walks);
    found.sort_by(|a, b| like[b.2].total_cmp(&like[a.2]).then(a.0.total_cmp(&b.0)));
    let distinct: FxHashSet<&[u8]> = found.iter().map(|f| f.1.as_slice()).collect();
    let candidates = distinct.len();
    if candidates != 1 { return None; }
    // the most similar walk whose fill differs from the chosen one
    let runner_up = found
        .iter()
        .find(|f| f.1 != found[0].1)
        .map(|f| (f.2, like[f.2]));
    let (_, sequence, walk, reverse) = found.swap_remove(0);
    Some(ReferenceFill {
        walk,
        reverse,
        length: sequence.len(),
        candidates,
        similarity: like[walk],
        runner_up,
        sequence,
    })
}

/// Share of the molecule's sampled `K`-mers (every `SIMILARITY_STEP`th, either strand) each walk
/// holds: how much the individual is like each reference run.
fn similarity(molecule: &[u8], walks: &[(Vec<u8>, bool)]) -> Vec<f64> {
    let sample: FxHashSet<u64> = CanonicalKmers::new(molecule, K)
        .filter(|(p, _)| p % SIMILARITY_STEP == 0)
        .map(|(_, x)| x)
        .collect();
    walks
        .iter()
        .map(|(seq, _)| {
            let held: FxHashSet<u64> = CanonicalKmers::new(seq, K)
                .map(|(_, x)| x)
                .filter(|x| sample.contains(x))
                .collect();
            held.len() as f64 / sample.len().max(1) as f64
        })
        .collect()
}

/// One sampled k-mer in this many for `similarity`.
const SIMILARITY_STEP: usize = 16;

/// A piece placed on the walk: its span there differs from its length by at most this many
/// bases, or this share of its length.
const PLACE_TOLERANCE_BP: usize = 100;
const PLACE_TOLERANCE_SHARE: f64 = 0.05;
/// Share of the pieces' bases that must place, and of the walk they must cover, for a scaffold.
const MIN_PLACED_SHARE: f64 = 0.9;
/// How far into a piece its first (last) K-mer the walk holds is looked for: the ends of a
/// piece are where the reads stopped, and can differ from every walk (a quarter of the piece at
/// most).
const SCAFFOLD_SEARCH: usize = 2000;
const MIN_WALK_COVER: f64 = 0.5;

#[derive(Debug, Clone, Serialize)]
pub struct ScaffoldPlacement {
    pub piece: usize,
    pub reverse: bool,
    /// Span on the walk (end may pass the walk's length: across its origin).
    pub walk_start: usize,
    pub walk_end: usize,
    /// Where its bases stand in the scaffold, and how many of its first bases were left out
    /// because the placement before covered them; `contained`: none of it was used.
    pub start: usize,
    pub end: usize,
    pub trimmed: usize,
    pub contained: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct ScaffoldGap {
    /// Placements (indices) on either side, the walk's bases taken, and where they stand.
    pub after: usize,
    pub before: usize,
    pub walk_start: usize,
    pub length: usize,
    pub start: usize,
    /// Filled with the individual's own bases: the reverse complement of the inverted copy the
    /// scaffold holds from its pieces (the other inverted repeat), not the walk's.
    pub mirrored: bool,
    /// Of its bases, those the walk gave (all of them unless mirrored).
    pub walk_bp: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct Scaffold {
    pub walk: usize,
    pub similarity: f64,
    pub walk_length: usize,
    pub placements: Vec<ScaffoldPlacement>,
    pub gaps: Vec<ScaffoldGap>,
    pub unplaced: Vec<usize>,
    /// Bases of the walk's own in the scaffold (gaps not mirrored), and of the individual's
    /// inverted copies filled into gaps.
    pub fill_bp: usize,
    pub mirrored_bp: usize,
    #[serde(skip)]
    pub sequence: Vec<u8>,
}

/// A circular molecule's pieces laid onto the circular reference walk most like them (`walks`:
/// sequence, circular), for a molecule as conserved in structure as a plastid, whose reads leave
/// it in pieces (1001 Genomes Arabidopsis at 100 bp: nine linear pieces of 127 kb, one copy of
/// the inverted repeat among them, for a 154 kb circle). Each piece is placed where the walk
/// holds its first and last `K`-mers at a distance its length allows (moved inwards up to `SCAFFOLD_SEARCH` bases), on either strand and as often
/// as the walk holds it: a collapsed inverted
/// repeat places twice, once each way, so the circle has both. In walk order, a placement
/// overlapping the one before loses the bases it repeats; the walk's own bases fill between
/// them. Those bases are the reference's, not the reads': the gaps are reported. None when no
/// walk holds the pieces, less than `MIN_PLACED_SHARE` of their bases place, or the
/// placements cover less than `MIN_WALK_COVER` of the walk.
pub fn panel_scaffold(pieces: &[Vec<u8>], walks: &[(Vec<u8>, bool)]) -> Option<Scaffold> {
    let pieces: Vec<Vec<u8>> = pieces.iter().map(|p| p.to_ascii_uppercase()).collect();
    let all: Vec<u8> = pieces.iter().flat_map(|p| p.iter().copied().chain(std::iter::once(b'N'))).collect();
    let like = similarity(&all, walks);
    // a panel run walk of a circular molecule spells the whole circle whether or not it is
    // marked circular (GFA P lines of a merged panel are not)
    let w = (0..walks.len()).filter(|&i| like[i] > 0.0).max_by(|&a, &b| like[a].total_cmp(&like[b]).then(b.cmp(&a)))?;
    let walk = walks[w].0.to_ascii_uppercase();
    let l = walk.len();
    if l < 2 * K {
        return None;
    }
    let doubled: Vec<u8> = [walk.as_slice(), walk.as_slice()].concat();
    let mut at: FxHashMap<&[u8], Vec<usize>> = FxHashMap::default();
    for (p, x) in doubled.windows(K).enumerate() {
        at.entry(x).or_default().push(p);
    }
    let mut placed: Vec<(usize, bool, usize, usize)> = Vec::new();
    let mut placed_bp = 0usize;
    let mut unplaced = Vec::new();
    for (i, piece) in pieces.iter().enumerate() {
        let n = piece.len();
        let mut any = false;
        if n >= 2 * (K + MAX_TRIM) {
            let search = SCAFFOLD_SEARCH.min(n / 4 - K);
            let tolerance = PLACE_TOLERANCE_BP.max((PLACE_TOLERANCE_SHARE * n as f64) as usize);
            for reverse in [false, true] {
                let s = if reverse { revcomp(piece) } else { piece.clone() };
                // the first start K-mer the walk holds, and the first end K-mer
                let start = (0..=search).find_map(|d| at.get(&s[d..d + K]).map(|v| (d, v)));
                let end = (0..=search).find_map(|e| at.get(&s[n - K - e..n - e]).map(|v| (e, v)));
                let (Some((d, starts)), Some((e, ends))) = (start, end) else { continue };
                for &p in starts.iter().filter(|&&p| p >= d && p - d < l) {
                    let ws = p - d;
                    let best = ends
                        .iter()
                        .map(|&q| q + K + e)
                        .filter(|&we| we > ws && we - ws <= 2 * l)
                        .min_by_key(|&we| (we - ws).abs_diff(n));
                    if let Some(we) = best.filter(|&we| (we - ws).abs_diff(n) <= tolerance) {
                        placed.push((i, reverse, ws, we));
                        any = true;
                    }
                }
            }
        }
        if any {
            placed_bp += n;
        } else {
            unplaced.push(i);
        }
    }
    let total: usize = pieces.iter().map(Vec::len).sum();
    if placed.is_empty() || (placed_bp as f64) < MIN_PLACED_SHARE * total as f64 {
        return None;
    }
    placed.sort_by_key(|&(i, r, ws, we)| (ws, std::cmp::Reverse(we), i, r));
    placed.dedup_by_key(|x| (x.2, x.3));
    let origin = placed[0].2;
    let mut sequence: Vec<u8> = Vec::new();
    let mut placements: Vec<ScaffoldPlacement> = Vec::new();
    let mut gaps: Vec<ScaffoldGap> = Vec::new();
    let mut cursor = 0usize; // walk bases covered, from the origin
    let mut last_used: Option<usize> = None;
    for &(i, reverse, ws, we) in &placed {
        let (rs, re) = (ws - origin, we - origin);
        if re <= cursor {
            placements.push(ScaffoldPlacement { piece: i, reverse, walk_start: ws, walk_end: we, start: sequence.len(), end: sequence.len(), trimmed: 0, contained: true });
            continue;
        }
        if rs > cursor {
            let from = (origin + cursor) % l;
            let fill: Vec<u8> = doubled[from..from + (rs - cursor)].to_vec();
            gaps.push(ScaffoldGap { after: last_used.unwrap_or(usize::MAX), before: placements.len(), walk_start: from, length: fill.len(), start: sequence.len(), mirrored: false, walk_bp: 0 });
            sequence.extend_from_slice(&fill);
        }
        let trimmed = cursor.saturating_sub(rs).min(pieces[i].len());
        let s = if reverse { revcomp(&pieces[i]) } else { pieces[i].clone() };
        let start = sequence.len();
        sequence.extend_from_slice(&s[trimmed..]);
        placements.push(ScaffoldPlacement { piece: i, reverse, walk_start: ws, walk_end: we, start, end: sequence.len(), trimmed, contained: false });
        last_used = Some(placements.len() - 1);
        cursor = re;
    }
    let covered: usize = cursor.min(l);
    // around the origin: the walk's bases up to the first placement, or the last placement's
    // bases that the first already holds left out
    if cursor < l {
        let from = (origin + cursor) % l;
        let fill: Vec<u8> = doubled[from..from + (l - cursor)].to_vec();
        gaps.push(ScaffoldGap { after: last_used.unwrap_or(usize::MAX), before: 0, walk_start: from, length: fill.len(), start: sequence.len(), mirrored: false, walk_bp: 0 });
        sequence.extend_from_slice(&fill);
    } else if cursor > l {
        let over = (cursor - l).min(sequence.len());
        sequence.truncate(sequence.len() - over);
    }
    let fill_bp: usize = gaps.iter().map(|g| g.length).sum();
    if ((covered.saturating_sub(fill_bp.min(covered))) as f64) < MIN_WALK_COVER * l as f64 {
        return None;
    }
    gaps.iter_mut().for_each(|g| if g.after == usize::MAX { g.after = last_used.unwrap_or(0) });
    let mirrored_bp = mirror_gaps(&mut sequence, &mut gaps);
    let fill_bp = gaps.iter().map(|g| g.walk_bp).sum();
    Some(Scaffold { walk: w, similarity: like[w], walk_length: l, placements, gaps, unplaced, fill_bp, mirrored_bp, sequence })
}

/// Gaps of a scaffold whose walk bases are the reverse complement of a stretch the individual's
/// own pieces hold elsewhere in it (the other copy of an inverted repeat, which the reads gave
/// once): filled with that stretch's reverse complement instead. The stretch is found by the
/// reverse complements of the gap's last and first `K`-mers it holds (moved inwards up to `SCAFFOLD_SEARCH`
/// bases), at the gap's length (within the placement tolerance), wholly outside every gap.
/// Returns the bases so filled; gap starts are moved to the new coordinates.
fn mirror_gaps(sequence: &mut Vec<u8>, gaps: &mut [ScaffoldGap]) -> usize {
    let n = sequence.len();
    let mut at: FxHashMap<&[u8], Vec<usize>> = FxHashMap::default();
    for (p, x) in sequence.windows(K).enumerate() {
        at.entry(x).or_default().push(p);
    }
    let in_gap = |a: usize, b: usize| gaps.iter().any(|g| a < g.start + g.length && g.start < b);
    let mut fills: Vec<(usize, Vec<u8>, usize)> = Vec::new();
    for (gi, g) in gaps.iter().enumerate() {
        let len = g.length;
        if len < 2 * (K + MAX_TRIM) || g.start + len > n {
            continue;
        }
        let gap = &sequence[g.start..g.start + len];
        let rcg = revcomp(gap);
        let search = SCAFFOLD_SEARCH.min(len / 4 - K);
        // the gap's ends can be single-copy sequence beside the inverted repeat: the copy is
        // looked for from the first and last K-mers it holds, inwards
        let first = (0..=search).find_map(|d| at.get(&rcg[d..d + K]).map(|v| (d, v.clone())));
        let last = (0..=search).find_map(|e| at.get(&rcg[len - K - e..len - e]).map(|v| (e, v.clone())));
        let (Some((d, firsts)), Some((e, lasts))) = (first, last) else { continue };
        let core = len - d - e;
        let tolerance = PLACE_TOLERANCE_BP.max((PLACE_TOLERANCE_SHARE * core as f64) as usize);
        let found = firsts.iter().find_map(|&a| {
            lasts
                .iter()
                .map(|&q| q + K)
                .filter(|&b| b > a && b <= n && (b - a).abs_diff(core) <= tolerance && !in_gap(a, b))
                .min_by_key(|&b| (b - a).abs_diff(core))
                .map(|b| (a, b))
        });
        if let Some((a, b)) = found {
            // rcg[d..len-e] is the reverse complement of gap[e..len-d]: that core takes the
            // reverse complement of the copy, the walk's bases stay either side of it
            let fill: Vec<u8> = gap[..e].iter().copied().chain(revcomp(&sequence[a..b])).chain(gap[len - d..].iter().copied()).collect();
            fills.push((gi, fill, b - a));
        }
    }
    drop(at);
    let mut mirrored = 0usize;
    // replace from the last gap back, so earlier starts stay; then move the starts
    fills.sort_by_key(|f| std::cmp::Reverse(gaps[f.0].start));
    for (gi, fill, copied) in &fills {
        let g = &gaps[*gi];
        sequence.splice(g.start..g.start + g.length, fill.iter().copied());
        mirrored += copied;
    }
    let walk_left: FxHashMap<usize, usize> = fills.iter().map(|(gi, f, c)| (*gi, f.len() - c)).collect();
    let mut delta: i64 = 0;
    let mut order: Vec<usize> = (0..gaps.len()).collect();
    order.sort_by_key(|&i| gaps[i].start);
    let filled: FxHashMap<usize, usize> = fills.iter().map(|(gi, f, _)| (*gi, f.len())).collect();
    for i in order {
        let g = &mut gaps[i];
        g.start = (g.start as i64 + delta) as usize;
        g.walk_bp = g.length;
        if let Some(&fl) = filled.get(&i) {
            delta += fl as i64 - g.length as i64;
            g.length = fl;
            g.mirrored = true;
            g.walk_bp = walk_left[&i];
        }
    }
    mirrored
}

#[cfg(test)]
mod scaffold_tests {
    use super::*;

    fn random(n: usize, seed: &mut u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                *seed ^= *seed << 13;
                *seed ^= *seed >> 7;
                *seed ^= *seed << 17;
                b"ACGT"[(*seed % 4) as usize]
            })
            .collect()
    }

    /// A plastid: LSC, IRb, SSC, IRa (IRa the reverse complement of IRb).
    fn plastid(seed: u64) -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>) {
        let mut s = seed;
        let (lsc, ir, ssc) = (random(6000, &mut s), random(1500, &mut s), random(1200, &mut s));
        let circle = [lsc.clone(), ir.clone(), ssc.clone(), revcomp(&ir)].concat();
        (circle, lsc, ir, ssc)
    }

    #[test]
    fn pieces_with_one_copy_of_the_inverted_repeat_scaffold_into_the_whole_circle() {
        let (relative, lsc, ir, ssc) = plastid(5);
        // the individual: a few substitutions, its pieces short of their ends by 40-60 bases,
        // the inverted repeat once
        let mut ind_lsc = lsc.clone();
        ind_lsc[3000] = if ind_lsc[3000] == b'A' { b'C' } else { b'A' };
        let pieces = vec![ind_lsc[50..5950].to_vec(), revcomp(&ir[40..1460]), ssc[60..1140].to_vec()];
        let walks = vec![(relative.clone(), true), (random(9000, &mut 9), true)];
        let s = panel_scaffold(&pieces, &walks).expect("a scaffold");
        assert_eq!(s.walk, 0);
        assert!(s.unplaced.is_empty(), "{s:?}");
        assert_eq!(s.sequence.len(), relative.len(), "{s:?}");
        // the inverted repeat placed twice, once each way; the substitution kept
        assert_eq!(s.placements.iter().filter(|p| p.piece == 1 && !p.contained).count(), 2, "{s:?}");
        let doubled = [s.sequence.clone(), s.sequence.clone()].concat();
        let probe = &ind_lsc[2980..3020];
        assert!(doubled.windows(probe.len()).any(|w| w == probe) || doubled.windows(probe.len()).any(|w| w == revcomp(probe).as_slice()));
        // the reference's bases: the trimmed ends only
        assert_eq!(s.fill_bp, relative.len() - (5900 + 2 * 1420 + 1080), "{:?}", s.gaps);
        assert_eq!(s.gaps.len(), 4);
    }

    #[test]
    fn the_inverted_repeat_given_once_fills_its_other_copy_with_the_individuals_bases() {
        let (relative, lsc, ir, ssc) = plastid(9);
        // the individual differs from the relative inside the inverted repeat (both copies)
        let mut ind_ir = ir.clone();
        ind_ir[700] = if ind_ir[700] == b'G' { b'T' } else { b'G' };
        // one piece IRb + SSC (so it places once), and the LSC
        let pieces = vec![lsc[20..5980].to_vec(), [ind_ir.clone(), ssc[..1100].to_vec()].concat()];
        let s = panel_scaffold(&pieces, &[(relative.clone(), true)]).expect("a scaffold");
        assert_eq!(s.sequence.len(), relative.len(), "{s:?}");
        assert!(s.gaps.iter().any(|g| g.mirrored && g.length > 1000), "{:?}", s.gaps);
        // both copies carry the individual's base: IRb forward, IRa as its reverse complement
        let probe = &ind_ir[680..720];
        let doubled = [s.sequence.clone(), s.sequence.clone()].concat();
        let count = |x: &[u8]| doubled[..s.sequence.len() + 40].windows(x.len()).filter(|w| *w == x).count();
        assert_eq!(count(probe) + count(&revcomp(probe)), 2, "{:?}", s.gaps);
        assert!(s.fill_bp < 300, "{s:?}");
    }

    #[test]
    fn pieces_the_walks_do_not_hold_are_not_scaffolded() {
        let (relative, _, _, _) = plastid(5);
        let mut s = 77u64;
        let pieces = vec![random(3000, &mut s), random(2000, &mut s)];
        assert!(panel_scaffold(&pieces, &[(relative, true)]).is_none());
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn an_end_repeated_at_the_start_with_its_homopolymer_at_another_length_closes() {
        // circle: X + A(17) + U + Y, assembled from just after a poly(A) at the start (21 As,
        // the run cut) to the end of U at the end (17 As, both sides there)
        let mut x = 7u64;
        let mut rnd = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    x ^= x << 13;
                    x ^= x >> 7;
                    x ^= x << 17;
                    b"ACGT"[(x % 4) as usize]
                })
                .collect()
        };
        let u = b"GAATTACATTAAATATAGAATATA".to_vec();
        let (body, gatt) = (rnd(5000), b"GATT".to_vec());
        let start = [vec![b'A'; 21], u.clone()].concat();
        let end = [gatt.clone(), vec![b'A'; 17], u.clone()].concat();
        let molecule = [start.clone(), body.clone(), end.clone()].concat();
        let j = overlap_join(&molecule, -53.0).expect("closes");
        let circle = &molecule[j.head..molecule.len() - j.tail];
        // one copy of U, the poly(A) at the end's length, U's last base once
        assert_eq!(circle.len(), body.len() + gatt.len() + 17 + u.len());
        assert_eq!(&circle[..1], &u[u.len() - 1..]);
        assert!(circle.ends_with(&[gatt.as_slice(), &[b'A'; 17], &u[..u.len() - 1]].concat()));
        // the pairs say the ends are far apart: no closure
        assert!(overlap_join(&molecule, -400.0).is_none());
        assert!(overlap_join(&molecule, 30.0).is_none());
        // nothing repeated: no closure
        assert!(overlap_join(&[rnd(700), rnd(5000)].concat(), -50.0).is_none());
    }

    use super::*;
    use std::fmt::Write as _;

    fn random(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = x;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                b"ACGT"[((z ^ (z >> 31)) >> 62) as usize]
            })
            .collect()
    }

    /// Pairs of 100 bp mates of 330 bp fragments every `step` bases round a circular genome,
    /// skipping reads that would lie wholly inside `[hole.0, hole.1)` (sequencing dropout).
    fn pairs(genome: &[u8], step: usize, hole: (usize, usize), path: &std::path::Path) {
        let n = genome.len();
        let circ = |a: usize, len: usize| -> Vec<u8> { (0..len).map(|i| genome[(a + i) % n]).collect() };
        let inside = |a: usize| a >= hole.0 && a + 100 <= hole.1;
        let mut fq = String::new();
        for (f, a) in (0..n).step_by(step).enumerate() {
            let m1 = circ(a, 100);
            let m2 = revcomp(&circ(a + 230, 100));
            for (m, s, pos) in [(1, m1, a), (2, m2, a + 230)] {
                if inside(pos % n) {
                    continue;
                }
                let s = String::from_utf8(s).unwrap();
                writeln!(fq, "@f{f} {m}\n{s}\n+\n{}", "I".repeat(s.len())).unwrap();
            }
        }
        std::fs::write(path, fq).unwrap();
    }

    #[test]
    fn read_pairs_spell_the_join_of_a_circle() {
        let genome = random(5000, 3);
        let dir = std::env::temp_dir().join(format!("ovasm-close-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let reads = dir.join("reads.fq");
        pairs(&genome, 4, (4900, 5000), &reads);
        // the de novo molecule stops 120 bp short of its own start
        let molecule = &genome[..4880];
        let c = close(molecule, &[reads.clone()]).unwrap();
        assert!(c.joining_pairs >= MIN_JOINING_PAIRS, "{c:?}");
        assert_eq!(c.fill, Some(120), "{c:?}");
        assert_eq!(c.sequence, genome[4880..].to_vec());
        // a molecule whose ends no fragment joins stays linear
        let half = [genome[..2500].to_vec(), random(2500, 4)].concat();
        let c = close(&half, &[reads]).unwrap();
        assert_eq!(c.fill, None);
        assert!(c.joining_pairs < MIN_JOINING_PAIRS, "{c:?}");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn slash_mate_names_preserve_closure_and_pair_evidence() {
        let genome = random(5000, 103);
        let dir = std::env::temp_dir().join(format!("ovasm-close-suffixed-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let original = dir.join("original.fq");
        let suffixed = dir.join("suffixed.fq");
        pairs(&genome, 4, (4900, 5000), &original);
        let input = std::fs::read_to_string(&original).unwrap();
        let mut rewritten = String::new();
        for (i, line) in input.lines().enumerate() {
            if i % 4 == 0 {
                let (id, mate) = line.split_once(' ').unwrap();
                writeln!(rewritten, "{id}/{mate} instrument").unwrap();
            } else {
                writeln!(rewritten, "{line}").unwrap();
            }
        }
        std::fs::write(&suffixed, rewritten).unwrap();
        let ordinary = close(&genome[..4880], &[original]).unwrap();
        let normalized = close(&genome[..4880], &[suffixed]).unwrap();
        assert!(ordinary.joining_pairs >= MIN_JOINING_PAIRS, "{ordinary:?}");
        assert_eq!(normalized.joining_pairs, ordinary.joining_pairs);
        assert_eq!(normalized.fill, Some(120));
        assert_eq!(normalized.sequence, genome[4880..].to_vec());
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn duplicate_read_names_do_not_multiply_closing_support() {
        let genome = random(5000, 123);
        let dir = std::env::temp_dir().join(format!("ovasm-close-duplicate-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let reads = dir.join("reads.fq");
        let m1 = &genome[4800..4900];
        let m2 = revcomp(&genome[30..130]);
        let mut fq = String::new();
        for _ in 0..5 {
            for (mate, seq) in [(1, m1), (2, m2.as_slice())] {
                writeln!(fq, "@same/{mate}\n{}\n+\n{}", String::from_utf8_lossy(seq), "I".repeat(seq.len())).unwrap();
            }
        }
        std::fs::write(&reads, fq).unwrap();
        let result = close(&genome[..4880], &[reads]).unwrap();
        assert_eq!(result.joining_pairs, 0, "duplicate IDs cannot establish independent pairs: {result:?}");
        assert_eq!(result.anchored_names, 0, "ambiguous IDs must not supply extension reads");
        assert_eq!(result.fill, None);
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn a_join_between_the_arms_of_an_inverted_repeat_is_spelled() {
        // the circle is S core X H with S the reverse complement of X (60 bp, longer than K):
        // the molecule S core X lacks H, unread, between the repeat's arms (a hairpin)
        let x = random(60, 21);
        let core = random(4800, 22);
        let h = random(70, 23);
        let genome = [revcomp(&x), core.clone(), x.clone(), h.clone()].concat();
        let dir = std::env::temp_dir().join(format!("ovasm-close-ir-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let reads = dir.join("reads.fq");
        let n = genome.len();
        pairs(&genome, 3, (n - 70, n), &reads);
        let molecule = &genome[..n - 70];
        let c = close(molecule, &[reads]).unwrap();
        assert!(c.joining_pairs >= MIN_JOINING_PAIRS, "{c:?}");
        assert_eq!(c.sequence, h, "{c:?}");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn the_pairs_measure_the_gap() {
        let genome = random(5000, 5);
        let dir = std::env::temp_dir().join(format!("ovasm-close-gap-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let reads = dir.join("reads.fq");
        // fragments of 330 bp; nothing is read inside the last 150 bp of the circle
        pairs(&genome, 4, (4850, 5100), &reads);
        let c = close(&genome[..4850], &[reads]).unwrap();
        assert_eq!(c.insert_median, Some(330.0), "{c:?}");
        let gap = c.gap_estimate.unwrap();
        assert!((gap - 150.0).abs() <= 10.0, "{c:?}");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn a_reference_walk_fills_a_gap_of_the_measured_length() {
        let genome = random(6000, 7);
        // the individual differs from the walk in its last 20 bases before the gap: trimmed back
        let mut molecule = genome[..5800].to_vec();
        let n = molecule.len();
        molecule[n - 10] = if molecule[n - 10] == b'A' { b'C' } else { b'A' };
        let walks = vec![(genome.clone(), true)];
        let f = reference_fill(&molecule, &walks, 200.0);
        assert!(f.is_none(), "an end that differs from the walk is not spelled by it");
        let molecule = &genome[..5800];
        let f = reference_fill(molecule, &walks, 190.0).unwrap();
        assert_eq!((f.length, f.candidates), (200, 1));
        assert_eq!(f.sequence, genome[5800..].to_vec());
        // the walk's distance must agree with the pairs'
        assert!(reference_fill(molecule, &walks, 900.0).is_none());
        // reverse strand of a walk
        let f = reference_fill(molecule, &[(revcomp(&genome), true)], 200.0).unwrap();
        assert!(f.reverse && f.sequence == genome[5800..].to_vec());
    }

    #[test]
    fn conflicting_reference_fills_are_not_resolved_by_global_similarity() {
        let individual = random(6000, 11);
        // run 0 differs from the individual at every 50th base (away from the gap's flanks) and
        // holds another 200 bp in the gap; run 1 is the individual's own sequence
        let mut far = individual.clone();
        for i in (300..5000).step_by(50) {
            far[i] = if far[i] == b'A' { b'C' } else { b'A' };
        }
        let other = random(200, 12);
        far[5800..].copy_from_slice(&other);
        let near = individual.clone();
        let walks = vec![(far, true), (near, true)];
        assert!(reference_fill(&individual[..5800], &walks, 200.0).is_none(),
            "global similarity must not resolve conflicting local gap alleles");
    }
}
