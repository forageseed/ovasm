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
    id.split(|&b| b == b' ' || b == b'\t').next().unwrap_or(id).to_vec()
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
    c.anchored_names = where_.len();
    let median = |mut v: Vec<i64>| -> Option<f64> {
        v.sort_unstable();
        (!v.is_empty()).then(|| v[v.len() / 2] as f64)
    };
    let mut inserts = Vec::new();
    let mut spans = Vec::new();
    for v in where_.values() {
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

#[cfg(test)]
mod tests {
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
    fn the_fill_comes_from_the_reference_most_like_the_individual() {
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
        let f = reference_fill(&individual[..5800], &walks, 200.0).unwrap();
        assert_eq!(f.walk, 1);
        assert_eq!(f.sequence, individual[5800..].to_vec());
        assert_eq!(f.candidates, 2);
        let (w, s) = f.runner_up.unwrap();
        assert!(w == 0 && s < f.similarity, "{f:?}");
    }
}
