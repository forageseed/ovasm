//! Frontier extension for reference reconstruction: sequence of the target individual that no
//! reference path holds, found from the WGS k-mer counts before recruitment.
//!
//! Recruitment takes a short read only when its bait covers 80% of it, so each pass moves an edge
//! of the bait by about a fifth of a read, and the saturation rule stops it after the second pass.
//! The ZM4 individual's mitochondrion holds 15 kb that none of the eight reference runs has
//! (two pieces, 8.5 and 6.6 kb, at the same ~370x as the sequence on either side); recruitment
//! stopped at its edges. Extending pass after pass instead (`--saturation 0`, 60 passes) walked
//! into the nucleus at the fifth pass (0.94 M to 4.6 M reads) and closed nothing.
//!
//! So the extension is walked on the k-mer counts of the whole library: one pass counts every
//! 31-mer in a count-min table, the target's depth is the median count of its reference
//! k-mers, and from each place where the reference stops being present in the individual
//! (a reference k-mer at depth followed by one that is not), the walk takes the next k-mers
//! at the target's depth until it rejoins the reference, stops, or runs out of budget. Nuclear
//! sequence is far shallower than an organelle (ZM4: nuclear k-mers at most ~30x against the
//! mitochondrion's ~350x), so the walk cannot step into it as recruitment did; a plastid is far
//! deeper than a mitochondrion, so a mitochondrial walk stops at plastid-like sequence too.
//!
//! The walked sequence only joins the bait. Reads still enter the target by the recruitment
//! rules, and the individual graph still holds only what reads support (projection and de novo
//! assembly); nothing walked here is emitted as sequence.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU16, Ordering};
use std::time::Instant;

use anyhow::Result;
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use serde::Serialize;

use crate::kmer::CanonicalKmers;

/// Word length of the walk: long enough that the short repeats of a plant mitochondrion (most
/// under 30 bp) do not branch it, and a 2-bit word still fits a `u64`.
pub const K: usize = 31;
const MASK: u64 = (1u64 << (2 * K)) - 1;

/// What the walk needs besides the counts.
#[derive(Debug, Clone, Serialize)]
pub struct Params {
    /// The count table holds 2^`slots_log2` 16-bit counters.
    pub slots_log2: u32,
    /// A k-mer is the target's from this share of its depth up. A quarter let the walk into the
    /// nucleus: ZM4's nuclear k-mers at a NUMT's edge count ~30, a collision in the table adds
    /// tens, and a quarter of the mitochondrion's 271 is 68.
    pub floor: f64,
    /// And up to this many times its depth (`recruit::HIGH_COPY_FOLD`): an organelle's own repeat
    /// is two or three times as deep, a plastid twenty times a mitochondrion.
    pub ceiling: f64,
    /// Of two next k-mers at depth, one less than this share of the other's count is taken for a
    /// sequencing error, a collision in the table or the nucleus at a NUMT's edge, not a branch:
    /// a branch of the organelle's own (a repeat's two exits) has both at its depth.
    pub dominance: f64,
    /// Most new k-mers walked from one place.
    pub max_per_end: usize,
    /// Most new k-mers walked for a target, as a multiple of its stored reference bases.
    pub max_fold: f64,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            slots_log2: 31,
            floor: 0.4,
            ceiling: 10.0,
            dominance: 0.5,
            max_per_end: 50_000,
            max_fold: 2.0,
        }
    }
}

/// Counters of one cache line: a k-mer's four counters are in one line, so a count costs one
/// memory access, not four (`KmerCounts`).
#[repr(align(64))]
struct Line([AtomicU16; 32]);

/// Count-min table of canonical `K`-mers with conservative update: a k-mer increments only its
/// smallest counters, so a count overstates by little more than the k-mers sharing all four of
/// them. ZM4: 23 Gb of reads, about 21 G k-mers over 2^31 counters.
pub struct KmerCounts {
    lines: Vec<Line>,
    line_bits: u32,
    pub reads: u64,
    pub bases: u64,
    pub seconds: f64,
}

#[inline]
fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

impl KmerCounts {
    pub fn new(slots_log2: u32) -> Self {
        let line_bits = slots_log2.max(5) - 5;
        Self {
            lines: (0..1usize << line_bits)
                .map(|_| Line(std::array::from_fn(|_| AtomicU16::new(0))))
                .collect(),
            line_bits,
            reads: 0,
            bases: 0,
            seconds: 0.0,
        }
    }

    /// Count every read of `inputs` in one pass.
    pub fn from_reads(inputs: &[PathBuf], slots_log2: u32) -> Result<Self> {
        let start = Instant::now();
        let mut table = Self::new(slots_log2);
        let (batches, reader) = crate::recruit::read_batches(inputs);
        for batch in batches {
            table.reads += batch.len() as u64;
            table.bases += batch.iter().map(|r| r.seq.len() as u64).sum::<u64>();
            batch.par_iter().for_each(|r| table.add_seq(&r.seq));
        }
        reader
            .join()
            .map_err(|_| anyhow::anyhow!("k-mer counting reader panicked"))??;
        table.seconds = start.elapsed().as_secs_f64();
        Ok(table)
    }

    fn add_seq(&self, seq: &[u8]) {
        for (_, key) in CanonicalKmers::new(seq, K) {
            self.add(key);
        }
    }

    /// The line and the four distinct counters of a k-mer.
    #[inline]
    fn slots(&self, key: u64) -> (&[AtomicU16; 32], [usize; 4]) {
        let h = mix64(key);
        let line = &self.lines[(h >> (64 - self.line_bits.max(1))) as usize & (self.lines.len() - 1)].0;
        let (mut slots, mut n, mut bits, mut round) = ([0usize; 4], 0, h, 0u64);
        while n < 4 {
            if bits == 0 {
                round += 1;
                bits = mix64(h ^ round);
            }
            let s = (bits & 31) as usize;
            bits >>= 5;
            if !slots[..n].contains(&s) {
                slots[n] = s;
                n += 1;
            }
        }
        (line, slots)
    }

    #[inline]
    fn add(&self, key: u64) {
        let (line, slots) = self.slots(key);
        let values = slots.map(|s| line[s].load(Ordering::Relaxed));
        let min = *values.iter().min().unwrap();
        if min == u16::MAX {
            return;
        }
        for (&s, &v) in slots.iter().zip(&values) {
            if v == min {
                // a lost race undercounts by one, beside thousands
                let _ = line[s].compare_exchange(min, min + 1, Ordering::Relaxed, Ordering::Relaxed);
            }
        }
    }

    /// Count of a canonical k-mer (at least its occurrences, up to 65,535).
    #[inline]
    pub fn get(&self, key: u64) -> u32 {
        let (line, slots) = self.slots(key);
        u32::from(slots.iter().map(|&s| line[s].load(Ordering::Relaxed)).min().unwrap())
    }
}

#[inline]
fn encode(base: u8) -> Option<u64> {
    match base {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Reverse complement of a forward `K`-mer, in the encoding of `kmer::CanonicalKmers`.
#[inline]
fn rc(x: u64) -> u64 {
    let mut y = !x;
    y = ((y >> 2) & 0x3333_3333_3333_3333) | ((y & 0x3333_3333_3333_3333) << 2);
    y = ((y >> 4) & 0x0F0F_0F0F_0F0F_0F0F) | ((y & 0x0F0F_0F0F_0F0F_0F0F) << 4);
    y.swap_bytes() >> (64 - 2 * K)
}

#[inline]
fn canonical(x: u64) -> u64 {
    x.min(rc(x))
}

/// Forward `K`-mers of a sequence, as (start, word); windows with a non-ACGT base are skipped.
fn forward_kmers(seq: &[u8]) -> Vec<(usize, u64)> {
    let (mut out, mut word, mut valid) = (Vec::new(), 0u64, 0usize);
    for (i, &b) in seq.iter().enumerate() {
        match encode(b) {
            Some(c) => {
                word = ((word << 2) | c) & MASK;
                valid += 1;
                if valid >= K {
                    out.push((i + 1 - K, word));
                }
            }
            None => valid = 0,
        }
    }
    out
}

fn decode(word: u64) -> Vec<u8> {
    (0..K)
        .map(|i| b"ACGT"[((word >> (2 * (K - 1 - i))) & 3) as usize])
        .collect()
}

#[inline]
fn successors(word: u64) -> [u64; 4] {
    std::array::from_fn(|b| ((word << 2) | b as u64) & MASK)
}

#[inline]
fn predecessors(word: u64) -> [u64; 4] {
    std::array::from_fn(|b| (word >> 2) | ((b as u64) << (2 * (K - 1))))
}

/// The frontier extension of one target, for `reconstruct.json`.
#[derive(Debug, Clone, Serialize, Default)]
pub struct Extension {
    pub reference_kmers: usize,
    /// Median count of the reference k-mers: the target's depth in the table.
    pub depth: f64,
    pub floor: f64,
    pub ceiling: f64,
    /// Reference k-mers at the target's depth.
    pub present_reference_kmers: usize,
    /// Places where the reference stops being present, each walked from.
    pub ends: usize,
    /// Ends whose walk found at least one new k-mer.
    pub ends_extended: usize,
    /// Steps of the walk back onto reference k-mers at depth.
    pub joins: usize,
    /// Steps where two next k-mers were at depth and comparable: both walked.
    pub branches: usize,
    /// Steps with no next k-mer at depth.
    pub dead_ends: usize,
    /// Ends whose walk stopped at `max_per_end`.
    pub capped_ends: usize,
    /// Walking stopped at the target's budget (`max_fold`).
    pub budget_reached: bool,
    pub new_kmers: usize,
    /// Walked pieces left out as tips (`MIN_TIP`).
    pub tips_removed: usize,
    pub unitigs: usize,
    pub unitig_bases: usize,
    /// None, or why the target was not walked.
    pub skipped: Option<&'static str>,
    #[serde(skip)]
    pub sequences: Vec<Vec<u8>>,
}

/// Fewest median count for a target to be walked: below it, depth cannot tell it from the
/// nucleus, and its reference is mostly absent from the individual.
const MIN_DEPTH: f64 = 10.0;

/// Reference k-mers in a row a walk takes before it counts as rejoined (about two reads). A
/// short stretch the reference holds elsewhere does not end it: ZM4's walk into 8.5 kb no
/// reference run has stopped 1,480 bp short at ~30 bp the reference holds, from both sides, and
/// those reads were never recruited.
const JOIN_RUN: usize = 300;

/// Shortest walked piece kept when one of its ends leads nowhere (about a read).
const MIN_TIP: usize = 150;

/// Walk the frontier of a target whose reference is `reference` (its walk and node sequences,
/// as baited) and holds `stored_bases` (for the budget).
pub fn extend(
    counts: &KmerCounts,
    reference: &[Vec<u8>],
    stored_bases: u64,
    p: &Params,
) -> Extension {
    let mut known: FxHashSet<u64> = FxHashSet::default();
    let words: Vec<Vec<(usize, u64)>> = reference.iter().map(|s| forward_kmers(s)).collect();
    for w in &words {
        known.extend(w.iter().map(|&(_, x)| canonical(x)));
    }
    let mut ext = Extension {
        reference_kmers: known.len(),
        ..Default::default()
    };
    if known.is_empty() {
        ext.skipped = Some("empty_reference");
        return ext;
    }
    let mut depths: Vec<u32> = known.par_iter().map(|&x| counts.get(x)).collect();
    depths.sort_unstable();
    ext.depth = f64::from(depths[depths.len() / 2]);
    let (floor, ceiling) = (p.floor * ext.depth, p.ceiling * ext.depth);
    (ext.floor, ext.ceiling) = (floor, ceiling);
    if ext.depth < MIN_DEPTH {
        ext.skipped = Some("too_shallow");
        return ext;
    }
    let at_depth = move |c: u32| f64::from(c) >= floor && f64::from(c) <= ceiling;
    ext.present_reference_kmers = depths.iter().filter(|&&c| at_depth(c)).count();
    // Ends: a present reference k-mer whose next one along the sequence (either strand) is not,
    // or that ends the sequence. Node sequences end where the graph goes on; their walk rejoins
    // the next node at once and adds nothing.
    let present = |x: u64| at_depth(counts.get(canonical(x)));
    let mut ends: Vec<u64> = words
        .par_iter()
        .flat_map_iter(|w| {
            let mut out = Vec::new();
            for strand in [false, true] {
                let seq: Vec<u64> = if strand {
                    w.iter().rev().map(|&(_, x)| rc(x)).collect()
                } else {
                    w.iter().map(|&(_, x)| x).collect()
                };
                for i in 0..seq.len() {
                    if present(seq[i]) && (i + 1 == seq.len() || !present(seq[i + 1])) {
                        out.push(seq[i]);
                    }
                }
            }
            out
        })
        .collect();
    ends.sort_unstable();
    ends.dedup();
    ext.ends = ends.len();
    let budget = (p.max_fold * stored_bases as f64) as usize;
    let mut new: FxHashSet<u64> = FxHashSet::default();
    'ends: for &end in &ends {
        // (k-mer, how many reference k-mers in a row the walk has just taken to reach it)
        let (mut stack, mut walked) = (vec![(end, 0usize)], 0usize);
        // reference k-mers this end's walk went through (`JOIN_RUN`), each once
        let mut through: FxHashSet<u64> = FxHashSet::default();
        while let Some((word, run)) = stack.pop() {
            let next: Vec<(u64, u32)> = successors(word)
                .into_iter()
                .map(|x| (x, counts.get(canonical(x))))
                .filter(|&(_, c)| at_depth(c))
                .collect();
            let Some(top) = next.iter().map(|n| n.1).max() else {
                ext.dead_ends += 1;
                continue;
            };
            let taken: Vec<u64> = next
                .iter()
                .filter(|n| f64::from(n.1) >= p.dominance * f64::from(top))
                .map(|n| n.0)
                .collect();
            if taken.len() > 1 {
                ext.branches += 1;
            }
            for x in taken {
                let key = canonical(x);
                if known.contains(&key) {
                    if run + 1 >= JOIN_RUN {
                        ext.joins += 1;
                    } else if through.insert(key) {
                        stack.push((x, run + 1));
                    }
                    continue;
                }
                if !new.insert(key) {
                    continue;
                }
                walked += 1;
                if new.len() >= budget {
                    ext.budget_reached = true;
                    ext.ends_extended += 1;
                    break 'ends;
                }
                if walked >= p.max_per_end {
                    ext.capped_ends += 1;
                    break;
                }
                stack.push((x, 0));
            }
            if walked >= p.max_per_end {
                break;
            }
        }
        ext.ends_extended += usize::from(walked > 0);
    }
    ext.new_kmers = new.len();
    // A piece that leads nowhere at one end and is shorter than a read is a sequencing error or
    // a collision in the table followed for a few steps, not sequence of the target.
    let anchored = |x: u64| known.contains(&canonical(x)) || new.contains(&canonical(x));
    let (kept, tips): (Vec<Vec<u8>>, Vec<Vec<u8>>) = unitigs(&new).into_iter().partition(|s| {
        let words = forward_kmers(s);
        let (first, last) = (words[0].1, words[words.len() - 1].1);
        let left = predecessors(first).into_iter().any(anchored);
        let right = successors(last).into_iter().any(anchored);
        (left && right) || s.len() >= MIN_TIP
    });
    ext.tips_removed = tips.len();
    ext.sequences = kept;
    ext.unitigs = ext.sequences.len();
    ext.unitig_bases = ext.sequences.iter().map(Vec::len).sum();
    ext
}

/// Maximal nonbranching chains of a set of canonical k-mers, as sequences.
fn unitigs(set: &FxHashSet<u64>) -> Vec<Vec<u8>> {
    let has = |x: u64| set.contains(&canonical(x));
    let out_of = |x: u64| successors(x).into_iter().filter(|&y| has(y)).collect::<Vec<_>>();
    let into = |x: u64| predecessors(x).into_iter().filter(|&y| has(y)).count();
    // the one next k-mer of `x` in the set that has `x` as its only predecessor
    let next = |x: u64| match out_of(x)[..] {
        [y] if into(y) == 1 => Some(y),
        _ => None,
    };
    let mut keys: Vec<u64> = set.iter().copied().collect();
    keys.sort_unstable();
    let mut used: FxHashSet<u64> = FxHashSet::default();
    let mut result = Vec::new();
    for key in keys {
        if used.contains(&key) {
            continue;
        }
        // back to the chain's first k-mer: forward along the reverse complement (a cycle stops
        // where it comes round)
        let (mut first, mut back) = (key, FxHashSet::default());
        back.insert(key);
        while let Some(y) = next(rc(first)).map(rc) {
            if used.contains(&canonical(y)) || !back.insert(canonical(y)) {
                break;
            }
            first = y;
        }
        used.insert(canonical(first));
        let mut seq = decode(first);
        let mut word = first;
        while let Some(y) = next(word) {
            if !used.insert(canonical(y)) {
                break;
            }
            seq.push(b"ACGT"[(y & 3) as usize]);
            word = y;
        }
        result.push(seq);
    }
    result
}

/// Write the walked sequences of a target as FASTA (diagnostic: they are bait, not output).
pub fn write_fasta(path: &Path, name: &str, ext: &Extension) -> Result<()> {
    let mut text = String::new();
    for (i, s) in ext.sequences.iter().enumerate() {
        text.push_str(&format!(">{name}.frontier.{} length={}\n", i + 1, s.len()));
        text.push_str(std::str::from_utf8(s)?);
        text.push('\n');
    }
    std::fs::write(path, text)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // high bits of a mixed counter: xorshift's low bits repeat 31-mers in a few kb
    fn random(n: usize, seed: u64) -> Vec<u8> {
        (0..n as u64)
            .map(|i| b"ACGT"[(mix64(seed << 32 | i) >> 62) as usize])
            .collect()
    }

    fn revcomp(s: &[u8]) -> Vec<u8> {
        crate::evidence::revcomp(s)
    }

    /// Reads of `genome` every `step` bases, both strands, `len` long.
    fn tile(table: &KmerCounts, genome: &[u8], len: usize, step: usize, copies: usize) {
        for _ in 0..copies {
            for a in (0..genome.len().saturating_sub(len)).step_by(step) {
                table.add_seq(&genome[a..a + len]);
                table.add_seq(&revcomp(&genome[a..a + len]));
            }
        }
    }

    #[test]
    fn words_match_the_canonical_kmers_of_the_bait() {
        let s = random(200, 7);
        let ours: Vec<u64> = forward_kmers(&s).iter().map(|&(_, x)| canonical(x)).collect();
        let theirs: Vec<u64> = CanonicalKmers::new(&s, K).map(|(_, x)| x).collect();
        assert_eq!(ours, theirs);
        let first = forward_kmers(&s)[0].1;
        assert_eq!(decode(first), s[..K].to_vec());
        assert_eq!(rc(first), forward_kmers(&revcomp(&s[..K]))[0].1);
    }

    #[test]
    fn counts_are_at_least_the_occurrences() {
        let table = KmerCounts::new(16);
        let s = random(400, 3);
        for _ in 0..7 {
            table.add_seq(&s);
        }
        table.add_seq(&revcomp(&s));
        for (_, x) in CanonicalKmers::new(&s, K) {
            assert_eq!(table.get(x), 8);
        }
        assert_eq!(table.get(canonical(forward_kmers(&random(40, 99))[0].1)), 0);
    }

    #[test]
    fn walks_across_sequence_no_reference_holds_and_rejoins() {
        let genome = random(4000, 11);
        let table = KmerCounts::new(22);
        tile(&table, &genome, 150, 3, 1);
        // the reference lacks 1,000-1,800 and holds another allele there
        let reference = [genome[..1000].to_vec(), random(300, 5), genome[1800..].to_vec()].concat();
        let ext = extend(&table, &[reference], 3300, &Params::default());
        assert!(ext.skipped.is_none());
        assert!(ext.joins >= 1, "{ext:?}");
        let missing = &genome[1000 - K + 1..1800 + K - 1];
        let found = ext
            .sequences
            .iter()
            .any(|s| s.windows(missing.len()).any(|w| w == missing) || {
                let r = revcomp(s);
                r.windows(missing.len()).any(|w| w == missing)
            });
        assert!(found, "{ext:?}");
        assert_eq!(ext.new_kmers, 800 + K - 1);
    }

    #[test]
    fn a_short_stretch_the_reference_holds_elsewhere_does_not_stop_the_walk() {
        let genome = random(4000, 13);
        let table = KmerCounts::new(22);
        tile(&table, &genome, 150, 3, 1);
        // the reference lacks 1,000-2,000 but holds 1,400-1,440 elsewhere (another run's copy)
        let reference = [
            genome[..1000].to_vec(),
            random(200, 6),
            genome[1400..1440].to_vec(),
            random(200, 8),
            genome[2000..].to_vec(),
        ]
        .concat();
        let ext = extend(&table, &[reference], 3500, &Params::default());
        // every k-mer between the two pieces but those of the 40 bp the reference holds
        assert_eq!(ext.new_kmers, 1000 + K - 1 - (40 - K + 1), "{ext:?}");
    }

    #[test]
    fn a_short_dead_end_off_the_walk_is_left_out() {
        let genome = random(4000, 17);
        let table = KmerCounts::new(22);
        tile(&table, &genome, 150, 3, 1);
        // reads as deep as the genome's that leave it at 1,200 for 40 bases and stop
        let spur = [genome[1080..1200].to_vec(), random(40, 9)].concat();
        for _ in 0..40 {
            table.add_seq(&spur);
            table.add_seq(&revcomp(&spur));
        }
        // a reference of the first 1,000 bases only: the one walk passes the branch at 1,200
        let ext = extend(&table, &[genome[..1000].to_vec()], 3400, &Params::default());
        assert_eq!(ext.branches, 1, "{ext:?}");
        assert!(ext.tips_removed >= 1, "{ext:?}");
        let tail = &spur[120..];
        assert!(ext.sequences.iter().all(|s| {
            !s.windows(20).any(|w| w == &tail[..20]) && !revcomp(s).windows(20).any(|w| w == &tail[..20])
        }));
        // the walk on from 1,000 stays, in two pieces either side of the branch
        let held: FxHashSet<u64> = ext
            .sequences
            .iter()
            .flat_map(|s| forward_kmers(s).into_iter().map(|(_, x)| canonical(x)))
            .collect();
        let missing = &genome[1000 - K + 1..1600 + K - 1];
        assert!(forward_kmers(missing).iter().all(|&(_, x)| held.contains(&canonical(x))));
    }

    #[test]
    fn shallow_nuclear_sequence_is_not_walked() {
        let organelle = random(3000, 21);
        let nucleus = [random(1500, 22), organelle[1200..1500].to_vec(), random(1500, 23)].concat();
        let table = KmerCounts::new(22);
        tile(&table, &organelle, 150, 3, 1); // ~100x
        tile(&table, &nucleus, 150, 30, 1); // ~10x, with a NUMT of the organelle
        let ext = extend(&table, &[organelle.clone()], 3000, &Params::default());
        assert_eq!(ext.new_kmers, 0, "{ext:?}");
        // and an organelle reference missing a piece still finds only that piece
        let partial = [organelle[..1000].to_vec(), organelle[1100..].to_vec()].concat();
        let ext = extend(&table, &[partial], 2900, &Params::default());
        assert_eq!(ext.new_kmers, 100 + K - 1, "{ext:?}");
    }

    #[test]
    fn budget_stops_the_walk() {
        let genome = random(4000, 31);
        let table = KmerCounts::new(22);
        tile(&table, &genome, 150, 3, 1);
        let ext = extend(&table, &[genome[..500].to_vec()], 50, &Params::default());
        assert!(ext.budget_reached);
        assert_eq!(ext.new_kmers, 100);
    }
}
