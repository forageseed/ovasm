//! Circular plastid molecules in the conventional layout: LSC, IRb, SSC, IRa, starting at the
//! first base of the large single-copy region (the IRa/LSC junction), on the strand that holds
//! rbcL forward and psbA reverse, with the small single-copy region turned so that ndhF is on
//! the reverse strand (the orientation of most GenBank plastomes, Arabidopsis NC_000932 and
//! tobacco among them).
//!
//! The assembler starts a circle wherever its walk starts and reads it on either strand, and the
//! SSC exists in both orientations in vivo (flip-flop recombination across the inverted repeat),
//! so two runs of one sample, or two samples of one species, came out rotated, reversed or with
//! the SSC turned: unlike files for the same plastome. The layout is decided from the sequence
//! alone, so it is the same whatever the graph and wherever the walk began, and applying it twice
//! changes nothing.
//!
//! The inverted repeat is found by k-mers: a position is in it when the reverse complement of its
//! k-mer occurs on the molecule. The two longest such runs (gaps of up to `MAX_GAP` bases closed,
//! for the few differences between the copies) are the copies when each is at least `MIN_IR`
//! bases and the shorter is at least `IR_RATIO` of the longer; of the two single-copy regions
//! between them the longer is the LSC. Genes come from the embedded protein panel scanned in six
//! frames (`panel::gene_strands`). A molecule without an inverted repeat (the IR-lacking legumes,
//! conifers) is only oriented, not rotated; one without rbcL and psbA keeps its strand; one
//! without ndhF (it is lost in some lineages) keeps its SSC orientation. What was done is written
//! into the FASTA header (`layout=`, `rotated_by=`, `reverse_complemented=`, `ssc_turned=`): the
//! header's `path=` is the walk before these changes.

use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use rustc_hash::FxHashMap;
use serde_json::{json, Value};

use crate::evidence::revcomp;
use crate::panel::{gene_strands, Panel, PanelParams};

/// k-mer length of the inverted-repeat search.
const K: usize = 31;
/// Shortest inverted-repeat copy (plastid IRs are 10-30 kb; IRs of a few kb remain in some
/// lineages; shorter inverted repeats elsewhere on a plastome are a few hundred bases at most).
pub const MIN_IR: usize = 1000;
/// Gap (in k-mer starts) closed inside one IR copy, for the SNPs and small indels between copies.
const MAX_GAP: usize = 300;
/// The shorter IR copy is at least this share of the longer one.
const IR_RATIO: f64 = 0.9;

fn base_code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// (forward, reverse-complement) k-mer starting at each position of the circle; None where the
/// k-mer holds a base other than ACGT.
fn circle_kmers(seq: &[u8]) -> Vec<Option<(u64, u64)>> {
    let n = seq.len();
    let mut out = vec![None; n];
    if n < K {
        return out;
    }
    let mask = (1u64 << (2 * K)) - 1;
    let (mut f, mut r, mut valid) = (0u64, 0u64, 0usize);
    for i in 0..n + K - 1 {
        match base_code(seq[i % n]) {
            Some(c) => {
                f = ((f << 2) | c) & mask;
                r = (r >> 2) | ((3 - c) << (2 * (K - 1)));
                valid += 1;
            }
            None => valid = 0,
        }
        if i + 1 >= K && valid >= K {
            out[i + 1 - K] = Some((f, r));
        }
    }
    out
}

/// The two inverted-repeat copies as (start, length) on the circle, in circle order, or None.
fn inverted_repeat(seq: &[u8]) -> Option<[(usize, usize); 2]> {
    let n = seq.len();
    let km = circle_kmers(seq);
    let mut present: FxHashMap<u64, u32> = FxHashMap::default();
    for &(f, _) in km.iter().flatten() {
        *present.entry(f).or_default() += 1;
    }
    let ir: Vec<bool> = km
        .iter()
        .map(|x| matches!(x, Some((f, r)) if f != r && present.contains_key(r)))
        .collect();
    // runs around the circle, unrolled from a position outside every run
    let start = (0..n).find(|&i| !ir[i])?;
    let mut runs: Vec<(usize, usize)> = Vec::new(); // first and last k-mer start, unrolled
    for j in 1..=n {
        let p = start + j;
        if !ir[p % n] {
            continue;
        }
        match runs.last_mut() {
            Some(r) if p - r.1 <= MAX_GAP => r.1 = p,
            _ => runs.push((p, p)),
        }
    }
    let mut spans: Vec<(usize, usize)> = runs.iter().map(|&(a, b)| (a, b + K - a)).collect();
    spans.sort_by(|x, y| y.1.cmp(&x.1).then(x.0.cmp(&y.0)));
    if spans.len() < 2 {
        return None;
    }
    let (mut a, mut b) = (spans[0], spans[1]);
    if b.1 < MIN_IR || (b.1 as f64) < IR_RATIO * a.1 as f64 {
        return None;
    }
    if b.0 < a.0 {
        std::mem::swap(&mut a, &mut b);
    }
    // the copies may not overlap (a palindrome is no inverted repeat with regions between)
    if a.0 + a.1 >= b.0 || b.0 + b.1 >= a.0 + n {
        return None;
    }
    Some([(a.0 % n, a.1), (b.0 % n, b.1)])
}

/// Region lengths of a molecule in the layout: LSC, IRb, SSC, IRa.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Regions {
    lsc: usize,
    irb: usize,
    ssc: usize,
    ira: usize,
}

/// Rotate to the first base of the LSC (the longer single-copy region): LSC, IR, SSC, IR.
fn rotate_to_lsc(seq: &[u8], ir: [(usize, usize); 2]) -> (Vec<u8>, usize, Regions) {
    let n = seq.len();
    let [(a, la), (b, lb)] = ir;
    let b = if b < a { b + n } else { b };
    let s1 = b - (a + la); // after copy a, before copy b
    let s2 = a + n - (b + lb); // after copy b, before copy a
    let (start, r) = if s1 >= s2 {
        ((a + la) % n, Regions { lsc: s1, irb: lb, ssc: s2, ira: la })
    } else {
        ((b + lb) % n, Regions { lsc: s2, irb: la, ssc: s1, ira: lb })
    };
    let mut out = seq.to_vec();
    out.rotate_left(start);
    (out, start, r)
}

/// The strand a gene lies on within [from, to): the strand of most of its hits' bases there.
fn strand_of(hits: &[(String, usize, usize, bool)], gene: &str, from: usize, to: usize) -> Option<bool> {
    let (mut fwd, mut rev) = (0usize, 0usize);
    for (name, s, e, forward) in hits {
        let mid = (s + e) / 2;
        if name.eq_ignore_ascii_case(gene) && mid >= from && mid < to {
            if *forward {
                fwd += e - s
            } else {
                rev += e - s
            }
        }
    }
    match fwd.cmp(&rev) {
        std::cmp::Ordering::Greater => Some(true),
        std::cmp::Ordering::Less => Some(false),
        std::cmp::Ordering::Equal => None,
    }
}

fn panel_hits(seq: &[u8], panel: &Panel) -> Vec<(String, usize, usize, bool)> {
    let p = PanelParams::default();
    gene_strands(seq, panel, &p)
        .into_iter()
        .map(|(g, s, e, f)| (panel.gene_name(g).to_string(), s, e, f))
        .collect()
}

/// Whether the strand is the conventional one: rbcL forward and psbA reverse within [from, to).
/// None when neither gene says.
fn conventional_strand(hits: &[(String, usize, usize, bool)], from: usize, to: usize) -> Option<bool> {
    let votes: Vec<bool> = [("rbcl", true), ("psba", false)]
        .iter()
        .filter_map(|&(g, want)| strand_of(hits, g, from, to).map(|s| s == want))
        .collect();
    match (votes.iter().filter(|&&v| v).count(), votes.iter().filter(|&&v| !v).count()) {
        (0, 0) => None,
        (y, n) if y > n => Some(true),
        (y, n) if n > y => Some(false),
        _ => None, // the genes disagree (a rearranged plastome): leave the strand
    }
}

/// The conventional layout of one circular plastid molecule and what was done to reach it.
pub fn layout(seq: &[u8], panel: &Panel) -> (Vec<u8>, Value) {
    let n = seq.len();
    let Some(ir) = inverted_repeat(seq) else {
        // no inverted repeat: orient only
        let hits = panel_hits(seq, panel);
        let flip = conventional_strand(&hits, 0, n) == Some(false);
        let out = if flip { revcomp(seq) } else { seq.to_vec() };
        return (out, json!({"layout": "no_inverted_repeat", "rotated_by": 0, "reverse_complemented": flip, "ssc_turned": false}));
    };
    // `rotated`: where the LSC began on the input, on the input's strand
    let (mut s, rotated, mut r) = rotate_to_lsc(seq, ir);
    let mut hits = panel_hits(&s, panel);
    let flip = conventional_strand(&hits, 0, r.lsc) == Some(false);
    if flip {
        // rc(LSC IRb SSC IRa) = rc(IRa) rc(SSC) rc(IRb) rc(LSC): rotate to rc(LSC)
        s = revcomp(&s);
        s.rotate_left(n - r.lsc);
        r = Regions { lsc: r.lsc, irb: r.ira, ssc: r.ssc, ira: r.irb };
        hits = panel_hits(&s, panel);
    }
    let (ssc_from, ssc_to) = (r.lsc + r.irb, r.lsc + r.irb + r.ssc);
    let turn = strand_of(&hits, "ndhf", ssc_from, ssc_to) == Some(true);
    if turn {
        let t = revcomp(&s[ssc_from..ssc_to]);
        s[ssc_from..ssc_to].copy_from_slice(&t);
    }
    (
        s,
        json!({
            "layout": "LSC-IRb-SSC-IRa",
            "lsc": r.lsc, "irb": r.irb, "ssc": r.ssc, "ira": r.ira,
            "rotated_by": rotated,
            "reverse_complemented": flip,
            "ssc_turned": turn,
        }),
    )
}

/// Rewrite the circular records of a plastid FASTA in the conventional layout, in place; linear
/// records are left as they are. Returns one report per record.
pub fn rewrite_fasta(path: &Path, panel: &Panel) -> Result<Vec<Value>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            records.push((h.to_string(), Vec::new()));
        } else if let Some(r) = records.last_mut() {
            r.1.extend_from_slice(line.trim().as_bytes());
        }
    }
    let mut out = String::with_capacity(text.len() + 256);
    let mut reports = Vec::new();
    for (header, seq) in records {
        let circular = header.split_whitespace().any(|t| t == "circular=true");
        let mut fields: Vec<String> = header
            .split_whitespace()
            .filter(|t| !["layout=", "rotated_by=", "reverse_complemented=", "ssc_turned="].iter().any(|p| t.starts_with(p)))
            .map(str::to_string)
            .collect();
        let seq = if circular {
            let (s, rep) = layout(&seq, panel);
            fields.push(format!("layout={}", rep["layout"].as_str().unwrap_or("")));
            fields.push(format!("rotated_by={}", rep["rotated_by"]));
            fields.push(format!("reverse_complemented={}", rep["reverse_complemented"]));
            fields.push(format!("ssc_turned={}", rep["ssc_turned"]));
            let mut rep = rep;
            rep["name"] = json!(fields.first().cloned().unwrap_or_default());
            reports.push(rep);
            s
        } else {
            seq
        };
        out.push('>');
        out.push_str(&fields.join(" "));
        out.push('\n');
        for chunk in seq.chunks(80) {
            out.push_str(std::str::from_utf8(chunk).unwrap_or(""));
            out.push('\n');
        }
    }
    fs::write(path, out).with_context(|| format!("writing {}", path.display()))?;
    Ok(reports)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A deterministic pseudo-random sequence (no long repeats by chance).
    fn random(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize]
            })
            .collect()
    }

    fn plastome(lsc: &[u8], ir: &[u8], ssc: &[u8]) -> Vec<u8> {
        [lsc, ir, ssc, &revcomp(ir)].concat()
    }

    #[test]
    fn rotation_strand_and_ssc_turn_come_back_to_one_layout() {
        let (lsc, ir, ssc) = (random(8000, 1), random(2500, 2), random(1500, 3));
        let canon = plastome(&lsc, &ir, &ssc);
        let panel = Panel::embedded(6, &[]).unwrap();
        let (want, rep) = layout(&canon, &panel);
        assert_eq!(rep["layout"], "LSC-IRb-SSC-IRa");
        // the repeat is the longest exact inverted repeat: bases at the single-copy ends that
        // happen to pair with each other join it (here a few at the SSC), as in real plastomes
        let (l, s) = (rep["lsc"].as_u64().unwrap(), rep["ssc"].as_u64().unwrap());
        assert!((7980..=8000).contains(&l) && (1480..=1500).contains(&s), "{rep}");
        assert_eq!(rep["lsc"].as_u64().unwrap() + rep["irb"].as_u64().unwrap() + s + rep["ira"].as_u64().unwrap(), 14_500);
        // without genes the strand stays; any rotation of either strand gives one of two sequences
        let rc_want = {
            let mut s = revcomp(&want);
            s.rotate_left(want.len() - l as usize);
            s
        };
        for start in [0, 123, 8000, 9100, 12_400] {
            let mut s = canon.clone();
            s.rotate_left(start);
            assert_eq!(layout(&s, &panel).0, want, "rotation {start}");
            let (r, _) = layout(&revcomp(&s), &panel);
            assert_eq!(r, rc_want, "reverse complement, rotation {start}");
        }
        // applying it twice changes nothing
        assert_eq!(layout(&want, &panel).0, want);
    }

    #[test]
    fn a_molecule_without_an_inverted_repeat_is_left_where_it_starts() {
        let s = random(12_000, 9);
        let panel = Panel::embedded(6, &[]).unwrap();
        let (out, rep) = layout(&s, &panel);
        assert_eq!(rep["layout"], "no_inverted_repeat");
        assert_eq!(out, s);
    }

    #[test]
    fn copies_with_a_few_differences_are_one_repeat() {
        let (lsc, ir, ssc) = (random(8000, 4), random(3000, 5), random(1500, 6));
        let mut ira = revcomp(&ir);
        for i in [400, 1100, 2000] {
            ira[i] = if ira[i] == b'A' { b'C' } else { b'A' };
        }
        let s = [lsc.as_slice(), &ir, &ssc, &ira].concat();
        let rep = layout(&s, &Panel::embedded(6, &[]).unwrap()).1;
        assert_eq!(rep["layout"], "LSC-IRb-SSC-IRa");
        assert_eq!(rep["lsc"].as_u64(), Some(8000));
    }
}
