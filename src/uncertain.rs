//! Stretches of a short-read molecule whose length the reads do not fix: tandem arrays and
//! microsatellites whose copies a de Bruijn graph collapses, longer than a read with a unique
//! anchor on either side.
//!
//! ZM4's plastid circle lacks 266 bp of the held-out HiFi plastid at four AT-rich arrays of
//! 128-247 bp with their anchors (36 and 21 bp units, two and three copies fewer). No 150 bp
//! read spans one; read pairs across them, against pairs across nearby places, put the circle
//! 20-30 bp short at three of them. So such a stretch is marked rather than trusted: where it
//! is, how many reads span it against how many would at its assembled length, and how much
//! longer the pairs across it say it is (an estimate: inserts vary with the place by +-20 bp).

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::Result;
use rustc_hash::FxHashMap;
use serde_json::{json, Value};

use crate::evidence::revcomp;
use crate::kmer::CanonicalKmers;
use crate::recruit::for_each_rec;

/// Word length of repeats and anchors.
const W: usize = 20;
/// Two copies of a word this close are a local repeat (a tandem array or microsatellite).
const LOCAL: usize = 500;
/// Repeated words this close together are one stretch.
const MERGE_GAP: usize = 20;
/// Bases either side of a stretch where control pairs are taken (`run`).
const CONTROL_OFFSET: usize = 1500;
/// A stretch is marked when its reads span it at under this share of the expected rate. ZM4's
/// plastid: 7 of ~185 reads (0.04) at a collapsed 21 bp array, 0.13-0.2 at stretches its HiFi
/// assembly has at the same length (AT-rich, thinly read).
const SPANNING_SHARE: f64 = 0.1;

struct Stretch {
    start: usize,
    end: usize,
    /// The unique words on either side: the left one ends at `left_end`, the right starts at
    /// `right_start`.
    left_end: usize,
    right_start: usize,
}

fn name_of(id: &[u8]) -> Vec<u8> {
    let name = id.split(|&b| b == b' ' || b == b'\t').next().unwrap_or(id);
    name.strip_suffix(b"/1").or_else(|| name.strip_suffix(b"/2")).unwrap_or(name).to_vec()
}

fn median(mut v: Vec<i64>) -> Option<f64> {
    v.sort_unstable();
    (!v.is_empty()).then(|| v[v.len() / 2] as f64)
}

/// Local repeat stretches of a molecule, with their nearest unique anchors.
fn stretches(seq: &[u8]) -> Vec<Stretch> {
    let mut at: FxHashMap<u64, Vec<usize>> = FxHashMap::default();
    for (p, x) in CanonicalKmers::new(seq, W) {
        at.entry(x).or_default().push(p);
    }
    let n = seq.len().saturating_sub(W - 1);
    let mut unique = vec![false; n];
    let mut local = vec![false; n];
    for ps in at.values() {
        for (k, &p) in ps.iter().enumerate() {
            if p >= n {
                continue;
            }
            unique[p] = ps.len() == 1;
            let near = (k > 0 && p - ps[k - 1] <= LOCAL) || (k + 1 < ps.len() && ps[k + 1] - p <= LOCAL);
            local[p] = near;
        }
    }
    let mut out = Vec::new();
    let mut p = 0;
    while p < n {
        if !local[p] {
            p += 1;
            continue;
        }
        let start = p;
        let mut end = p + W;
        let mut q = p + 1;
        while q < n && q <= end + MERGE_GAP {
            if local[q] {
                end = q + W;
            }
            q += 1;
        }
        let left = (0..start.saturating_sub(W - 1)).rev().find(|&i| unique[i]);
        let right = (end.min(n)..n).find(|&i| unique[i]);
        if let (Some(l), Some(r)) = (left, right) {
            out.push(Stretch {
                start,
                end: end.min(seq.len()),
                left_end: l + W,
                right_start: r,
            });
        }
        p = q;
    }
    out
}

/// Pick the middle of the longest circular run of unique words. Ties are
/// resolved by sequence, never input coordinates. Only k-mer indexing wraps;
/// analysis uses one copy of the molecule and scans each read once.
fn circular_cut(seq: &[u8]) -> usize {
    let n = seq.len();
    if n < W { return 0; }
    let extended = [seq, &seq[..W - 1]].concat();
    let words: Vec<_> = CanonicalKmers::new(&extended, W).take(n).collect();
    if words.len() != n { return 0; }
    let mut counts = FxHashMap::default();
    for &(_, key) in &words { *counts.entry(key).or_insert(0usize) += 1; }
    let unique: Vec<_> = words.iter().map(|(_, key)| counts[key] == 1).collect();
    let Some(start) = unique.iter().position(|&u| !u) else {
        return (0..n).min_by_key(|&i| &extended[i..i+W]).unwrap_or(0);
    };
    let (mut run, mut best_len, mut best) = (0usize, 0usize, 0usize);
    for step in 1..=n {
        let i = (start + step) % n;
        if unique[i] { run += 1; } else if run > 0 {
            let mid = (i + n - run + run / 2 + W / 2) % n;
            if run > best_len || (run == best_len && extended[mid..mid+W] < extended[best..best+W]) {
                best_len = run; best = mid;
            }
            run = 0;
        }
    }
    best
}

/// Mark the stretches of `molecules` (name, sequence, circular; those of at least `min_len`)
/// whose length the reads (`reads`, mates sharing their name's first word) do not fix. Writes
/// `uncertain_regions.tsv` into `out`.
pub fn run(molecules: &[(String, Vec<u8>, bool)], reads: &[PathBuf], min_len: usize, out: &Path) -> Result<Value> {
    let mut tsv = String::from(
        "molecule\tstart\tend\tlength\tanchor_distance\tspanning_reads\texpected_spanning\tspanning_pairs\tpair_longer_bp\tpair_se\tmarked\treason\twraps_origin\tanalysis_rotation\n",
    );
    let mut marked_total = 0usize;
    let mut listed = Vec::new();
    for (name, seq, circular) in molecules.iter().filter(|m| m.1.len() >= min_len) {
        let mut seq = seq.to_ascii_uppercase();
        let rotation = if *circular { circular_cut(&seq) } else { 0 };
        seq.rotate_left(rotation);
        let st = stretches(&seq);
        if st.is_empty() {
            continue;
        }
        // reads: anchors of the stretches, and placement by unique words for the pairs
        let mut anchors: FxHashMap<Vec<u8>, Vec<(usize, bool)>> = FxHashMap::default();
        for (i, s) in st.iter().enumerate() {
            anchors.entry(seq[s.left_end - W..s.left_end].to_vec()).or_default().push((i, true));
            anchors.entry(seq[s.right_start..s.right_start + W].to_vec()).or_default().push((i, false));
        }
        let mut place: FxHashMap<Vec<u8>, i64> = FxHashMap::default();
        let mut counts: FxHashMap<&[u8], u32> = FxHashMap::default();
        for x in seq.windows(W) {
            *counts.entry(x).or_insert(0) += 1;
        }
        let rc_seq = revcomp(&seq);
        for x in rc_seq.windows(W) {
            *counts.entry(x).or_insert(0) += 1;
        }
        for (p, x) in seq.windows(W).enumerate() {
            if counts[x] == 1 {
                place.insert(x.to_vec(), p as i64);
            }
        }
        let mut spans: Vec<Vec<i64>> = vec![Vec::new(); st.len()];
        let mut placed: FxHashMap<Vec<u8>, Vec<(bool, i64, i64)>> = FxHashMap::default();
        let (mut reads_n, mut read_bases) = (0u64, 0u64);
        for path in reads {
            for_each_rec(path, |r| {
                let s = r.seq.to_ascii_uppercase();
                reads_n += 1;
                read_bases += s.len() as u64;
                let rs = revcomp(&s);
                for strand in [&s, &rs] {
                    let mut left: BTreeMap<usize, usize> = BTreeMap::new();
                    for (p, x) in strand.windows(W).enumerate() {
                        for &(i, is_left) in anchors.get(x).into_iter().flatten() {
                            if is_left {
                                left.entry(i).or_insert(p + W);
                            } else if let Some(&e) = left.get(&i) {
                                if p >= e { spans[i].push((p - e) as i64); }
                            }
                        }
                    }
                }
                let n = s.len() as i64;
                if let Some((i, p)) = s.windows(W).enumerate().find_map(|(i, x)| place.get(x).map(|&p| (i, p))) {
                    placed.entry(name_of(&r.id)).or_default().push((true, p - i as i64, p - i as i64 + n));
                } else if let Some((i, p)) = rs.windows(W).enumerate().find_map(|(i, x)| place.get(x).map(|&p| (i, p))) {
                    placed.entry(name_of(&r.id)).or_default().push((false, p - i as i64, p - i as i64 + n));
                }
            })?;
        }
        // proper pairs: (leftmost, rightmost) on the molecule
        let pairs: Vec<(i64, i64)> = placed
            .values()
            .filter_map(|v| {
                let f = v.iter().find(|x| x.0)?;
                let r = v.iter().find(|x| !x.0)?;
                (r.2 > f.1 && r.2 - f.1 < 2000).then_some((f.1, r.2))
            })
            .collect();
        let read_len = if reads_n > 0 { (read_bases / reads_n) as f64 } else { 0.0 };
        let starts_per_base = reads_n as f64 / seq.len() as f64 / 2.0; // both strands of each read were scanned
        let across = |a: i64, b: i64| -> Vec<i64> {
            pairs.iter().filter(|&&(s, e)| s < a && e > b).map(|&(s, e)| e - s).collect()
        };
        for (i, s) in st.iter().enumerate() {
            let d = (s.right_start - s.left_end) as f64;
            // a read holds both anchors when it starts within read_len - d - 2W bases
            let room = (read_len - d - 2.0 * W as f64 + 1.0).max(0.0);
            let expected = 2.0 * starts_per_base * room;
            let observed = spans[i].len();
            let over = across(s.left_end as i64, s.right_start as i64);
            let (lo, hi) = (s.left_end.saturating_sub(CONTROL_OFFSET), s.right_start + CONTROL_OFFSET);
            let mut control = Vec::new();
            for c in [lo, hi.min(seq.len().saturating_sub(1))] {
                control.extend(across(c as i64 - 10, c as i64 + 10));
            }
            let (shift, se) = match (median(over.clone()), median(control.clone())) {
                (Some(m), Some(c)) if over.len() >= 10 => {
                    let mean = over.iter().sum::<i64>() as f64 / over.len() as f64;
                    let sd = (over.iter().map(|&x| (x as f64 - mean).powi(2)).sum::<f64>() / over.len() as f64).sqrt();
                    (Some(c - m), Some(1.2533 * sd / (over.len() as f64).sqrt()))
                }
                _ => (None, None),
            };
            let discordant = spans[i].iter().filter(|&&length| length != d as i64).count();
            let reason = if observed >= 3 && discordant * 5 >= observed * 4 {
                Some("spanning reads disagree with assembled repeat length")
            // Do not call a repeat resolved when both its assembled length and a
            // different length have substantial spanning-read support. This is
            // a conservative warning, not a heteroplasmy or allele call.
            } else if discordant >= 3 && discordant * 5 >= observed && observed - discordant >= 3 {
                Some("spanning reads support multiple repeat lengths")
            } else if room <= 0.0 {
                Some("no read can span it")
            } else if (observed as f64) < SPANNING_SHARE * expected {
                Some("far fewer reads span it than its length allows")
            } else {
                // the pairs' estimate is reported, not a reason: inserts across a place vary
                // by +-20 bp with the place alone (ZM4's plastid, two controls: 10 and 18 bp)
                None
            };
            let start = (s.start + rotation) % seq.len();
            let end = if *circular { (s.end + rotation) % seq.len() } else { s.end };
            let wraps = *circular && start > end;
            if reason.is_some() {
                marked_total += 1;
                listed.push(json!({"molecule": name, "start": start, "end": end, "wraps_origin": wraps, "analysis_rotation": rotation, "anchor_distance": d,
                    "spanning_reads": observed, "expected_spanning": expected, "pair_longer_bp": shift,
                    "pair_se": se, "reason": reason}));
            }
            writeln!(
                tsv,
                "{name}\t{}\t{}\t{}\t{d}\t{observed}\t{expected:.1}\t{}\t{}\t{}\t{}\t{}\t{wraps}\t{rotation}",
                start,
                end,
                s.end - s.start,
                over.len(),
                shift.map_or("-".to_string(), |x| format!("{x:.0}")),
                se.map_or("-".to_string(), |x| format!("{x:.1}")),
                if reason.is_some() { "yes" } else { "no" },
                reason.unwrap_or("-"),
            )?;
        }
    }
    fs::write(out.join("uncertain_regions.tsv"), tsv)?;
    Ok(json!({"marked": marked_total, "regions": listed, "table": out.join("uncertain_regions.tsv"),
        "coordinates": "zero-based half-open; wraps_origin=true means [start,length) plus [0,end)",
        "basis": "local repeats (20-mers repeated within 500 bp) between unique anchors: marked when no read can span them, far fewer do than expected, or spanning reads support different lengths; pair_longer_bp (pairs across against pairs 1.5 kb away) is an estimate only; lengths not corrected"}))
}

#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn circular_repeat_detection_is_rotation_invariant() {
        let unit = random(21, 13);
        let original = [random(1000, 1), unit.repeat(8), random(1300, 2)].concat();
        let mut canonical = original.clone();
        canonical.rotate_left(circular_cut(&original));
        let expected = stretches(&canonical);
        assert!(!expected.is_empty());
        for offset in [1, 990, 1050, 1100, original.len()-1] {
            let mut rotated = original.clone(); rotated.rotate_left(offset);
            let cut = circular_cut(&rotated); rotated.rotate_left(cut);
            assert_eq!(rotated, canonical, "offset {offset}");
            let observed = stretches(&rotated);
            assert_eq!(observed.len(), expected.len());
            for (a,b) in observed.iter().zip(&expected) {
                assert_eq!((a.start,a.end,a.left_end,a.right_start),(b.start,b.end,b.left_end,b.right_start));
            }
        }
    }

    #[test]
    fn mixed_spanning_lengths_are_not_called_certain() {
        let unit = random(21, 33);
        let molecule = [random(700, 11), unit.repeat(3), random(700, 12)].concat();
        let st = stretches(&molecule);
        assert_eq!(st.len(), 1);
        let s = &st[0];
        let exact = molecule[s.left_end-W..s.right_start+W].to_vec();
        let variant = [molecule[s.left_end-W..s.start+21].to_vec(), unit.clone(),
                       molecule[s.start+21..s.right_start+W].to_vec()].concat();
        let dir = std::env::temp_dir().join(format!("ovasm-mixed-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let mut fq = String::new();
        for (i, read) in [exact, variant].iter().cycle().take(12).enumerate() {
            writeln!(fq, "@r{i}\n{}\n+\n{}", String::from_utf8_lossy(read), "I".repeat(read.len())).unwrap();
        }
        let reads = dir.join("reads.fq");
        fs::write(&reads, fq).unwrap();
        let result = run(&[("m".to_string(), molecule, false)], &[reads], 1000, &dir).unwrap();
        assert_eq!(result["marked"], 1, "{result}");
        assert_eq!(result["regions"][0]["reason"], "spanning reads support multiple repeat lengths");
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn mate_suffixes_share_fragment_identity() {
        assert_eq!(name_of(b"fragment/1 comment"), name_of(b"fragment/2 comment"));
        assert_ne!(name_of(b"fragment/3"), name_of(b"fragment/1"));
    }

    #[test]
    fn a_collapsed_tandem_array_is_marked_with_its_missing_length() {
        // the DNA holds seven copies of a 21 bp unit; the molecule three
        let unit = random(21, 3);
        let (left, right) = (random(3000, 1), random(3000, 2));
        let truth = [left.clone(), unit.repeat(7), right.clone()].concat();
        let molecule = [left, unit.repeat(3), right].concat();
        let dir = std::env::temp_dir().join(format!("ovasm-uncertain-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        // pairs of 100 bp mates of 330 bp fragments every 2 bases
        let mut fq = String::new();
        for (f, a) in (0..truth.len() - 330).step_by(2).enumerate() {
            let m1 = &truth[a..a + 100];
            let m2 = revcomp(&truth[a + 230..a + 330]);
            for (m, s) in [(1, m1.to_vec()), (2, m2)] {
                writeln!(fq, "@f{f} {m}\n{}\n+\n{}", String::from_utf8(s).unwrap(), "I".repeat(100)).unwrap();
            }
        }
        let reads = dir.join("r.fq");
        fs::write(&reads, fq).unwrap();
        let s = run(&[("m".to_string(), molecule, false)], &[reads], 1000, &dir).unwrap();
        assert_eq!(s["marked"], 1, "{s}");
        let shift = s["regions"][0]["pair_longer_bp"].as_f64().unwrap();
        assert!((shift - 84.0).abs() <= 12.0, "{s}");
        fs::remove_dir_all(dir).unwrap();
    }
}
