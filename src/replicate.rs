//! Downsampled replicates of an assembly and the junctions they reproduce (`ovasm run --mode
//! standard`).
//!
//! A junction is a link of the assembly graph with a branching end: where a unitig can be left
//! in more than one way, or entered from more than one. The link that closes a circular unitig
//! and any other link between two unbranched ends is not a junction (where a circle is cut is
//! arbitrary).
//!
//! **Junction signature.** For a link `a -> b` with an overlap of `ov` bases, take the junction
//! span: the `FLANK` bases of `a` before the overlap, the overlap itself, and the `FLANK` bases
//! of `b` after it (a flank is cut short at the end of its segment), in upper case. The
//! signature is a 64-bit hash (FNV-1a, then the splitmix64 finaliser) of whichever of the span
//! and its reverse complement sorts first. The reverse complement of the span of `a -> b` is the
//! span of `b' -> a'`, so the signature does not depend on the strand a graph happens to spell,
//! nor on segment names, and two graphs give the same signature for a junction they both hold at
//! the same position. (With overlaps that are not exact, which `ovasm assemble` never writes,
//! the two strands could differ inside the overlap.)
//!
//! **Reproduced.** A replicate reproduces a junction of the representative graph when it holds
//! a junction with the same signature (`J`), or when every canonical `KMER`-mer of the span is
//! spelled by the replicate's graph, inside a segment or across a link (`C`: the same sequence
//! adjacency, without a branch there, as when the replicate's k spans a repeat the
//! representative graph collapses, or a boundary sits a few bases away). Otherwise it is absent
//! (`-`). The reproduction rate is the share of replicates with `J` or `C`.
//!
//! **Adaptive stop.** Replicates are assembled at one k, so their junction sets (the sets of
//! signatures) are comparable. After the requested number, replicates stop as soon as the set
//! was the same in the last three (unchanged over two consecutive replicates), or at the
//! maximum.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use rustc_hash::FxHashSet;

use crate::evidence::{flip, revcomp, Graph, Oriented};
use crate::kmer::CanonicalKmers;
use crate::recruit::{read_batches, unit_hash};

/// Bases taken from each side of a junction for its signature.
pub const FLANK: usize = 100;
/// k-mer length of the "spelled by the graph" test.
pub const KMER: usize = 31;
/// Most replicates of a run.
pub const MAX_REPLICATES: usize = 5;
/// Fewest replicates before the stopping rule can apply (three equal sets).
pub const MIN_REPLICATES_TO_STOP: usize = 3;

/// The junction span of the link `a -> b` (see the module documentation).
pub fn junction_span(g: &Graph, a: Oriented, b: Oriented, ov: usize, flank: usize) -> Vec<u8> {
    let (sa, sb) = (&g.seqs[a.0], &g.seqs[b.0]);
    let ov_a = ov.min(sa.len());
    let ov_b = ov.min(sb.len());
    // `a` read in its orientation, up to and including the overlap
    let left_len = flank.min(sa.len() - ov_a) + ov_a;
    let mut span = if a.1 {
        sa[sa.len() - left_len..].to_vec()
    } else {
        revcomp(&sa[..left_len])
    };
    // `b` read in its orientation, after the overlap
    let right_len = flank.min(sb.len() - ov_b);
    if b.1 {
        span.extend_from_slice(&sb[ov_b..ov_b + right_len]);
    } else {
        span.extend(revcomp(&sb[sb.len() - ov_b - right_len..sb.len() - ov_b]));
    }
    span.make_ascii_uppercase();
    span
}

/// Orientation-independent 64-bit signature of a junction span.
pub fn signature(span: &[u8]) -> u64 {
    let rc = revcomp(span);
    let canonical = if span <= rc.as_slice() { span } else { &rc };
    let mut h: u64 = 0xCBF2_9CE4_8422_2325; // FNV-1a
    for &b in canonical {
        h ^= b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01B3);
    }
    // splitmix64 finaliser
    h = (h ^ (h >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    h = (h ^ (h >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    h ^ (h >> 31)
}

pub fn signature_hex(sig: u64) -> String {
    format!("{sig:016x}")
}

/// A link of a graph, one per adjacency (a link and its reverse complement are one).
#[derive(Clone, Debug)]
pub struct Junction {
    /// Position of the link in the graph's link list (and so in an evidence report's).
    pub index: usize,
    pub from: Oriented,
    pub to: Oriented,
    pub overlap: usize,
    /// One of its ends has another link: a real junction (see the module documentation).
    pub branching: bool,
    pub span: Vec<u8>,
    pub signature: u64,
}

fn canonical_link(a: Oriented, b: Oriented) -> (Oriented, Oriented) {
    let rc = (flip(b), flip(a));
    if (a, b) <= rc {
        (a, b)
    } else {
        rc
    }
}

/// Every link of the graph once, in file order, with its span, signature and whether it branches.
pub fn links(g: &Graph) -> Vec<Junction> {
    let mut seen: BTreeSet<(Oriented, Oriented)> = BTreeSet::new();
    let mut unique: Vec<(usize, Oriented, Oriented, usize)> = Vec::new();
    for (i, (&(a, b), &ov)) in g.links.iter().zip(&g.overlaps).enumerate() {
        if seen.insert(canonical_link(a, b)) {
            unique.push((i, a, b, ov));
        }
    }
    // ways out of each oriented end: a -> b leaves a, and (read backwards) leaves b'
    let mut out: BTreeMap<Oriented, usize> = BTreeMap::new();
    for &(_, a, b, _) in &unique {
        *out.entry(a).or_default() += 1;
        if (flip(b), flip(a)) != (a, b) {
            *out.entry(flip(b)).or_default() += 1;
        }
    }
    unique
        .into_iter()
        .map(|(index, a, b, ov)| {
            let span = junction_span(g, a, b, ov, FLANK);
            Junction {
                index,
                from: a,
                to: b,
                overlap: ov,
                branching: out[&a] > 1 || out[&flip(b)] > 1,
                signature: signature(&span),
                span,
            }
        })
        .collect()
}

/// The junction set of a graph: signatures of its branching links.
pub fn junction_set(links: &[Junction]) -> BTreeSet<u64> {
    links
        .iter()
        .filter(|j| j.branching)
        .map(|j| j.signature)
        .collect()
}

/// Canonical `KMER`-mers a graph spells: those of every segment and of every link's span.
pub fn spelled_kmers(g: &Graph, links: &[Junction]) -> FxHashSet<u64> {
    let mut set = FxHashSet::default();
    for s in &g.seqs {
        set.extend(CanonicalKmers::new(s, KMER).map(|(_, k)| k));
    }
    for j in links {
        set.extend(CanonicalKmers::new(&j.span, KMER).map(|(_, k)| k));
    }
    set
}

/// Whether every `KMER`-mer of a span is in `kmers`; `None` when the span holds none (shorter
/// than `KMER`, or broken by non-ACGT bases), so nothing can be said.
pub fn spelled(span: &[u8], kmers: &FxHashSet<u64>) -> Option<bool> {
    let mut n = 0usize;
    let mut all = true;
    for (_, k) in CanonicalKmers::new(span, KMER) {
        n += 1;
        all &= kmers.contains(&k);
    }
    // a span with a non-ACGT base yields fewer k-mers than its length allows: not judged
    (n > 0 && n == span.len() + 1 - KMER).then_some(all)
}

/// How one replicate holds a junction of the representative graph.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Reproduced {
    /// A junction with the same signature.
    Junction,
    /// The same sequence adjacency, spelled without a branch at that place.
    Contiguous,
    Absent,
}

impl Reproduced {
    pub fn letter(self) -> char {
        match self {
            Reproduced::Junction => 'J',
            Reproduced::Contiguous => 'C',
            Reproduced::Absent => '-',
        }
    }

    pub fn present(self) -> bool {
        self != Reproduced::Absent
    }
}

/// What a replicate's graph says about a junction span of another graph.
pub fn reproduced(
    sig: u64,
    span: &[u8],
    replicate_set: &BTreeSet<u64>,
    replicate_kmers: &FxHashSet<u64>,
) -> Reproduced {
    if replicate_set.contains(&sig) {
        Reproduced::Junction
    } else if spelled(span, replicate_kmers) == Some(true) {
        Reproduced::Contiguous
    } else {
        Reproduced::Absent
    }
}

/// Whether the junction set was the same in the last three replicates, i.e. it did not change
/// over two consecutive replicates.
pub fn stable(sets: &[BTreeSet<u64>]) -> bool {
    let n = sets.len();
    n >= MIN_REPLICATES_TO_STOP && sets[n - 1] == sets[n - 2] && sets[n - 2] == sets[n - 3]
}

/// Why no further replicate is drawn after those in `sets`, or `None` to draw another.
/// At least `requested` replicates are drawn and at most `max`.
pub fn stop_reason(sets: &[BTreeSet<u64>], requested: usize, max: usize) -> Option<&'static str> {
    if sets.len() < requested.min(max) {
        None
    } else if stable(sets) {
        Some("junction set unchanged over two consecutive replicates")
    } else if sets.len() >= max {
        Some("maximum number of replicates reached; junction set still changing")
    } else {
        None
    }
}

/// Write the reads of `inputs` whose name hashes under `fraction` with this salt (the same
/// deterministic draw as `ovasm recruit --target-depth --salt`) to `out` as FASTQ, read by
/// read: nothing is held in memory. Returns reads and bases written.
pub fn downsample(inputs: &[PathBuf], fraction: f64, salt: u64, out: &Path) -> Result<(u64, u64)> {
    let (rx, reader) = read_batches(inputs);
    let mut w = BufWriter::with_capacity(
        1 << 20,
        File::create(out).with_context(|| format!("cannot write {}", out.display()))?,
    );
    let (mut reads, mut bases) = (0u64, 0u64);
    for batch in rx {
        for r in &batch {
            if fraction < 1.0 && unit_hash(&r.id, salt) >= fraction {
                continue;
            }
            write_fastq_record(&mut w, &r.id, &r.seq, r.qual.as_deref())?;
            reads += 1;
            bases += r.seq.len() as u64;
        }
    }
    reader.join().expect("reader thread panicked")?;
    w.flush()?;
    Ok((reads, bases))
}

/// One FASTQ record; FASTA input gets a neutral quality, as in `ovasm recruit`.
pub fn write_fastq_record(
    w: &mut impl Write,
    id: &[u8],
    seq: &[u8],
    qual: Option<&[u8]>,
) -> Result<()> {
    w.write_all(b"@")?;
    w.write_all(id)?;
    w.write_all(b"\n")?;
    w.write_all(seq)?;
    w.write_all(b"\n+\n")?;
    match qual {
        Some(q) => w.write_all(q)?,
        None => w.write_all(&vec![b'I'; seq.len()])?,
    }
    w.write_all(b"\n")?;
    Ok(())
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

    fn graph(segs: &[&[u8]], links: &[(Oriented, Oriented, usize)]) -> Graph {
        Graph {
            names: (0..segs.len()).map(|i| format!("s{i}")).collect(),
            seqs: segs.iter().map(|s| s.to_vec()).collect(),
            links: links.iter().map(|&(a, b, _)| (a, b)).collect(),
            overlaps: links.iter().map(|&(_, _, ov)| ov).collect(),
        }
    }

    #[test]
    fn the_signature_is_the_same_on_either_strand_and_under_renaming() {
        let (a, b) = (seq(400, 1), seq(300, 2));
        // a+ -> b+ in one graph
        let g1 = graph(&[&a, &b], &[((0, true), (1, true), 0)]);
        // the same adjacency in a graph that stores b first and a reverse-complemented: b+ is
        // followed by ... no: a+ -> b+ read backwards is b- -> a-, and with a stored as its
        // reverse complement (call it A), a- is A+: the link is b- -> A+
        let ra = revcomp(&a);
        let g2 = graph(&[&b, &ra], &[((0, false), (1, true), 0)]);
        // and written as the reverse link A- -> b+
        let g3 = graph(&[&b, &ra], &[((1, false), (0, true), 0)]);
        let s1 = links(&g1)[0].signature;
        assert_eq!(s1, links(&g2)[0].signature);
        assert_eq!(s1, links(&g3)[0].signature);
        // the span itself is the junction's 2 x FLANK bases
        let span = junction_span(&g1, (0, true), (1, true), 0, FLANK);
        assert_eq!(span, [&a[300..], &b[..100]].concat());
        assert_eq!(signature(&span), signature(&revcomp(&span)));
        // a different junction (a+ -> b-) has another signature
        let g4 = graph(&[&a, &b], &[((0, true), (1, false), 0)]);
        assert_ne!(s1, links(&g4)[0].signature);
    }

    #[test]
    fn an_overlap_is_counted_once_and_either_strand_gives_the_same_span() {
        // a's last 50 bases are b's first 50
        let a = seq(400, 3);
        let b = [&a[350..], &seq(300, 4)[..]].concat();
        let g = graph(&[&a, &b], &[((0, true), (1, true), 50)]);
        let fwd = junction_span(&g, (0, true), (1, true), 50, FLANK);
        assert_eq!(fwd, [&a[250..], &b[50..150]].concat());
        assert_eq!(fwd.len(), 2 * FLANK + 50);
        let back = junction_span(&g, (1, false), (0, false), 50, FLANK);
        assert_eq!(back, revcomp(&fwd));
        // the same adjacency in a blunt graph (overlap removed from b) has the same middle, but
        // its flanks reach less far: a different span length, so it is matched as contiguous
        // sequence instead (see `reproduced`)
        let blunt = graph(&[&a, &b[50..]], &[((0, true), (1, true), 0)]);
        let lk = links(&blunt);
        let kmers = spelled_kmers(&g, &links(&g));
        assert_eq!(spelled(&lk[0].span, &kmers), Some(true));
    }

    #[test]
    fn short_segments_cut_the_flanks_without_panicking() {
        let (a, b) = (seq(40, 5), seq(30, 6));
        let g = graph(&[&a, &b], &[((0, true), (1, false), 0)]);
        let span = junction_span(&g, (0, true), (1, false), 0, FLANK);
        assert_eq!(span, [&a[..], &revcomp(&b)[..]].concat());
        let back = junction_span(&g, (1, true), (0, false), 0, FLANK);
        assert_eq!(back, revcomp(&span));
    }

    #[test]
    fn only_links_with_a_branching_end_are_junctions() {
        // a repeat r entered from x and y, left to z and w; c is a circle closed on itself
        let (x, y, r, z, w, c) = (
            seq(300, 10),
            seq(300, 11),
            seq(300, 12),
            seq(300, 13),
            seq(300, 14),
            seq(300, 15),
        );
        let g = graph(
            &[&x, &y, &r, &z, &w, &c],
            &[
                ((0, true), (2, true), 0),
                ((1, true), (2, true), 0),
                ((2, true), (3, true), 0),
                ((2, true), (4, true), 0),
                ((5, true), (5, true), 0),
                // the first link again, written from the other strand: the same adjacency
                ((2, false), (0, false), 0),
            ],
        );
        let lk = links(&g);
        assert_eq!(lk.len(), 5, "the reverse-complement duplicate is one link");
        assert_eq!(
            lk.iter().map(|j| j.branching).collect::<Vec<_>>(),
            [true, true, true, true, false]
        );
        assert_eq!(junction_set(&lk).len(), 4);
    }

    #[test]
    fn a_junction_resolved_inside_a_longer_segment_is_contiguous_and_a_missing_one_absent() {
        let (x, r, z, other) = (seq(500, 20), seq(500, 21), seq(500, 22), seq(500, 23));
        // representative graph: x -> r and r -> z, r also reached from `other`
        let rep = graph(
            &[&x, &r, &z, &other],
            &[
                ((0, true), (1, true), 0),
                ((3, true), (1, true), 0),
                ((1, true), (2, true), 0),
            ],
        );
        let rep_links = links(&rep);
        // replicate: x r z in one piece (its k spans r); nothing about `other`
        let xrz = [&x[..], &r[..], &z[..]].concat();
        let replicate = graph(&[&xrz], &[]);
        let rl = links(&replicate);
        let (set, kmers) = (junction_set(&rl), spelled_kmers(&replicate, &rl));
        let got: Vec<Reproduced> = rep_links
            .iter()
            .map(|j| reproduced(j.signature, &j.span, &set, &kmers))
            .collect();
        assert_eq!(
            got,
            [
                Reproduced::Contiguous,
                Reproduced::Absent,
                Reproduced::Contiguous
            ]
        );
        // the representative graph against itself: every junction as a junction
        let (set, kmers) = (junction_set(&rep_links), spelled_kmers(&rep, &rep_links));
        assert!(rep_links.iter().filter(|j| j.branching).all(|j| reproduced(
            j.signature,
            &j.span,
            &set,
            &kmers
        ) == Reproduced::Junction));
    }

    #[test]
    fn a_span_with_an_unknown_base_is_not_judged() {
        let mut s = seq(200, 30);
        let kmers: FxHashSet<u64> = CanonicalKmers::new(&s, KMER).map(|(_, k)| k).collect();
        assert_eq!(spelled(&s, &kmers), Some(true));
        assert_eq!(spelled(&s[..20], &kmers), None);
        s[100] = b'N';
        assert_eq!(spelled(&s, &kmers), None);
    }

    fn set(v: &[u64]) -> BTreeSet<u64> {
        v.iter().copied().collect()
    }

    #[test]
    fn replicates_stop_when_three_in_a_row_agree_or_at_the_maximum() {
        let (a, b) = (set(&[1, 2, 3]), set(&[1, 2]));
        // fewer than requested: go on, whatever the sets
        assert_eq!(stop_reason(&[a.clone(), a.clone()], 3, 5), None);
        // three equal sets: unchanged over two consecutive replicates
        let three = [a.clone(), a.clone(), a.clone()];
        assert!(stable(&three));
        assert!(stop_reason(&three, 3, 5).unwrap().contains("unchanged"));
        // the third differs: a fourth is drawn
        let changed = [a.clone(), a.clone(), b.clone()];
        assert!(!stable(&changed));
        assert_eq!(stop_reason(&changed, 3, 5), None);
        // b, b after it: two transitions without change are needed, so a fifth
        let four = [a.clone(), a.clone(), b.clone(), b.clone()];
        assert_eq!(stop_reason(&four, 3, 5), None);
        let five = [a.clone(), a.clone(), b.clone(), b.clone(), b.clone()];
        assert!(stop_reason(&five, 3, 5).unwrap().contains("unchanged"));
        // never settling: stops at the maximum and says so
        let restless = [a.clone(), b.clone(), a.clone(), b.clone(), a.clone()];
        assert_eq!(stop_reason(&restless[..4], 3, 5), None);
        assert!(stop_reason(&restless, 3, 5).unwrap().contains("maximum"));
        // more requested than three: no early stop before them
        assert_eq!(stop_reason(&[a.clone(), a.clone(), a.clone()], 4, 5), None);
        // a maximum of two can never show stability
        assert!(stop_reason(&[a.clone(), a.clone()], 2, 2)
            .unwrap()
            .contains("maximum"));
        // empty junction sets (a single circle) agree as well
        let none = [set(&[]), set(&[]), set(&[])];
        assert!(stop_reason(&none, 3, 5).unwrap().contains("unchanged"));
    }

    #[test]
    fn downsampling_is_reproducible_and_salted() {
        let dir = std::env::temp_dir().join(format!("ovasm-replicate-ds-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let reads = dir.join("reads.fa");
        let text: String = (0..4000)
            .map(|i| format!(">r{i} extra\nACGTACGTAC\n"))
            .collect();
        std::fs::write(&reads, text).unwrap();
        let names = |p: &Path| -> Vec<String> {
            std::fs::read_to_string(p)
                .unwrap()
                .lines()
                .filter(|l| l.starts_with('@'))
                .map(String::from)
                .collect()
        };
        let draw = |salt: u64, name: &str| {
            let out = dir.join(name);
            let (n, bases) = downsample(std::slice::from_ref(&reads), 0.5, salt, &out).unwrap();
            let got = names(&out);
            assert_eq!(got.len() as u64, n);
            assert_eq!(bases, n * 10);
            got
        };
        let (a, a2, b) = (draw(1, "a.fq"), draw(1, "a2.fq"), draw(2, "b.fq"));
        assert_eq!(a, a2, "the same salt draws the same reads");
        assert_ne!(a, b, "another salt draws others");
        assert!((a.len() as f64 / 4000.0 - 0.5).abs() < 0.05, "{}", a.len());
        let both = a.iter().filter(|x| b.contains(x)).count() as f64 / 4000.0;
        assert!((both - 0.25).abs() < 0.05, "independent draws: {both}");
        // a fraction of one keeps every read, in order
        let all = dir.join("all.fq");
        assert_eq!(downsample(&[reads], 1.0, 9, &all).unwrap().0, 4000);
        std::fs::remove_dir_all(&dir).ok();
    }
}
