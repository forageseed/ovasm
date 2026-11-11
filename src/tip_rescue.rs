//! Reads at the open ends of an assembled graph that seed recruitment set aside as NUMT
//! candidates (`find recruit`'s partial reads), taken back to assemble across the open ends.
//!
//! Recruitment takes a read when its bait covers most of it, and extends its bait only by reads
//! at about the organelle's depth. A mitochondrial stretch that the plastid also holds (an MTPT)
//! is many times deeper than the mitochondrion in the reads, so the reads running into it from
//! either side stay partial and recruitment stops at both of its edges; the graph then has two
//! open ends with the missing stretch between them. Filipendula ulmaria (DToL HiFi, Oatk
//! reference 865,372 bp): 11,475 bp of its mitochondrion lie between two open ends, an MTPT
//! whose reads, at about the mitochondrion's depth, were all in the NUMT candidate pool (207
//! reads over it, 10 across it with mitochondrial sequence on both sides); with them the graph
//! closed into one molecule holding the whole reference.
//!
//! The rescue walks from each open end: the last `TIP_WINDOW` bases of the unitig at an end are
//! the first frontier, a pool read sharing at least `MIN_SHARED` k-mers with the frontier is
//! taken, and the k-mers of a taken read that the graph lacks join the frontier, so the next
//! round reaches the reads beyond it. k-mers of the other organelle's reads (at its depth
//! among them, `other_kmers`) never join the frontier: the walk crosses an MTPT only on reads holding
//! mitochondrial sequence beside it, not from one plastid-like stretch to every other (with
//! them, Filipendula's walk took 696 reads, 11.1 Mb, in three rounds and hit the cap). Nor do
//! k-mers that more pool reads hold than `POOL_DEPTH_FOLD` times the organelle's read depth:
//! the pool also holds thousands of reads of high-copy sequence (nuclear plastid insertions,
//! the reads at Filipendula's other MTPTs: 1,000-2,900 pool reads deep against 134x), and a
//! read reaching one of them took 2,833 reads in one round of the replayed walk. And a new
//! k-mer joins only when more of the round's reads hold it than `SUPPORT_NUCLEAR_FOLD` times
//! the nuclear depth: the pool holds NUMT reads, one nuclear locus to a few of them, whose
//! nuclear flanks lead to more nuclear reads (the replayed Filipendula walk: 91 reads that
//! map nowhere on the mitochondrion by its third round, 1,097 reads in its fourth); a gap of
//! the organelle has the organelle's reads across it. It stops when a round takes nothing, after `MAX_ROUNDS`,
//! or at `MAX_SHARE` of the recruited bases. Whether the reads are kept is decided by the
//! assembly they give (`pipeline`): fewer open ends, else the first assembly stays.

use crate::kmer::CanonicalKmers;
use anyhow::{Context, Result};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

/// k-mer length of the rescue.
pub const RESCUE_K: usize = 31;
/// Bases at each open end of the graph the walk starts from.
pub const TIP_WINDOW: usize = 5000;
/// k-mers a read must share with the frontier to be taken: a 1 kb block, as a NUMT candidate's
/// least organellar block (`find recruit --numt-min-block`).
pub const MIN_SHARED: usize = 1000;
/// Most rounds of the walk.
pub const MAX_ROUNDS: usize = 20;
/// Most bases the rescue takes, as a share of the recruited reads' bases: a walk that keeps
/// growing past this is in repeats of the nucleus, not across a gap of the organelle.
pub const MAX_SHARE: f64 = 0.1;
/// k-mers held by more pool reads than this many times the organelle's read depth are
/// high-copy sequence, not the way across a gap (`collect`).
pub const POOL_DEPTH_FOLD: f64 = 2.0;
/// Least count among the other organelle's reads for a k-mer to be that organelle's.
pub const OTHER_MIN_COUNT: u32 = 3;
/// A new k-mer joins the frontier when more than this many times the nuclear depth of the
/// round's reads hold it (`collect`), and at least `MIN_SUPPORT` of them.
pub const SUPPORT_NUCLEAR_FOLD: f64 = 1.5;
pub const MIN_SUPPORT: u32 = 3;
/// Share of the other organelle's depth among its reads (the count-weighted median count of
/// their k-mers) a k-mer must reach there to be that organelle's: its reads also hold a few of
/// this organelle's reads, at this organelle's depth. Convolvulus arvensis: every k-mer of both
/// open ends of the mitochondrion stood at least 3 times among the plastid's 3,943 reads, the
/// plastid at about 250x, and the walk had nothing to start from.
pub const OTHER_DEPTH_SHARE: f64 = 0.25;
/// Reads held in memory at a time while the pool is scanned.
const CHUNK: usize = 4096;

#[derive(Debug, Default, Clone, PartialEq)]
pub struct Rescue {
    /// Open ends of the graph the walk started from.
    pub open_ends: usize,
    pub rounds: usize,
    pub reads: usize,
    pub bases: u64,
    /// The walk reached `MAX_SHARE` of the recruited bases (or `MAX_ROUNDS`) and stopped.
    pub capped: bool,
}

/// The sequence at each open end of a GFA (last `window` bases towards the end) and every
/// k-mer of its segments.
fn open_end_windows(gfa: &str, window: usize) -> (Vec<Vec<u8>>, FxHashSet<u64>) {
    let mut segs: Vec<(String, Vec<u8>)> = Vec::new();
    let mut used: FxHashSet<(String, bool)> = FxHashSet::default();
    for line in gfa.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        match f.first() {
            Some(&"S") if f.len() > 2 && f[2] != "*" => {
                segs.push((f[1].to_string(), f[2].as_bytes().to_vec()));
            }
            // L a o1 b o2: a's right end when o1 is +, b's left end when o2 is +
            Some(&"L") if f.len() > 4 => {
                used.insert((f[1].to_string(), f[2] == "+"));
                used.insert((f[3].to_string(), f[4] != "+"));
            }
            _ => {}
        }
    }
    let mut graph = FxHashSet::default();
    let mut windows = Vec::new();
    for (name, seq) in &segs {
        graph.extend(CanonicalKmers::new(seq, RESCUE_K).map(|(_, km)| km));
        let w = window.min(seq.len());
        if !used.contains(&(name.clone(), false)) {
            windows.push(seq[..w].to_vec());
        }
        if !used.contains(&(name.clone(), true)) {
            windows.push(seq[seq.len() - w..].to_vec());
        }
    }
    (windows, graph)
}

/// Remove from `kmers` those that more than `max` reads of `pool` hold.
fn drop_pool_deep(pool: &Path, kmers: &mut FxHashSet<u64>, max: u32) -> Result<()> {
    let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
    let mut reader = needletail::parse_fastx_file(pool)
        .with_context(|| format!("cannot read {}", pool.display()))?;
    while let Some(rec) = reader.next() {
        let rec = rec.context("bad record in the NUMT candidate pool")?;
        let seen: FxHashSet<u64> = CanonicalKmers::new(&rec.seq(), RESCUE_K)
            .map(|(_, km)| km)
            .filter(|km| kmers.contains(km))
            .collect();
        for km in seen {
            *counts.entry(km).or_insert(0) += 1;
        }
    }
    kmers.retain(|km| counts.get(km).map_or(true, |&c| c <= max));
    Ok(())
}

fn shared(seq: &[u8], frontier: &FxHashSet<u64>) -> usize {
    CanonicalKmers::new(seq, RESCUE_K)
        .filter(|(_, km)| frontier.contains(km))
        .count()
}

/// The k-mers of the reads in `files` at their organelle's depth: counted at least
/// `OTHER_DEPTH_SHARE` of the count-weighted median count, and `OTHER_MIN_COUNT`.
fn other_kmers(files: &[PathBuf]) -> Result<FxHashSet<u64>> {
    let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
    for f in files {
        let mut reader = needletail::parse_fastx_file(f)
            .with_context(|| format!("cannot read {}", f.display()))?;
        while let Some(rec) = reader.next() {
            let rec = rec.context("bad record")?;
            for (_, km) in CanonicalKmers::new(&rec.seq(), RESCUE_K) {
                *counts.entry(km).or_insert(0) += 1;
            }
        }
    }
    let min = OTHER_MIN_COUNT.max((OTHER_DEPTH_SHARE * weighted_median(&counts) as f64).ceil() as u32);
    Ok(counts.into_iter().filter(|c| c.1 >= min).map(|c| c.0).collect())
}

/// The count at which half the k-mer occurrences stand (0 when there are none).
fn weighted_median(counts: &FxHashMap<u64, u32>) -> u32 {
    let mut c: Vec<u32> = counts.values().copied().collect();
    c.sort_unstable();
    let total: u64 = c.iter().map(|&x| u64::from(x)).sum();
    let mut seen = 0u64;
    for x in c {
        seen += u64::from(x);
        if 2 * seen >= total {
            return x;
        }
    }
    0
}

/// Walk from the open ends of `gfa` through the reads of `pool` (FASTQ/FASTA), writing the
/// reads taken to `out` (FASTQ). The k-mers of the reads in `other` (the other organelle's)
/// are left out of the walk, and with `read_depth` (the organelle's) so are k-mers more pool
/// reads hold than `POOL_DEPTH_FOLD` times it; with `nuclear_depth`, a new k-mer joins only
/// when more of the round's reads hold it than `SUPPORT_NUCLEAR_FOLD` times that (else
/// `MIN_SUPPORT`); `max_bases` bounds what is taken.
pub fn collect(
    gfa: &Path,
    pool: &Path,
    other: &[PathBuf],
    read_depth: Option<f64>,
    nuclear_depth: Option<f64>,
    out: &Path,
    max_bases: u64,
) -> Result<Rescue> {
    let text = fs::read_to_string(gfa).with_context(|| format!("cannot read {}", gfa.display()))?;
    let (windows, mut graph) = open_end_windows(&text, TIP_WINDOW);
    let mut rescue = Rescue { open_ends: windows.len(), ..Default::default() };
    if windows.is_empty() {
        File::create(out).with_context(|| format!("cannot write {}", out.display()))?;
        return Ok(rescue);
    }
    // the other organelle's k-mers: neither frontier nor added to it
    let other = other_kmers(other)?;
    let mut frontier: FxHashSet<u64> = FxHashSet::default();
    for w in &windows {
        frontier.extend(
            CanonicalKmers::new(w, RESCUE_K)
                .map(|(_, km)| km)
                .filter(|km| !other.contains(km)),
        );
    }
    graph.extend(other);
    let mut taken_idx: FxHashSet<usize> = FxHashSet::default();
    let mut taken: Vec<(Vec<u8>, Vec<u8>, Option<Vec<u8>>)> = Vec::new();
    while !frontier.is_empty() && rescue.rounds < MAX_ROUNDS && !rescue.capped {
        rescue.rounds += 1;
        let mut round: Vec<(usize, Vec<u8>, Vec<u8>, Option<Vec<u8>>)> = Vec::new();
        let mut reader = needletail::parse_fastx_file(pool)
            .with_context(|| format!("cannot read {}", pool.display()))?;
        let mut chunk: Vec<(usize, Vec<u8>, Vec<u8>, Option<Vec<u8>>)> = Vec::with_capacity(CHUNK);
        let mut i = 0usize;
        let scan = |chunk: &mut Vec<(usize, Vec<u8>, Vec<u8>, Option<Vec<u8>>)>,
                        round: &mut Vec<(usize, Vec<u8>, Vec<u8>, Option<Vec<u8>>)>| {
            let hits: Vec<bool> = chunk
                .par_iter()
                .map(|r| shared(&r.2, &frontier) >= MIN_SHARED)
                .collect();
            for (r, h) in chunk.drain(..).zip(hits) {
                if h {
                    round.push(r);
                }
            }
        };
        while let Some(rec) = reader.next() {
            let rec = rec.context("bad record in the NUMT candidate pool")?;
            if !taken_idx.contains(&i) {
                chunk.push((i, rec.id().to_vec(), rec.seq().to_vec(), rec.qual().map(<[u8]>::to_vec)));
                if chunk.len() == CHUNK {
                    scan(&mut chunk, &mut round);
                }
            }
            i += 1;
        }
        scan(&mut chunk, &mut round);
        if round.is_empty() {
            break;
        }
        let support = nuclear_depth
            .map_or(MIN_SUPPORT, |d| MIN_SUPPORT.max((SUPPORT_NUCLEAR_FOLD * d).floor() as u32 + 1));
        let mut held: FxHashMap<u64, u32> = FxHashMap::default();
        for (idx, id, seq, qual) in round {
            if rescue.bases + seq.len() as u64 > max_bases {
                rescue.capped = true;
                break;
            }
            rescue.bases += seq.len() as u64;
            rescue.reads += 1;
            taken_idx.insert(idx);
            let novel: FxHashSet<u64> = CanonicalKmers::new(&seq, RESCUE_K)
                .map(|(_, km)| km)
                .filter(|km| !graph.contains(km))
                .collect();
            for km in novel {
                *held.entry(km).or_insert(0) += 1;
            }
            taken.push((id, seq, qual));
        }
        let mut grown: FxHashSet<u64> =
            held.into_iter().filter(|h| h.1 >= support).map(|h| h.0).collect();
        if let Some(d) = read_depth {
            drop_pool_deep(pool, &mut grown, (POOL_DEPTH_FOLD * d).ceil() as u32)?;
        }
        frontier.extend(grown);
    }
    if rescue.rounds == MAX_ROUNDS {
        rescue.capped = true;
    }
    let mut w = BufWriter::new(File::create(out).with_context(|| format!("cannot write {}", out.display()))?);
    for (id, seq, qual) in &taken {
        w.write_all(b"@")?;
        w.write_all(id)?;
        w.write_all(b"\n")?;
        w.write_all(seq)?;
        w.write_all(b"\n+\n")?;
        match qual {
            Some(q) if q.len() == seq.len() => w.write_all(q)?,
            _ => w.write_all(&vec![b'I'; seq.len()])?,
        }
        w.write_all(b"\n")?;
    }
    w.flush()?;
    Ok(rescue)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic pseudo-random DNA.
    fn dna(n: usize, mut s: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                b"ACGT"[(s % 4) as usize]
            })
            .collect()
    }

    fn fastq(reads: &[&[u8]]) -> String {
        reads
            .iter()
            .enumerate()
            .map(|(i, r)| {
                format!("@r{i}\n{}\n+\n{}\n", String::from_utf8_lossy(r), "I".repeat(r.len()))
            })
            .collect()
    }

    #[test]
    fn walks_from_both_open_ends_across_the_gap_and_leaves_unrelated_reads() {
        // a linear graph: one unitig of 20 kb, the genome 30 kb with a 10 kb gap after it
        let genome = dna(30_000, 7);
        let unitig = &genome[..20_000];
        let gfa = format!("H\tVN:Z:1.0\nS\tu0\t{}\n", String::from_utf8_lossy(unitig));
        // pool: reads of 6 kb tiling 17-30 kb every 1.5 kb, and unrelated (nuclear) reads
        let reads: Vec<Vec<u8>> = (0..)
            .map(|i| 17_000 + i * 1500)
            .take_while(|&s| s + 6000 <= 30_000)
            .map(|s| genome[s..s + 6000].to_vec())
            .collect();
        // three reads at each place: the organelle's depth, above `MIN_SUPPORT`
        let mut reads: Vec<Vec<u8>> = reads.iter().flat_map(|r| [r.clone(), r.clone(), r.clone()]).collect();
        let gap_reads = reads.len();
        reads.push(dna(8000, 99));
        reads.push(dna(8000, 123));
        let dir = std::env::temp_dir().join(format!("ovasm-tip-rescue-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let (g, p, o) = (dir.join("a.gfa"), dir.join("pool.fq"), dir.join("out.fq"));
        fs::write(&g, gfa).unwrap();
        let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
        fs::write(&p, fastq(&refs)).unwrap();
        let r = collect(&g, &p, &[], None, None, &o, u64::MAX).unwrap();
        assert_eq!(r.open_ends, 2);
        assert_eq!(r.reads, gap_reads, "{r:?}");
        assert!(!r.capped);
        assert!(r.rounds >= 2, "the walk needs more than one round: {r:?}");
        // the gap's middle held by the other organelle too (an MTPT, 20.5-26.5 kb, longer than a read less two blocks): the walk
        // crosses it only on reads with sequence of its own beyond it, and stops there
        let other = dir.join("other.fq");
        let mtpt = &genome[20_500..26_500];
        // the other organelle's reads: 20 of the MTPT, and 3 of this organelle's own reads
        // (one at the open end), which do not make its sequence the other's
        let mut theirs: Vec<&[u8]> = vec![mtpt; 20];
        theirs.extend([&genome[16_000..22_000], &genome[2_000..8_000], &genome[9_000..15_000]]);
        theirs.extend([&genome[16_000..22_000], &genome[16_000..22_000]]);
        fs::write(&other, fastq(&theirs)).unwrap();
        let m = collect(&g, &p, &[other], None, None, &o, u64::MAX).unwrap();
        assert!(m.reads < gap_reads && m.reads >= 2, "{m:?}");
        // the cap stops the walk
        let capped = collect(&g, &p, &[], None, None, &o, 7000).unwrap();
        assert!(capped.capped && capped.reads == 1, "{capped:?}");
        // a closed graph has no open end and takes nothing
        let circle = format!(
            "S\tu0\t{}\nL\tu0\t+\tu0\t+\t0M\n",
            String::from_utf8_lossy(unitig)
        );
        fs::write(&g, circle).unwrap();
        // a high-copy stretch beyond the end (here 40 pool reads of it, against a read depth
        // of 4): the walk takes the reads at the end but goes no further through it
        let deep = &genome[20_000..23_000];
        let mut with_deep = reads.clone();
        with_deep.extend((0..40).map(|_| deep.to_vec()));
        let refs: Vec<&[u8]> = with_deep.iter().map(Vec::as_slice).collect();
        fs::write(&p, fastq(&refs)).unwrap();
        let d = collect(&g, &p, &[], Some(4.0), None, &o, u64::MAX).unwrap();
        assert!(d.reads < gap_reads + 40, "{d:?}");
        let none = collect(&g, &p, &[], None, None, &o, u64::MAX).unwrap();
        assert_eq!((none.open_ends, none.reads), (0, 0));
        fs::remove_dir_all(&dir).ok();
    }
}
