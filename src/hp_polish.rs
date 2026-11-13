//! Homopolymer run lengths of noisy-read assemblies from the uncorrected reads
//! (`ovasm assemble polish`; `run ont` and `run hybrid` polish their molecules with it).
//!
//! Nanopore reads call homopolymer runs short far more often than long, at every strand alike, and
//! the k-mer self-correction keeps a read's short call where that call is itself a recurring,
//! solid k-mer: the molecules of Nipponbare ONT came out at 0.998 identity, 84% of the corrected
//! reads' remaining indels homopolymer deletions. Yet over many reads the most frequent call is the
//! true one: against the Nipponbare HiFi plastid, the mode of the raw ONT calls gave the true run
//! length at 99.94% of 13,916 runs read 10 times or more (length 8: 60% of calls exact, 34%
//! short; length 10: 32% exact, 54% short spread over several lengths). The rounded mean, as
//! `hpc::restore` takes, is biased short (length 10: mean 9.1).
//!
//! Each molecule is homopolymer-compressed. Every raw read, compressed too, is placed on it by the
//! exact compressed k-mers it shares with it (both orientations), each word extended along its
//! diagonal to the longest exact stretch. Compressed sequence is of low complexity and a lone
//! 19-word recurs: words alone put calls of other places on a run (Nipponbare plastid, A×5 at
//! 49,465: 35 of 70 word calls were "1", none by alignment); only a stretch of `STRETCH` exact
//! compressed bases places a read. Inside it the read's run at each position (the stretch's two
//! end runs aside) is an observation of the molecule's run there, once per read and position; a
//! read's first and last run are cut by its ends and do not count.
//!
//! A run read `MIN_OBSERVATIONS` times or more is judged. The most frequent call is still short at
//! the long runs (plastid T×7 at 7,223: 32 calls of 6, 30 of 7), so the calls are taken by
//! likelihood with the reads' own error model: P(call | true length), for A/T and G/C runs apart,
//! learned by hard EM from the runs' current lengths (first the mode), each run then given the
//! length of highest prior × likelihood, until no run changes. Short runs, read nearly always
//! exactly, anchor the model: the plastid's model gave length 7 read as 7 at 56%, as 6 at 21%;
//! homopolymer indels against the HiFi plastid went 236 (unpolished) → 110 (mode) → 51 (EM).
//! The compressed sequence itself is never changed.

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Result};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;

use crate::recruit::read_batches;

/// Length of the exact compressed words that place a read (a compressed nanopore read of 91%
/// identity keeps about half of its 19-mers exact).
pub const WORD: usize = 19;
/// Fewest observations for a run to be judged.
pub const MIN_OBSERVATIONS: u32 = 10;
/// Reads are used until their bases reach this many times the molecules' length.
pub const MAX_DEPTH: f64 = 200.0;
/// Longest run length told apart; longer calls count here.
const MAX_RUN: usize = 63;
/// A read's exact word places it only with this many words (itself included) on diagonals at
/// most `DIAGONAL_SLACK` compressed bases apart (indels other than homopolymer ones shift it).
const MIN_WORDS_ON_DIAGONAL: usize = 3;
const DIAGONAL_SLACK: i64 = 20;
/// Fewest exact compressed bases in a row for a read's stretch to place it (40 and 80 placed the
/// plastid's runs alike; a lone 19-word misplaces).
pub const STRETCH: usize = 40;
/// Most rounds of the error model's EM.
const EM_ROUNDS: usize = 12;

#[derive(Default, Serialize, Clone)]
pub struct MoleculeReport {
    pub name: String,
    pub length_before: usize,
    pub length_after: usize,
    pub runs: usize,
    pub runs_judged: usize,
    pub runs_lengthened: usize,
    pub runs_shortened: usize,
    pub bases_added: usize,
    pub bases_removed: usize,
}

#[derive(Default, Serialize)]
pub struct PolishReport {
    pub reads: u64,
    pub read_bases: u64,
    pub reads_placed: u64,
    pub observations: u64,
    pub min_observations: u32,
    /// runs whose length the likelihood set apart from the most frequent call
    pub runs_moved_from_mode: usize,
    pub em_rounds: usize,
    /// P(call = L | true length L), A/T runs, lengths 1 to 12 (the learned model)
    pub at_exact_call_rate: Vec<f64>,
    pub molecules: Vec<MoleculeReport>,
    pub elapsed_seconds: f64,
}

struct Molecule {
    name: String,
    header: String,
    circular: bool,
    bases: Vec<u8>,
    runs: Vec<u32>,
    /// observation counts per position and run length
    seen: Vec<[u32; MAX_RUN + 1]>,
}

fn code(b: u8) -> Option<u64> {
    match b {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// 2-bit words of `WORD` compressed bases with their start positions (windows holding N skipped).
fn words(seq: &[u8]) -> Vec<(usize, u64)> {
    let mask = (1u64 << (2 * WORD)) - 1;
    let mut out = Vec::new();
    let (mut v, mut valid) = (0u64, 0usize);
    for (i, &b) in seq.iter().enumerate() {
        match code(b) {
            Some(c) => {
                v = ((v << 2) | c) & mask;
                valid += 1;
                if valid >= WORD {
                    out.push((i + 1 - WORD, v));
                }
            }
            None => valid = 0,
        }
    }
    out
}

fn revcomp_hpc(bases: &[u8], runs: &[u32]) -> (Vec<u8>, Vec<u32>) {
    (crate::evidence::revcomp(bases), runs.iter().rev().copied().collect())
}

/// Read a FASTA of molecules (header kept; `circular=true` in it marks a circle).
fn read_molecules(path: &Path) -> Result<Vec<Molecule>> {
    let text = fs::read_to_string(path)?;
    let mut out = Vec::new();
    let mut cur: Option<(String, Vec<u8>)> = None;
    let push = |h: String, s: Vec<u8>, out: &mut Vec<Molecule>| -> Result<()> {
        let (bases, runs) = crate::hpc::compress(&s)?;
        let n = bases.len();
        out.push(Molecule {
            name: h.split_whitespace().next().unwrap_or("").to_string(),
            circular: h.contains("circular=true"),
            header: h,
            bases,
            runs,
            seen: vec![[0; MAX_RUN + 1]; n],
        });
        Ok(())
    };
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if let Some((h0, s0)) = cur.take() {
                push(h0, s0, &mut out)?;
            }
            cur = Some((h.to_string(), Vec::new()));
        } else if let Some((_, s)) = cur.as_mut() {
            s.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some((h0, s0)) = cur.take() {
        push(h0, s0, &mut out)?;
    }
    Ok(out)
}

/// The run length a position takes: the most frequent call when read `min_obs` times or more.
fn judge(current: u32, seen: &[u32; MAX_RUN + 1], min_obs: u32) -> Option<u32> {
    let total: u32 = seen.iter().sum();
    if total < min_obs {
        return None;
    }
    let top = *seen.iter().max().unwrap();
    if (current as usize) <= MAX_RUN && seen[current as usize] == top {
        return Some(current);
    }
    // the longest of the most frequent: nanopore errors are short calls
    let best = (1..=MAX_RUN).rev().find(|&l| seen[l] == top)?;
    Some(best as u32)
}

/// Run lengths by likelihood under the reads' own error model (hard EM). `runs` holds each judged
/// run's class (0 A/T, 1 G/C), its calls and its starting length (the mode); returns the lengths,
/// the rounds taken, and the model's exact-call rate of A/T lengths 1 to 12.
fn em_lengths(runs: &[(usize, &[u32; MAX_RUN + 1])], start: &[u32]) -> (Vec<u32>, usize, Vec<f64>) {
    type Table = Vec<[f64; MAX_RUN + 1]>;
    let mut call = start.to_vec();
    let mut rounds = 0;
    let mut rate = Vec::new();
    let fit = |call: &[u32]| -> (Vec<Table>, Vec<Vec<f64>>, Vec<[f64; 13]>) {
        // counts of calls per (class, true length), and of lengths per class
        let mut counts = vec![vec![[0u64; MAX_RUN + 1]; MAX_RUN + 1]; 2];
        let mut prior = vec![vec![0u64; MAX_RUN + 1]; 2];
        for ((c, seen), &l) in runs.iter().zip(call) {
            let l = (l as usize).min(MAX_RUN);
            prior[*c][l] += 1;
            for (o, &v) in seen.iter().enumerate() {
                counts[*c][l][o] += v as u64;
            }
        }
        let mut exact = vec![[0.0; 13]; 2];
        let mut logp: Vec<Table> = vec![vec![[0.0; MAX_RUN + 1]; MAX_RUN + 1]; 2];
        let mut logprior = vec![vec![f64::NEG_INFINITY; MAX_RUN + 1]; 2];
        for c in 0..2 {
            let runs_c: u64 = prior[c].iter().sum();
            for l in 1..=MAX_RUN {
                if prior[c][l] > 0 {
                    logprior[c][l] = (prior[c][l] as f64 / runs_c as f64).ln();
                }
                // add-half smoothing: a call never seen at this length is unlikely, not impossible
                let n: u64 = counts[c][l].iter().sum();
                let denom = n as f64 + 0.5 * (MAX_RUN + 1) as f64;
                for o in 0..=MAX_RUN {
                    logp[c][l][o] = ((counts[c][l][o] as f64 + 0.5) / denom).ln();
                }
                if l <= 12 && n > 0 {
                    exact[c][l] = counts[c][l][l] as f64 / n as f64;
                }
            }
        }
        (logp, logprior, exact)
    };
    while rounds < EM_ROUNDS {
        let (logp, logprior, exact) = fit(&call);
        rate = exact[0][1..].to_vec();
        let next: Vec<u32> = runs
            .par_iter()
            .zip(&call)
            .map(|((c, seen), &cur)| {
                let lo = seen.iter().position(|&v| v > 0).unwrap_or(1).max(2) - 1;
                let hi = (seen.iter().rposition(|&v| v > 0).unwrap_or(MAX_RUN) + 1).min(MAX_RUN);
                let score = |l: usize| {
                    logprior[*c][l] + seen.iter().enumerate().map(|(o, &v)| v as f64 * logp[*c][l][o]).sum::<f64>()
                };
                // the current length wins a tie
                let mut best = (score(cur as usize), cur);
                for l in lo..=hi {
                    let s = score(l);
                    if s > best.0 {
                        best = (s, l as u32);
                    }
                }
                best.1
            })
            .collect();
        rounds += 1;
        let changed = next.iter().zip(&call).filter(|(a, b)| a != b).count();
        call = next;
        if changed == 0 {
            break;
        }
    }
    (call, rounds, rate)
}

/// Polish the homopolymer runs of the molecules in `fasta` with the uncorrected `reads`; write
/// them to `out` (same headers, `length=` updated) and return the report.
pub fn polish(fasta: &Path, reads: &[PathBuf], out: &Path, min_obs: u32) -> Result<PolishReport> {
    let t0 = Instant::now();
    let mut mols = read_molecules(fasta)?;
    if mols.is_empty() {
        bail!("{} holds no molecule", fasta.display());
    }
    let total_len: usize = mols.iter().map(|m| m.runs.iter().map(|&r| r as usize).sum::<usize>()).sum();
    // index: word -> (molecule, position); a circle's words run over its end
    let mut index: FxHashMap<u64, Vec<(u32, u32)>> = FxHashMap::default();
    for (mi, m) in mols.iter().enumerate() {
        let n = m.bases.len();
        let ext: Vec<u8> = if m.circular && n >= WORD {
            m.bases.iter().chain(&m.bases[..WORD - 1]).copied().collect()
        } else {
            m.bases.clone()
        };
        for (p, w) in words(&ext) {
            index.entry(w).or_default().push((mi as u32, p as u32));
        }
    }
    let lens: Vec<usize> = mols.iter().map(|m| m.bases.len()).collect();
    let mols_bases: Vec<Vec<u8>> = mols.iter().map(|m| m.bases.clone()).collect();
    let circular_mol: Vec<bool> = mols.iter().map(|m| m.circular).collect();
    // a circle opened inside a run: its first and last compressed bases are one run
    let circular_join: Vec<bool> = mols
        .iter()
        .map(|m| m.circular && m.bases.len() > 1 && m.bases.first() == m.bases.last())
        .collect();
    let budget = (MAX_DEPTH * total_len as f64) as u64;
    let mut rep = PolishReport { min_observations: min_obs, ..Default::default() };
    let (rx, reader) = read_batches(reads);
    let mut stopped = false;
    'batches: for batch in rx {
        // per read: (molecule, position, run length) observations, once per read and position
        let found: Vec<Vec<(u32, u32, u32)>> = batch
            .par_iter()
            .map(|r| {
                let Ok((b, runs)) = crate::hpc::compress(&r.seq) else {
                    return Vec::new();
                };
                let (rb, rruns) = revcomp_hpc(&b, &runs);
                let mut seen: FxHashSet<(u32, u32)> = FxHashSet::default();
                let mut obs = Vec::new();
                for (seq, rl) in [(&b, &runs), (&rb, &rruns)] {
                    let last = seq.len().saturating_sub(1);
                    // (molecule, diagonal, molecule position, read position) of every exact word
                    let mut hits: Vec<(u32, i64, u32, usize)> = Vec::new();
                    for (q, w) in words(seq) {
                        if let Some(at) = index.get(&w) {
                            for &(mi, p) in at {
                                let n = lens[mi as usize] as i64;
                                hits.push((mi, (p as i64 - q as i64).rem_euclid(n), p, q));
                            }
                        }
                    }
                    // a word counts only beside others on its diagonal: a lone exact word of a
                    // read is as likely chance (or a nuclear copy's) as the read's place
                    hits.sort_unstable();
                    // the stretch last taken: (molecule, diagonal, its end on the read); the
                    // words inside it are the same stretch again
                    let mut covered: Option<(u32, i64, usize)> = None;
                    for i in 0..hits.len() {
                        let (mi, d, p, q) = hits[i];
                        if matches!(covered, Some((cm, cd, end)) if cm == mi && cd == d && q < end) {
                            continue;
                        }
                        let near = |j: usize| hits[j].0 == mi && (hits[j].1 - d).abs() <= DIAGONAL_SLACK;
                        let mut support = 1;
                        let mut j = i;
                        while j > 0 && support < MIN_WORDS_ON_DIAGONAL && near(j - 1) {
                            j -= 1;
                            support += 1;
                        }
                        j = i + 1;
                        while j < hits.len() && support < MIN_WORDS_ON_DIAGONAL && near(j) {
                            j += 1;
                            support += 1;
                        }
                        if support < MIN_WORDS_ON_DIAGONAL {
                            continue;
                        }
                        let n = lens[mi as usize];
                        let m = &mols_bases[mi as usize];
                        let circular = circular_mol[mi as usize];
                        // the molecule's base `off` places from p (None off a linear end)
                        let at = |off: i64| -> Option<u8> {
                            let x = p as i64 + off;
                            if circular {
                                Some(m[x.rem_euclid(n as i64) as usize])
                            } else if (0..n as i64).contains(&x) {
                                Some(m[x as usize])
                            } else {
                                None
                            }
                        };
                        // extend the word both ways while read and molecule agree (a circle at
                        // most once round)
                        let (mut a, mut c) = (0usize, WORD);
                        while a + c < n && a < q && at(-(a as i64) - 1) == Some(seq[q - a - 1]) {
                            a += 1;
                        }
                        while a + c < n && q + c < seq.len() && at(c as i64) == Some(seq[q + c]) {
                            c += 1;
                        }
                        covered = Some((mi, d, q + c));
                        if a + c < STRETCH {
                            continue;
                        }
                        let circle_edge = circular_join[mi as usize];
                        // the stretch's own end runs may be an indel's edge
                        for k in (1 - a as i64)..(c as i64 - 1) {
                            let qq = (q as i64 + k) as usize;
                            if qq == 0 || qq == last {
                                continue; // cut by the read's end
                            }
                            let pos = (p as i64 + k).rem_euclid(n as i64) as usize;
                            if circle_edge && (pos == 0 || pos + 1 == n) {
                                continue; // one run cut in two by where the circle is opened
                            }
                            if seen.insert((mi, pos as u32)) {
                                obs.push((mi, pos as u32, rl[qq]));
                            }
                        }
                    }
                }
                obs
            })
            .collect();
        for (r, obs) in batch.iter().zip(found) {
            rep.reads += 1;
            rep.read_bases += r.seq.len() as u64;
            rep.reads_placed += u64::from(!obs.is_empty());
            rep.observations += obs.len() as u64;
            for (mi, pos, len) in obs {
                mols[mi as usize].seen[pos as usize][(len as usize).min(MAX_RUN)] += 1;
            }
            if rep.read_bases >= budget {
                stopped = true;
                break 'batches;
            }
        }
    }
    // stopped early, the reader ends on its own once the channel is gone
    if !stopped {
        reader.join().map_err(|_| anyhow::anyhow!("read parser panicked"))??;
    }
    if let Some(dir) = out.parent() {
        fs::create_dir_all(dir)?;
    }
    // the most frequent call of every judged run, then lengths by likelihood (all molecules share
    // one error model: the same reads)
    let mut judged: Vec<(usize, usize)> = Vec::new();
    let mut classed: Vec<(usize, &[u32; MAX_RUN + 1])> = Vec::new();
    let mut mode: Vec<u32> = Vec::new();
    for (mi, m) in mols.iter().enumerate() {
        for (i, (&b, &cur)) in m.bases.iter().zip(&m.runs).enumerate() {
            if b == b'N' {
                continue;
            }
            if let Some(l) = judge(cur, &m.seen[i], min_obs) {
                judged.push((mi, i));
                classed.push((usize::from(b == b'C' || b == b'G'), &m.seen[i]));
                mode.push(l);
            }
        }
    }
    let (lengths, rounds, rate) = em_lengths(&classed, &mode);
    rep.em_rounds = rounds;
    rep.at_exact_call_rate = rate;
    rep.runs_moved_from_mode = lengths.iter().zip(&mode).filter(|(a, b)| a != b).count();
    let mut call: FxHashMap<(usize, usize), u32> = FxHashMap::default();
    for (k, l) in judged.into_iter().zip(lengths) {
        call.insert(k, l);
    }
    let mut w = BufWriter::new(File::create(out)?);
    for (mi, m) in mols.iter().enumerate() {
        let mut r = MoleculeReport {
            name: m.name.clone(),
            runs: m.runs.len(),
            length_before: m.runs.iter().map(|&x| x as usize).sum(),
            ..Default::default()
        };
        let mut seq = Vec::with_capacity(r.length_before + 64);
        for (i, (&b, &cur)) in m.bases.iter().zip(&m.runs).enumerate() {
            let mut len = cur;
            if b != b'N' {
                if let Some(&l) = call.get(&(mi, i)) {
                    r.runs_judged += 1;
                    if l > cur {
                        r.runs_lengthened += 1;
                        r.bases_added += (l - cur) as usize;
                    } else if l < cur {
                        r.runs_shortened += 1;
                        r.bases_removed += (cur - l) as usize;
                    }
                    len = l;
                }
            }
            seq.extend(std::iter::repeat(b).take(len as usize));
        }
        r.length_after = seq.len();
        let header = header_with_length(&m.header, seq.len());
        writeln!(w, ">{header}")?;
        w.write_all(&seq)?;
        w.write_all(b"\n")?;
        rep.molecules.push(r);
    }
    w.flush()?;
    rep.elapsed_seconds = t0.elapsed().as_secs_f64();
    Ok(rep)
}

/// The header with its `length=` field (if any) set to `len`.
fn header_with_length(h: &str, len: usize) -> String {
    h.split(' ')
        .map(|f| if f.starts_with("length=") { format!("length={len}") } else { f.to_string() })
        .collect::<Vec<_>>()
        .join(" ")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (*seed >> 11) as f64 / (1u64 << 53) as f64
    }

    #[test]
    fn the_most_frequent_call_wins_and_ties_go_long() {
        let mut s = [0u32; MAX_RUN + 1];
        s[7] = 30;
        s[8] = 50;
        s[6] = 20;
        assert_eq!(judge(7, &s, 10), Some(8));
        s[7] = 50;
        assert_eq!(judge(7, &s, 10), Some(7), "the current length among the most frequent stays");
        assert_eq!(judge(6, &s, 10), Some(8), "a tie away from the current goes to the longer");
        let few = { let mut f = [0u32; MAX_RUN + 1]; f[3] = 9; f };
        assert_eq!(judge(4, &few, 10), None, "under the minimum nothing is judged");
    }

    #[test]
    fn the_learned_error_model_recovers_runs_whose_most_frequent_call_is_short() {
        // nanopore-like: length 7 read as 6 a little more often than as 7, so the mode of a true
        // 7 is often 6; lengths 1–6 mostly read right. 30 calls per run.
        let p = |l: usize| -> Vec<(usize, f64)> {
            match l {
                7 => vec![(5, 0.10), (6, 0.42), (7, 0.40), (8, 0.08)],
                6 => vec![(4, 0.08), (5, 0.20), (6, 0.64), (7, 0.08)],
                5 => vec![(4, 0.15), (5, 0.78), (6, 0.07)],
                _ => vec![(l - 1, 0.05), (l, 0.92), (l + 1, 0.03)],
            }
        };
        let mut seed = 11u64;
        let (mut truth, mut seen) = (Vec::new(), Vec::new());
        for i in 0..6000 {
            let l = [2, 2, 3, 3, 4, 5, 5, 6, 6, 7][i % 10];
            let mut s = [0u32; MAX_RUN + 1];
            for _ in 0..30 {
                let mut u = lcg(&mut seed);
                for (o, w) in p(l) {
                    if u < w {
                        s[o] += 1;
                        break;
                    }
                    u -= w;
                }
            }
            truth.push(l as u32);
            seen.push(s);
        }
        let start: Vec<u32> = seen.iter().map(|s| judge(1, s, 1).unwrap()).collect();
        let runs: Vec<(usize, &[u32; MAX_RUN + 1])> = seen.iter().map(|s| (0, s)).collect();
        let (call, rounds, rate) = em_lengths(&runs, &start);
        let right = |c: &[u32]| c.iter().zip(&truth).filter(|(a, b)| a == b).count();
        let sevens = |c: &[u32]| c.iter().zip(&truth).filter(|(a, b)| **b == 7 && a == b).count();
        assert!(sevens(&start) < 500, "the mode misses many 7s: {}", sevens(&start));
        assert!(sevens(&call) >= 570, "likelihood gets the 7s back: {} of 600", sevens(&call));
        assert!(right(&call) >= 5900, "{} of 6000", right(&call));
        assert!(rounds >= 2);
        assert!((rate[6] - 0.40).abs() < 0.05, "learned P(7|7) {}", rate[6]);
    }

    #[test]
    fn short_homopolymer_calls_are_lengthened_back_from_the_reads() {
        // a 3 kb circle with runs up to 9; the "assembly" has every run of 6 or more one short;
        // 120 reads call runs short at a rate growing with length (never the mode), both strands
        let mut seed = 7u64;
        let mut truth = Vec::new();
        while truth.len() < 3000 {
            let b = b"ACGT"[(lcg(&mut seed) * 4.0) as usize % 4];
            if truth.last() == Some(&b) {
                continue;
            }
            let len = if lcg(&mut seed) < 0.1 { 6 + (lcg(&mut seed) * 4.0) as usize } else { 1 + (lcg(&mut seed) * 3.0) as usize };
            truth.extend(std::iter::repeat(b).take(len));
        }
        if truth.first() == truth.last() {
            let other = *b"ACGT".iter().find(|&&c| Some(&c) != truth.first()).unwrap();
            truth.push(other);
        }
        let (tb, tr) = crate::hpc::compress(&truth).unwrap();
        let spell = |runs: &[u32]| -> Vec<u8> {
            tb.iter().zip(runs).flat_map(|(&b, &l)| std::iter::repeat(b).take(l as usize)).collect()
        };
        let damaged: Vec<u32> = tr.iter().map(|&l| if l >= 6 { l - 1 } else { l }).collect();
        let dir = std::env::temp_dir().join(format!("ovasm-hp-polish-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("m.fasta");
        fs::write(&fa, format!(">m circular=true length={}\n{}\n", damaged.iter().sum::<u32>(), String::from_utf8(spell(&damaged)).unwrap())).unwrap();
        let mut fq = String::new();
        for i in 0..120 {
            let start = (lcg(&mut seed) * tr.len() as f64) as usize;
            let runs: Vec<u32> = (0..600)
                .map(|j| {
                    let l = tr[(start + j) % tr.len()];
                    let short = lcg(&mut seed) < 0.05 * l as f64; // length 9: 45% short, still not the mode
                    if short && l > 1 { l - 1 } else { l }
                })
                .collect();
            let bases: Vec<u8> = (0..600).map(|j| tb[(start + j) % tb.len()]).collect();
            let mut read: Vec<u8> = bases.iter().zip(&runs).flat_map(|(&b, &l)| std::iter::repeat(b).take(l as usize)).collect();
            if i % 2 == 1 {
                read = crate::evidence::revcomp(&read);
            }
            fq.push_str(&format!("@r{i}\n{}\n+\n{}\n", String::from_utf8(read.clone()).unwrap(), "I".repeat(read.len())));
        }
        let rq = dir.join("r.fastq");
        fs::write(&rq, fq).unwrap();
        let out = dir.join("p.fasta");
        let rep = polish(&fa, &[rq], &out, MIN_OBSERVATIONS).unwrap();
        let text = fs::read_to_string(&out).unwrap();
        let seq: Vec<u8> = text.lines().nth(1).unwrap().bytes().collect();
        assert!(text.starts_with(&format!(">m circular=true length={}", truth.len())), "{}", &text[..40]);
        assert_eq!(seq, truth, "every run back to its true length");
        assert_eq!(rep.molecules[0].runs_shortened, 0);
        assert!(rep.molecules[0].runs_lengthened > 0);
        fs::remove_dir_all(&dir).unwrap();
    }
}
