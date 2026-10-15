//! Independent Rust homopolymer-coordinate assembly and run-length restoration.
//!
//! Uses OVASM's Rust sparse graph, not Syncasm/FFI. Exact compressed overlaps
//! identify shared graph positions. Original read positions contribute once per
//! shared coordinate, even when several overlapping anchors observe them. This
//! experimental backend does not implement read correction or repeat unzipping.

use crate::{
    assemble::{self, AssembleParams, AssembleReport},
    evidence::{parse_gfa, revcomp, Graph},
    recruit::{peak_rss_mb, read_batches},
};
use anyhow::{bail, Context, Result};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use std::{
    fmt::Write as _,
    fs,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};

fn compress(sequence: &[u8]) -> Result<(Vec<u8>, Vec<u32>)> {
    let mut bases = Vec::new();
    let mut lengths: Vec<u32> = Vec::new();
    for b in sequence.iter().map(u8::to_ascii_uppercase) {
        if matches!(b, b'A' | b'C' | b'G' | b'T') {
            if bases.last() == Some(&b) {
                let last = lengths.last_mut().unwrap();
                *last = last
                    .checked_add(1)
                    .context("homopolymer exceeds u32 length")?;
            } else {
                bases.push(b);
                lengths.push(1);
            }
        } else {
            // Every ambiguous base is a hard sampling boundary; never join the
            // ACGT runs either side of it or give it a consensus observation.
            bases.push(b'N');
            lengths.push(1);
        }
    }
    Ok((bases, lengths))
}

struct UnionFind(Vec<usize>);
impl UnionFind {
    fn find(&mut self, mut x: usize) -> usize {
        while self.0[x] != x {
            self.0[x] = self.0[self.0[x]];
            x = self.0[x];
        }
        x
    }
    fn union(&mut self, a: usize, b: usize) {
        let (a, b) = (self.find(a), self.find(b));
        // Stable representative; no orientation is needed for run lengths.
        if a < b {
            self.0[b] = a;
        } else {
            self.0[a] = b;
        }
    }
}

#[derive(Clone, Copy)]
struct Anchor {
    segment: usize,
    position: usize,
    forward: bool,
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Span {
    segment: usize,
    forward: bool,
    diagonal: i64,
    start: usize,
    end: usize,
}
impl Span {
    fn same_diagonal(&self, other: &Self) -> bool {
        (self.segment, self.forward, self.diagonal)
            == (other.segment, other.forward, other.diagonal)
    }
    fn graph_position(&self, q: usize) -> usize {
        (if self.forward {
            self.diagonal + q as i64
        } else {
            self.diagonal - q as i64
        }) as usize
    }
}

#[derive(Default, Serialize)]
pub struct ConsensusStats {
    pub reads: u64,
    pub bases: u64,
    pub exact_anchor_placements: u64,
    /// A read word can legitimately occur at several repeated graph positions.
    /// These placements are not evidence of a unique locus or physical origin.
    pub multiple_placement_words: u64,
    pub rejected_hash_collisions: u64,
    pub distinct_observations: u64,
    pub duplicate_anchor_observations: u64,
    pub censored_terminal_runs: u64,
    pub graph_run_positions: usize,
    pub shared_run_positions: usize,
    pub minimum_observations: u64,
    pub unknown_run_positions: usize,
}

struct Consensus {
    lengths: Vec<Vec<usize>>,
    stats: ConsensusStats,
}

fn restore(graph: &Graph, inputs: &[PathBuf], p: &AssembleParams) -> Result<Consensus> {
    let mut offsets = vec![0];
    for seq in &graph.seqs {
        offsets.push(offsets.last().unwrap() + seq.len());
    }
    let positions = *offsets.last().unwrap();
    let mut uf = UnionFind((0..positions).collect());
    for (i, &(a, b)) in graph.links.iter().enumerate() {
        let ov = graph.overlaps[i];
        let left = graph.oriented_seq(a);
        let right = graph.oriented_seq(b);
        if left[left.len() - ov..] != right[..ov] {
            bail!("inconsistent HPC overlap at link {i}");
        }
        for j in 0..ov {
            let ai = if a.1 { left.len() - ov + j } else { ov - 1 - j };
            let bi = if b.1 { j } else { right.len() - 1 - j };
            uf.union(offsets[a.0] + ai, offsets[b.0] + bi);
        }
    }
    let roots: Vec<_> = (0..positions).map(|i| uf.find(i)).collect();
    let mut index: FxHashMap<u64, Vec<Anchor>> = FxHashMap::default();
    for (segment, seq) in graph.seqs.iter().enumerate() {
        for (hash, forward, position) in assemble::sample_read(seq, p) {
            index.entry(hash).or_default().push(Anchor {
                segment,
                forward,
                position: position as usize,
            });
        }
    }
    let mut sums = vec![0u64; positions];
    let mut counts = vec![0u64; positions];
    let mut stats = ConsensusStats {
        graph_run_positions: positions,
        shared_run_positions: positions - roots.iter().copied().collect::<FxHashSet<_>>().len(),
        ..Default::default()
    };
    let (rx, reader) = read_batches(inputs);
    for batch in rx {
        for read in batch {
            stats.reads += 1;
            stats.bases += read.seq.len() as u64;
            let (seq, runs) = compress(&read.seq)?;
            let rc = revcomp(&seq);
            let mut spans = Vec::new();
            for (hash, forward, position) in assemble::sample_read(&seq, p) {
                let q = position as usize;
                let mut matches = 0;
                for a in index.get(&hash).map(Vec::as_slice).unwrap_or(&[]) {
                    let same = forward == a.forward;
                    let word = if same {
                        &seq[q..q + p.k]
                    } else {
                        &rc[seq.len() - q - p.k..seq.len() - q]
                    };
                    if word != &graph.seqs[a.segment][a.position..a.position + p.k] {
                        stats.rejected_hash_collisions += 1;
                        continue;
                    }
                    matches += 1;
                    spans.push(Span {
                        segment: a.segment,
                        forward: same,
                        diagonal: if same {
                            a.position as i64 - q as i64
                        } else {
                            (a.position + p.k - 1 + q) as i64
                        },
                        start: q,
                        end: q + p.k,
                    });
                }
                stats.exact_anchor_placements += matches;
                stats.multiple_placement_words += u64::from(matches > 1);
            }
            // Merge collinear overlapping anchor intervals first. Counting every
            // k-long anchor independently would overweight densely sampled runs.
            spans.sort_unstable();
            let mut merged: Vec<Span> = Vec::new();
            for s in spans {
                if let Some(last) = merged.last_mut() {
                    if last.same_diagonal(&s) && s.start <= last.end {
                        stats.duplicate_anchor_observations +=
                            last.end.min(s.end).saturating_sub(s.start) as u64;
                        last.end = last.end.max(s.end);
                        continue;
                    }
                }
                merged.push(s);
            }
            let mut observed = FxHashSet::default();
            for mut s in merged {
                // An exact anchor certifies a diagonal, not just its k bases.
                // Extend it to mismatches/read/segment boundaries. Arbitrary
                // GFA ends need not themselves be closed-syncmer starts.
                let matches = |q: usize| {
                    let gp = if s.forward {
                        s.diagonal + q as i64
                    } else {
                        s.diagonal - q as i64
                    };
                    let target = &graph.seqs[s.segment];
                    if gp < 0
                        || gp as usize >= target.len()
                        || !matches!(seq[q], b'A' | b'C' | b'G' | b'T')
                    {
                        return false;
                    }
                    let b = target[gp as usize];
                    seq[q]
                        == if s.forward {
                            b
                        } else {
                            match b {
                                b'A' => b'T',
                                b'C' => b'G',
                                b'G' => b'C',
                                b'T' => b'A',
                                _ => b'N',
                            }
                        }
                };
                while s.start > 0 && matches(s.start - 1) {
                    s.start -= 1;
                }
                while s.end < seq.len() && matches(s.end) {
                    s.end += 1;
                }
                for q in s.start..s.end {
                    // The first and last run may be truncated by the read ends.
                    if q == 0 || q + 1 == seq.len() {
                        stats.censored_terminal_runs += 1;
                        continue;
                    }
                    let root = roots[offsets[s.segment] + s.graph_position(q)];
                    if observed.insert((root, q)) {
                        sums[root] = sums[root]
                            .checked_add(runs[q] as u64)
                            .context("run sum overflow")?;
                        counts[root] += 1;
                        stats.distinct_observations += 1;
                    } else {
                        stats.duplicate_anchor_observations += 1;
                    }
                }
            }
        }
    }
    reader
        .join()
        .map_err(|_| anyhow::anyhow!("read parser panicked"))??;
    stats.minimum_observations = roots.iter().map(|&r| counts[r]).min().unwrap_or(0);
    stats.unknown_run_positions = roots.iter().filter(|&&r| counts[r] == 0).count();
    if stats.unknown_run_positions > 0 {
        bail!(
            "{} graph run positions lack an uncensored original observation; no base graph emitted",
            stats.unknown_run_positions
        );
    }
    let lengths = graph
        .seqs
        .iter()
        .enumerate()
        .map(|(s, seq)| {
            (0..seq.len())
                .map(|i| {
                    let r = roots[offsets[s] + i];
                    ((sums[r] + counts[r] / 2) / counts[r]) as usize
                })
                .collect()
        })
        .collect();
    Ok(Consensus { lengths, stats })
}

fn base_graph(graph: &Graph, consensus: &Consensus) -> Result<Graph> {
    let seqs: Vec<Vec<u8>> = graph
        .seqs
        .iter()
        .zip(&consensus.lengths)
        .map(|(seq, lens)| {
            seq.iter()
                .zip(lens)
                .flat_map(|(&b, &n)| std::iter::repeat(b).take(n))
                .collect()
        })
        .collect();
    let mut result = Graph {
        names: graph.names.clone(),
        seqs,
        links: graph.links.clone(),
        overlaps: Vec::new(),
    };
    for (i, &(a, b)) in graph.links.iter().enumerate() {
        let lens = &consensus.lengths[a.0];
        let n = graph.overlaps[i];
        let ov = if a.1 {
            lens[lens.len() - n..].iter().sum()
        } else {
            lens[..n].iter().sum()
        };
        let left = result.oriented_seq(a);
        let right = result.oriented_seq(b);
        if ov >= left.len() || ov >= right.len() || left[left.len() - ov..] != right[..ov] {
            bail!("restored-base overlap disagrees at link {i}; no graph emitted");
        }
        result.overlaps.push(ov);
    }
    Ok(result)
}

#[derive(Serialize)]
pub struct HpcReport {
    pub graph_engine: &'static str,
    pub experimental: bool,
    pub uses_syncasm: bool,
    pub homopolymer_compressed_k: bool,
    pub read_error_correction: bool,
    pub read_threaded_unzipping: bool,
    pub physical_topology_identified: bool,
    pub reads: u64,
    pub bases: u64,
    pub compressed_bases: u64,
    pub k: usize,
    pub s: usize,
    pub min_count: u32,
    pub unitigs: usize,
    pub links: usize,
    pub total_length: usize,
    pub consensus_estimator: &'static str,
    pub consensus: ConsensusStats,
    pub compressed_backbone: AssembleReport,
    pub intermediates: PathBuf,
    pub elapsed_seconds: f64,
    pub peak_rss_mb: Option<f64>,
}

pub fn run(
    inputs: &[PathBuf],
    k: usize,
    s: usize,
    min_count: u32,
    out: &Path,
    report: &Path,
) -> Result<HpcReport> {
    if inputs.is_empty() || k < 3 || k % 2 == 0 || s < 5 || s > 31 || s >= k || min_count == 0 {
        bail!("Rust HPC requires reads, odd k > s, s in 5..31 and positive --min-count");
    }
    if out.exists() || report.exists() || out == report {
        bail!("Rust HPC requires distinct new graph/report output paths");
    }
    let started = Instant::now();
    let fingerprints: Vec<_> = inputs
        .iter()
        .map(|p| -> Result<_> {
            let m = fs::metadata(p)?;
            Ok((m.len(), m.modified()?))
        })
        .collect::<Result<_>>()?;
    let work = out.with_extension("rust-hpc");
    if let Some(parent) = work.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::create_dir(&work)
        .context("Rust HPC intermediates already exist; do not reuse a failed run")?;
    let compressed = work.join("reads.hpc.fasta");
    let mut writer = BufWriter::new(fs::File::create(&compressed)?);
    let (mut reads, mut bases, mut compressed_bases) = (0u64, 0u64, 0u64);
    let (rx, reader) = read_batches(inputs);
    for batch in rx {
        for r in batch {
            let (seq, _) = compress(&r.seq)?;
            writeln!(writer, ">r{reads}")?;
            writer.write_all(&seq)?;
            writer.write_all(b"\n")?;
            reads += 1;
            bases += r.seq.len() as u64;
            compressed_bases += seq.len() as u64;
        }
    }
    reader
        .join()
        .map_err(|_| anyhow::anyhow!("read parser panicked"))??;
    writer.flush()?;
    drop(writer);
    if reads == 0 {
        bail!("no reads");
    }
    let p = AssembleParams {
        k,
        s,
        min_count: Some(min_count),
        min_count_fraction: None,
    };
    let hpc_graph = work.join("compressed.gfa");
    let backbone = assemble::run_backbone_layout(
        vec![compressed],
        p,
        &hpc_graph,
        &work.join("compressed.json"),
        true,
    )?;
    let graph = parse_gfa(&hpc_graph)?;
    if graph.seqs.is_empty() {
        bail!("Rust HPC retained no graph; no base graph emitted");
    }
    let consensus = restore(&graph, inputs, &p)?;
    if (reads, bases) != (consensus.stats.reads, consensus.stats.bases) {
        bail!("input statistics changed between compression and consensus");
    }
    for (path, before) in inputs.iter().zip(fingerprints) {
        let m = fs::metadata(path)?;
        if before != (m.len(), m.modified()?) {
            bail!("input changed during assembly: {}", path.display());
        }
    }
    let restored = base_graph(&graph, &consensus)?;
    let mut text = String::from("H\tVN:Z:1.0\tPG:Z:ovasm-rust-hpc-experimental\n");
    let source = fs::read_to_string(&hpc_graph)?;
    let tags: Vec<_> = source
        .lines()
        .filter(|l| l.starts_with("S\t"))
        .map(|l| {
            l.split('\t')
                .skip(3)
                .filter(|x| !x.starts_with("LN:i:"))
                .collect::<Vec<_>>()
                .join("\t")
        })
        .collect();
    for (i, seq) in restored.seqs.iter().enumerate() {
        writeln!(
            text,
            "S\t{}\t{}\tLN:i:{}\tHC:i:{}\t{}",
            restored.names[i],
            std::str::from_utf8(seq)?,
            seq.len(),
            graph.seqs[i].len(),
            tags[i]
        )?;
    }
    for (i, &(a, b)) in restored.links.iter().enumerate() {
        writeln!(
            text,
            "L\t{}\t{}\t{}\t{}\t{}M",
            restored.names[a.0],
            if a.1 { '+' } else { '-' },
            restored.names[b.0],
            if b.1 { '+' } else { '-' },
            restored.overlaps[i]
        )?;
    }
    let r = HpcReport { graph_engine: "ovasm-rust-hpc", experimental: true, uses_syncasm: false,
        homopolymer_compressed_k: true, read_error_correction: false, read_threaded_unzipping: false,
        physical_topology_identified: false, reads, bases, compressed_bases, k, s, min_count,
        unitigs: restored.seqs.len(), links: restored.links.len(), total_length: restored.seqs.iter().map(Vec::len).sum(),
        consensus_estimator: "rounded mean of distinct original read/run observations per shared graph position; terminal runs censored; multiple placements retained",
        consensus: consensus.stats, compressed_backbone: backbone, intermediates: work,
        elapsed_seconds: started.elapsed().as_secs_f64(), peak_rss_mb: peak_rss_mb() };
    if let Some(parent) = report.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(out, text)?;
    fs::write(report, serde_json::to_string_pretty(&r)?)?;
    Ok(r)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};
    static ID: AtomicUsize = AtomicUsize::new(0);
    fn dir() -> PathBuf {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join(".agent-work/rust-hpc/unit-tests");
        fs::create_dir_all(&root).unwrap();
        let p = root.join(format!(
            "ovasm-rust-hpc-{}-{}",
            std::process::id(),
            ID.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&p).unwrap();
        assert_eq!(
            p.parent().unwrap().canonicalize().unwrap(),
            root.canonicalize().unwrap()
        );
        p
    }
    fn random_runs(n: usize) -> Vec<u8> {
        random_runs_seed(n, 311)
    }
    fn random_runs_seed(n: usize, mut state: u64) -> Vec<u8> {
        let mut seq = Vec::new();
        for _ in 0..n {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let b = b"ACGT"[state as usize % 4];
            seq.push(if seq.last() == Some(&b) {
                b"ACGT"[(state as usize % 4 + 1) % 4]
            } else {
                b
            });
        }
        seq
    }
    #[test]
    fn compression_preserves_ambiguous_boundaries_and_long_runs() {
        let mut seq = vec![b'a'; 1000];
        seq.extend_from_slice(b"NnAAActTT?");
        let (bases, runs) = compress(&seq).unwrap();
        assert_eq!(
            bases,
            b"ANN ACTN"
                .iter()
                .copied()
                .filter(|&b| b != b' ')
                .collect::<Vec<_>>()
        );
        assert_eq!(runs, vec![1000, 1, 1, 3, 1, 3, 1]);
        assert_eq!(runs.iter().map(|&n| n as usize).sum::<usize>(), seq.len());
    }
    #[test]
    fn reverse_complement_has_reverse_run_lengths() {
        let seq = b"AAAACCTNNGGGATTTTT";
        let (a, n) = compress(seq).unwrap();
        let (b, m) = compress(&revcomp(seq)).unwrap();
        assert_eq!(b, revcomp(&a));
        assert_eq!(m, n.into_iter().rev().collect::<Vec<_>>());
    }
    #[test]
    fn invalid_backend_combinations_and_missing_parameters_fail() {
        use clap::Parser;
        for flag in [
            "--multichrom-graph",
            "--auto-k-depth",
            "--min-count-fraction",
            "--read-context",
        ] {
            let mut args = vec![
                "ovasm",
                "assemble",
                "--reads",
                "missing",
                "--rust-hpc",
                "--out",
                "g.gfa",
                "--out-json",
                "g.json",
                flag,
            ];
            if flag != "--multichrom-graph" {
                args.push("1");
            }
            assert!(crate::Cli::try_parse_from(args).is_err());
        }
        let cli = crate::Cli::try_parse_from([
            "ovasm",
            "assemble",
            "--reads",
            "missing",
            "--rust-hpc",
            "--out",
            "g.gfa",
            "--out-json",
            "g.json",
        ])
        .unwrap();
        assert!(crate::dispatch(cli.command)
            .unwrap_err()
            .to_string()
            .contains("explicit --k"));
        let d = dir();
        let input = d.join("reads.fa");
        fs::write(&input, ">x\nACGT\n").unwrap();
        #[cfg(not(syncasm_native))]
        {
            let error = crate::syncasm::run(
                vec![input],
                101,
                11,
                1,
                1,
                &d.join("g.gfa"),
                &d.join("g.json"),
            );
            assert!(error
                .err()
                .unwrap()
                .to_string()
                .contains("does not include"));
        }
        #[cfg(syncasm_native)]
        let _ = input;
        fs::remove_dir_all(d).unwrap();
    }
    #[test]
    fn shared_reverse_overlap_uses_distinct_observations() {
        let d = dir();
        let seq = random_runs(500);
        let graph = Graph {
            names: vec!["a".into(), "b".into()],
            seqs: vec![seq[..350].to_vec(), revcomp(&seq[250..])],
            links: vec![((0, true), (1, false))],
            overlaps: vec![100],
        };
        let expand = |n: usize| {
            seq.iter()
                .flat_map(|&b| std::iter::repeat(b).take(n))
                .collect::<Vec<_>>()
        };
        let input = d.join("reads.fa");
        fs::write(
            &input,
            format!(
                ">f\n{}\n>r\n{}\n",
                String::from_utf8(expand(2)).unwrap(),
                String::from_utf8(revcomp(&expand(4))).unwrap()
            ),
        )
        .unwrap();
        let p = AssembleParams {
            k: 51,
            s: 11,
            min_count: Some(1),
            min_count_fraction: None,
        };
        // A read-end-censored graph coordinate must fail instead of inventing a length.
        assert!(restore(&graph, &[input.clone()], &p).is_err());
        let front = [if seq[0] == b'A' { b'C' } else { b'A' }];
        let back = [if *seq.last().unwrap() == b'A' {
            b'C'
        } else {
            b'A'
        }];
        let padded = [front.as_slice(), expand(2).as_slice(), back.as_slice()].concat();
        let padded_rc =
            revcomp(&[front.as_slice(), expand(4).as_slice(), back.as_slice()].concat());
        fs::write(
            &input,
            format!(
                ">f\n{}\n>r\n{}\n",
                String::from_utf8(padded).unwrap(),
                String::from_utf8(padded_rc).unwrap()
            ),
        )
        .unwrap();
        let c = restore(&graph, &[input], &p).unwrap();
        assert!(c.lengths.iter().flatten().all(|&n| n == 3));
        assert_eq!(c.stats.distinct_observations, 1000);
        assert_eq!(c.stats.minimum_observations, 2);
        assert_eq!(c.stats.shared_run_positions, 100);
        let base = base_graph(&graph, &c).unwrap();
        assert_eq!(base.overlaps, vec![300]);
        fs::remove_dir_all(d).unwrap();
    }
    #[test]
    fn rust_graph_recovers_two_independent_circles() {
        let d = dir();
        let a = random_runs(1800);
        let b = revcomp(&random_runs_seed(2700, 557));
        let mut fasta = String::new();
        for (c, circle) in [&a, &b].into_iter().enumerate() {
            let expanded: Vec<_> = circle.iter().flat_map(|&x| [x, x]).collect();
            let doubled = expanded.repeat(3);
            for i in 0..16 {
                let start = i * 137 % expanded.len();
                let seq = &doubled[start..start + expanded.len() + 500];
                let seq = if i % 2 == 0 {
                    seq.to_vec()
                } else {
                    revcomp(seq)
                };
                writeln!(fasta, ">c{c}_{i}\n{}", String::from_utf8(seq).unwrap()).unwrap();
            }
        }
        let input = d.join("reads.fa");
        fs::write(&input, fasta).unwrap();
        let out = d.join("graph.gfa");
        let result = run(&[input], 101, 15, 3, &out, &d.join("report.json")).unwrap();
        assert_eq!((result.reads, result.unitigs, result.links), (32, 2, 2));
        assert!(!result.uses_syncasm && !result.read_threaded_unzipping);
        assert_eq!(result.consensus.unknown_run_positions, 0);
        let graph = parse_gfa(&out).unwrap();
        for circle in [&a, &b] {
            let expanded: Vec<_> = circle.iter().flat_map(|&x| [x, x]).collect();
            assert!(graph.seqs.iter().enumerate().any(|(i, seq)| {
                let ov = graph
                    .links
                    .iter()
                    .enumerate()
                    .find(|(_, &(x, y))| x.0 == i && x == y)
                    .map(|(j, _)| graph.overlaps[j])
                    .unwrap_or(0);
                let trimmed = &seq[..seq.len() - ov];
                trimmed.len() == expanded.len()
                    && [trimmed.to_vec(), revcomp(trimmed)]
                        .iter()
                        .any(|s| expanded.repeat(2).windows(s.len()).any(|w| w == s))
            }));
        }
        fs::remove_dir_all(d).unwrap();
    }
}
