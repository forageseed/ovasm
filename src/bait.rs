//! Bait k-mer sets: canonical k-mer -> bitmask of organelle targets it came from.

use std::path::Path;

use anyhow::{bail, Context, Result};
use rustc_hash::{FxHashMap, FxHashSet};

use crate::kmer::CanonicalKmers;

/// At most 8 targets (one bit each in the `u8` mask).
pub const MAX_TARGETS: usize = 8;

#[derive(Debug, Clone, serde::Serialize)]
pub struct Target {
    pub name: String,
    /// Summed length of the seed sequences, used as the genome-size estimate for depth.
    pub seed_length: u64,
    pub seed_files: Vec<String>,
    /// Seeded from reads of the sample (`add_seed_reads`): complete, so never extended.
    pub from_reads: bool,
}

#[derive(Clone)]
pub struct BaitSet {
    pub k: usize,
    map: FxHashMap<u64, u8>,
    pub targets: Vec<Target>,
}

impl BaitSet {
    pub fn new(k: usize) -> Self {
        Self {
            k,
            map: FxHashMap::default(),
            targets: Vec::new(),
        }
    }

    /// Register a target name, returning its index (existing names are reused).
    pub fn target(&mut self, name: &str) -> Result<usize> {
        if let Some(i) = self.targets.iter().position(|t| t.name == name) {
            return Ok(i);
        }
        if self.targets.len() == MAX_TARGETS {
            bail!("at most {MAX_TARGETS} targets are supported");
        }
        self.targets.push(Target {
            name: name.to_string(),
            seed_length: 0,
            seed_files: Vec::new(),
            from_reads: false,
        });
        Ok(self.targets.len() - 1)
    }

    /// Add every canonical k-mer of every record in a FASTA/FASTQ(.gz) seed file.
    pub fn add_seed_file(&mut self, target: usize, path: &Path) -> Result<()> {
        let mut reader = needletail::parse_fastx_file(path)
            .with_context(|| format!("cannot read seed file {}", path.display()))?;
        while let Some(record) = reader.next() {
            let record =
                record.with_context(|| format!("malformed record in {}", path.display()))?;
            let seq = record.seq();
            self.targets[target].seed_length += seq.len() as u64;
            self.add_sequence(target, &seq);
        }
        self.targets[target]
            .seed_files
            .push(path.display().to_string());
        Ok(())
    }

    /// Seed a target from reads instead of a reference: the k-mers solid in them (see
    /// `solid_kmers`), e.g. corrected long reads of the same sample, to recruit its short reads.
    /// A reference seed from another accession lacks what the sample has and it does not
    /// (Nipponbare: the NC_011033 seed held 211 bp of a 1.4 kb stretch, which no short read was
    /// recruited for). The genome-size estimate is the number of solid k-mers.
    pub fn add_seed_reads(&mut self, target: usize, path: &Path, min_count: u32) -> Result<usize> {
        let mut reader = needletail::parse_fastx_file(path)
            .with_context(|| format!("cannot read seed reads {}", path.display()))?;
        let mut reads: Vec<Vec<u8>> = Vec::new();
        while let Some(record) = reader.next() {
            let record =
                record.with_context(|| format!("malformed record in {}", path.display()))?;
            reads.push(record.seq().into_owned());
        }
        let solid = solid_kmers(
            self.k,
            reads.iter().map(Vec::as_slice),
            min_count,
            SEED_PEAK_FRACTION,
        );
        let bit = 1u8 << target;
        for &kmer in &solid {
            *self.map.entry(kmer).or_insert(0) |= bit;
        }
        self.targets[target].seed_length += solid.len() as u64;
        self.targets[target]
            .seed_files
            .push(format!("{} (reads)", path.display()));
        self.targets[target].from_reads = true;
        Ok(solid.len())
    }

    pub fn add_sequence(&mut self, target: usize, seq: &[u8]) {
        self.add_sequence_counted(target, seq);
    }

    /// `add_sequence`, returning how many k-mers were new to the target.
    pub fn add_sequence_counted(&mut self, target: usize, seq: &[u8]) -> usize {
        let bit = 1u8 << target;
        let mut added = 0;
        for (_, kmer) in CanonicalKmers::new(seq, self.k) {
            let entry = self.map.entry(kmer).or_insert(0);
            if *entry & bit == 0 {
                *entry |= bit;
                added += 1;
            }
        }
        added
    }

    /// Add k-mers seen in at least `min_count` recruited reads of a target (iterative baiting).
    ///
    /// Counting per read (not per occurrence) suppresses sequencing-error k-mers, which appear
    /// in only one read, while real organelle k-mers recur across hundreds of reads.
    pub fn extend_from_reads<'a, I>(&mut self, target: usize, reads: I, min_count: u32) -> usize
    where
        I: IntoIterator<Item = &'a [u8]>,
    {
        let bit = 1u8 << target;
        let mut added = 0;
        for kmer in solid_kmers(self.k, reads, min_count, 0.1) {
            let entry = self.map.entry(kmer).or_insert(0);
            if *entry & bit == 0 {
                *entry |= bit;
                added += 1;
            }
        }
        added
    }

    #[inline]
    pub fn get(&self, kmer: u64) -> u8 {
        self.map.get(&kmer).copied().unwrap_or(0)
    }

    pub fn len(&self) -> usize {
        self.map.len()
    }
}

/// Seed reads keep k-mers from a quarter of the depth peak. Long reads recruited for an
/// organelle include NUMT-bearing nuclear reads (Nipponbare ONT: 67-205 kb reads), whose flanks
/// reach a tenth of the peak; as bait those flanks recruited 16,079 nuclear short reads that
/// assembled into separate tangles. From 0.1 to 0.25 of the peak, gold 21-mer recall went from
/// 99.65% to 99.64% and non-gold k-mers from 17,200 to 10,800.
pub const SEED_PEAK_FRACTION: f64 = 0.25;

/// Canonical k-mers seen in at least `min_count` reads, and from `peak_fraction` of the depth
/// peak (a tenth for iterative extension).
///
/// Counting per read (not per occurrence) suppresses sequencing-error k-mers, which appear in only
/// one read, while real organelle k-mers recur across hundreds of reads. The depth-relative floor:
/// noisy reads (ONT, CLR) at hundreds of x produce millions of error k-mers seen a few times each,
/// and random short k-mers also occur in the nuclear genome, so a fixed small count lets extension
/// run away into nuclear reads. Keep k-mers from a tenth of the depth peak (the peak past the error
/// spike).
pub fn solid_kmers<'a, I>(k: usize, reads: I, min_count: u32, peak_fraction: f64) -> Vec<u64>
where
    I: IntoIterator<Item = &'a [u8]>,
{
    let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
    for seq in reads {
        let mut seen: FxHashSet<u64> = FxHashSet::default();
        for (_, kmer) in CanonicalKmers::new(seq, k) {
            if seen.insert(kmer) {
                *counts.entry(kmer).or_insert(0) += 1;
            }
        }
    }
    let mut hist = vec![0u64; 4098];
    for &n in counts.values() {
        hist[(n as usize).min(4097)] += 1;
    }
    let valley = (3..4097)
        .find(|&c| hist[c] <= hist[c - 1] && hist[c] < hist[c + 1])
        .unwrap_or(3);
    let peak = (valley..4097).max_by_key(|&c| hist[c]).unwrap_or(valley);
    let threshold = min_count.max((peak_fraction * peak as f64) as u32);
    let mut out: Vec<u64> = counts
        .into_iter()
        .filter(|&(_, n)| n >= threshold)
        .map(|(kmer, _)| kmer)
        .collect();
    out.sort_unstable();
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn genome(n: usize, seed: u64) -> Vec<u8> {
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

    #[test]
    fn a_seed_from_reads_keeps_the_genome_and_drops_their_errors() {
        // 30 kb circular genome, 3 kb reads every 50 bp (60x), each with one error.
        let g = genome(30_000, 7);
        let circ = [g.clone(), g[..3000].to_vec()].concat();
        let mut reads: Vec<Vec<u8>> = (0..g.len())
            .step_by(50)
            .map(|st| circ[st..st + 3000].to_vec())
            .collect();
        for (i, r) in reads.iter_mut().enumerate() {
            let p = (i * 997) % r.len();
            r[p] = if r[p] == b'A' { b'C' } else { b'A' };
        }
        let dir = std::env::temp_dir().join(format!("ovasm-bait-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("reads.fa");
        let text: String = reads
            .iter()
            .enumerate()
            .map(|(i, r)| format!(">r{i}\n{}\n", std::str::from_utf8(r).unwrap()))
            .collect();
        std::fs::write(&fa, text).unwrap();

        let mut bait = BaitSet::new(25);
        let t = bait.target("mitochondrion").unwrap();
        let n = bait.add_seed_reads(t, &fa, 3).unwrap();
        let genome_kmers: FxHashSet<u64> = CanonicalKmers::new(&circ, 25).map(|(_, k)| k).collect();
        let found = genome_kmers
            .iter()
            .filter(|&&k| bait.get(k) & (1 << t) != 0)
            .count();
        assert_eq!(found, genome_kmers.len(), "every genome k-mer is bait");
        assert_eq!(n, genome_kmers.len(), "no error k-mer is");
        let size = bait.targets[t].seed_length as f64;
        assert!((size - 30_000.0).abs() < 100.0, "genome size {size}");
        std::fs::remove_dir_all(&dir).ok();
    }
}
