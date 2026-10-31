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
    lookup: Option<LookupIndex>,
    changes: Option<FxHashSet<u64>>,
}

/// Read-only index of the bait used during one scan: a blocked Bloom filter (three bits in one
/// 64-byte block per key) in front of an exact open-addressing copy of the key-to-mask map. A
/// clear bit proves absence; a set one is settled by the table, so masks are exactly the map's.
///
/// Classification is bound by memory latency: a 10 kb HiFi read looks up ten thousand k-mers at
/// random places of a bait of ~10 million keys. The previous single-bit filter passed 12% of
/// absent keys on to the hash map (one more cache miss each), and the map's buckets could not be
/// prefetched. Here one block holds all of a key's bits (about 0.6% pass, 12 bits a key), the
/// table's slots can be prefetched (`BaitSet::scan_masks`) and a hit is a few probes in one line.
#[derive(Clone)]
struct LookupIndex {
    bloom: Vec<BloomBlock>,
    bloom_mask: usize,
    /// `EMPTY` for an unused slot; linear probing.
    keys: Vec<u64>,
    masks: Vec<u8>,
    table_mask: usize,
}

#[derive(Clone, Copy, Default)]
#[repr(align(64))]
struct BloomBlock([u64; 8]);

/// No canonical k-mer (at most 62 bits, `MAX_K` 31) is all ones.
const EMPTY: u64 = u64::MAX;
const BLOOM_BITS_PER_KEY: usize = 12;

/// Half of MurmurHash3's finalizer: k-mer codes are far from uniform in their low bits, and one
/// multiply keeps a lookup of a small (cache-resident) set as cheap as the bit test it is.
#[inline]
fn mix(key: u64) -> u64 {
    let h = (key ^ (key >> 33)).wrapping_mul(0xFF51_AFD7_ED55_8CCD);
    h ^ (h >> 33)
}

/// An index this small stays in cache: `scan_masks` looks it up without prefetching.
const PREFETCH_FROM_BYTES: usize = 4 << 20;

#[inline]
fn prefetch<T>(r: &T) {
    #[cfg(target_arch = "x86_64")]
    // SAFETY: a prefetch is only a hint; `r` is a valid reference.
    unsafe {
        core::arch::x86_64::_mm_prefetch(r as *const T as *const i8, core::arch::x86_64::_MM_HINT_T0);
    }
    #[cfg(not(target_arch = "x86_64"))]
    let _ = r;
}

impl LookupIndex {
    fn build(map: &FxHashMap<u64, u8>) -> Self {
        let n = map.len();
        let blocks = (n * BLOOM_BITS_PER_KEY).div_ceil(512).max(1).next_power_of_two();
        // load at most 0.8: hits take a few probes; misses come here only past the filter
        let slots = (n * 5 / 4 + 1).max(16).next_power_of_two();
        let mut index = Self {
            bloom: vec![BloomBlock::default(); blocks],
            bloom_mask: blocks - 1,
            keys: vec![EMPTY; slots],
            masks: vec![0; slots],
            table_mask: slots - 1,
        };
        for (&key, &mask) in map {
            let h = mix(key);
            let block = &mut index.bloom[h as usize & index.bloom_mask];
            for bit in Self::bits(h) {
                block.0[bit >> 6] |= 1u64 << (bit & 63);
            }
            let mut slot = Self::slot(h, index.table_mask);
            while index.keys[slot] != EMPTY {
                slot = (slot + 1) & index.table_mask;
            }
            index.keys[slot] = key;
            index.masks[slot] = mask;
        }
        index
    }

    /// Three bit positions in a 512-bit block, from bits the block index does not use.
    #[inline]
    fn bits(h: u64) -> [usize; 3] {
        [
            (h >> 37) as usize & 511,
            (h >> 46) as usize & 511,
            (h >> 55) as usize & 511,
        ]
    }

    #[inline]
    fn slot(h: u64, table_mask: usize) -> usize {
        h.rotate_left(32) as usize & table_mask
    }

    #[inline]
    fn may_contain(&self, h: u64) -> bool {
        let block = &self.bloom[h as usize & self.bloom_mask];
        Self::bits(h)
            .iter()
            .all(|&bit| block.0[bit >> 6] & (1u64 << (bit & 63)) != 0)
    }

    #[inline]
    fn exact(&self, key: u64, h: u64) -> u8 {
        let mut slot = Self::slot(h, self.table_mask);
        loop {
            let k = self.keys[slot];
            if k == key {
                return self.masks[slot];
            }
            if k == EMPTY {
                return 0;
            }
            slot = (slot + 1) & self.table_mask;
        }
    }

    #[inline]
    fn get(&self, key: u64) -> u8 {
        let h = mix(key);
        if self.may_contain(h) {
            self.exact(key, h)
        } else {
            0
        }
    }
}

/// How far ahead `scan_masks` prefetches: filter blocks `AHEAD` k-mers ahead, table slots of the
/// keys that passed `AHEAD / 2` ahead.
const AHEAD: usize = 16;

impl BaitSet {
    pub fn new(k: usize) -> Self {
        Self {
            k,
            map: FxHashMap::default(),
            targets: Vec::new(),
            lookup: None,
            changes: None,
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
        self.lookup = None;
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
            let entry = self.map.entry(kmer).or_insert(0);
            if *entry & bit == 0 {
                if let Some(changes) = self.changes.as_mut() {
                    changes.insert(kmer);
                }
            }
            *entry |= bit;
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
        self.lookup = None;
        let bit = 1u8 << target;
        let mut added = 0;
        for (_, kmer) in CanonicalKmers::new(seq, self.k) {
            let entry = self.map.entry(kmer).or_insert(0);
            if *entry & bit == 0 {
                if let Some(changes) = self.changes.as_mut() {
                    changes.insert(kmer);
                }
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
        self.lookup = None;
        let bit = 1u8 << target;
        let mut added = 0;
        for kmer in solid_kmers(self.k, reads, min_count, 0.1) {
            let entry = self.map.entry(kmer).or_insert(0);
            if *entry & bit == 0 {
                if let Some(changes) = self.changes.as_mut() {
                    changes.insert(kmer);
                }
                *entry |= bit;
                added += 1;
            }
        }
        added
    }

    #[inline]
    pub fn get(&self, kmer: u64) -> u8 {
        match &self.lookup {
            Some(index) => index.get(kmer),
            None => self.map.get(&kmer).copied().unwrap_or(0),
        }
    }

    /// `f(start, mask)` for every canonical k-mer of `seq` in order, as
    /// `CanonicalKmers::new(seq, k).map(|(pos, km)| (pos, self.get(km)))` gives them, with the
    /// prepared index's memory reads issued ahead so that their latencies overlap.
    #[inline]
    pub fn scan_masks(&self, seq: &[u8], mut f: impl FnMut(usize, u8)) {
        let index = match &self.lookup {
            Some(index) if index.bloom.len() * 64 + index.keys.len() * 9 >= PREFETCH_FROM_BYTES => index,
            _ => {
                for (pos, km) in CanonicalKmers::new(seq, self.k) {
                    f(pos, self.get(km));
                }
                return;
            }
        };
        // ring of (start, key, hash, passed the filter); entry i is filter-tested when i + AHEAD/2
        // arrives and resolved when i + AHEAD does
        let mut ring = [(0usize, 0u64, 0u64, false); AHEAD];
        let mut n = 0usize;
        let test = |e: &mut (usize, u64, u64, bool)| {
            e.3 = index.may_contain(e.2);
            if e.3 {
                prefetch(&index.keys[LookupIndex::slot(e.2, index.table_mask)]);
                prefetch(&index.masks[LookupIndex::slot(e.2, index.table_mask)]);
            }
        };
        for (pos, key) in CanonicalKmers::new(seq, self.k) {
            if n >= AHEAD {
                let e = ring[n % AHEAD];
                f(e.0, if e.3 { index.exact(e.1, e.2) } else { 0 });
            }
            let h = mix(key);
            prefetch(&index.bloom[h as usize & index.bloom_mask]);
            ring[n % AHEAD] = (pos, key, h, false);
            if n >= AHEAD / 2 {
                test(&mut ring[(n - AHEAD / 2) % AHEAD]);
            }
            n += 1;
        }
        // the last AHEAD entries: those not yet filter-tested first, then all in order
        for i in n.saturating_sub(AHEAD / 2)..n {
            test(&mut ring[i % AHEAD]);
        }
        for i in n.saturating_sub(AHEAD)..n {
            let e = ring[i % AHEAD];
            f(e.0, if e.3 { index.exact(e.1, e.2) } else { 0 });
        }
    }

    /// Whether any of `keys` is in the set, as `keys.any(|km| self.get(km) != 0)` says, with the
    /// lookups of a large prepared index prefetched like `scan_masks`. The replay passes ask this
    /// of every cached read against the k-mers the last extension changed.
    #[inline]
    pub fn any_key(&self, mut keys: impl Iterator<Item = u64>) -> bool {
        let index = match &self.lookup {
            Some(index) if index.bloom.len() * 64 + index.keys.len() * 9 >= PREFETCH_FROM_BYTES => index,
            _ => return keys.any(|km| self.get(km) != 0),
        };
        let mut ring = [(0u64, 0u64, false); AHEAD];
        let mut n = 0usize;
        let test = |e: &mut (u64, u64, bool)| {
            e.2 = index.may_contain(e.1);
            if e.2 {
                prefetch(&index.keys[LookupIndex::slot(e.1, index.table_mask)]);
            }
        };
        for key in keys.by_ref() {
            if n >= AHEAD {
                let e = ring[n % AHEAD];
                if e.2 && index.exact(e.0, e.1) != 0 {
                    return true;
                }
            }
            let h = mix(key);
            prefetch(&index.bloom[h as usize & index.bloom_mask]);
            ring[n % AHEAD] = (key, h, false);
            if n >= AHEAD / 2 {
                test(&mut ring[(n - AHEAD / 2) % AHEAD]);
            }
            n += 1;
        }
        for i in n.saturating_sub(AHEAD / 2)..n {
            test(&mut ring[i % AHEAD]);
        }
        (n.saturating_sub(AHEAD)..n).any(|i| {
            let e = ring[i % AHEAD];
            e.2 && index.exact(e.0, e.1) != 0
        })
    }

    /// The prepared index as `gpu_delta.cu` reads it: filter words (eight per block), block mask,
    /// table keys, table mask. Only the change set's are used there, whose masks are all 1.
    pub fn lookup_parts(&self) -> Option<(&[u64], u64, &[u64], u64)> {
        let index = self.lookup.as_ref()?;
        // SAFETY: a `BloomBlock` is eight `u64`s (64 bytes, no padding).
        let words = unsafe {
            std::slice::from_raw_parts(index.bloom.as_ptr().cast::<u64>(), index.bloom.len() * 8)
        };
        Some((words, index.bloom_mask as u64, &index.keys, index.table_mask as u64))
    }

    /// Rebuild after extension, before scanning reads (see `LookupIndex`).
    /// This changes lookup cost only: masks and all classification scores are exact.
    pub fn prepare_lookup(&mut self) {
        self.lookup = Some(LookupIndex::build(&self.map));
    }

    /// Only replay needs mutation tracking. A new target bit on an existing key
    /// matters just as much as a new key: either can change a read's classification.
    pub fn track_changes(&mut self) {
        self.changes = Some(FxHashSet::default());
    }

    pub fn take_changes(&mut self) -> Self {
        let keys = std::mem::take(self.changes.as_mut().expect("tracking enabled"));
        let mut delta = Self::new(self.k);
        delta.map = keys.into_iter().map(|key| (key, 1)).collect();
        delta.prepare_lookup();
        delta
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
    // The last bin holds every k-mer seen in more than 4,096 reads, and it is the peak of a
    // target deeper than that: left out, the peak fell in the error tail and recurrent errors
    // joined the bait (Medicago polymorpha plastid at ~10,000x: 1.15 M k-mers added to a
    // 125 kb genome in one pass).
    let peak = (valley..=4097).max_by_key(|&c| hist[c]).unwrap_or(valley);
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

    #[test]
    fn tracking_includes_existing_keys_that_gain_another_target() {
        let seq = genome(1000, 401);
        let extra = genome(1000, 402);
        let mut bait = BaitSet::new(21);
        bait.target("a").unwrap();
        bait.target("b").unwrap();
        bait.add_sequence(0, &seq);
        let before = bait.clone();
        bait.track_changes();
        bait.add_sequence(1, &seq);
        bait.extend_from_reads(0, [extra.as_slice(), extra.as_slice(), extra.as_slice()], 3);
        let delta = bait.take_changes();
        for (&key, &mask) in &bait.map {
            assert_eq!(delta.get(key) != 0, before.get(key) != mask);
        }
        assert_eq!(bait.take_changes().len(), 0);
        assert!(CanonicalKmers::new(&seq, 21).all(|(_, key)| delta.get(key) != 0));
    }

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
    fn prepared_lookup_preserves_masks_and_extension_invalidates_it() {
        let mut bait = BaitSet::new(21);
        for t in 0..MAX_TARGETS {
            let target = bait.target(&format!("t{t}")).unwrap();
            bait.add_sequence(target, &genome(5000, t as u64 + 1));
        }
        let shared = genome(1000, 81);
        for t in 0..MAX_TARGETS {
            bait.add_sequence(t, &shared);
        }
        let exact = bait.clone();
        bait.prepare_lookup();
        for seq in [genome(100_000, 19), shared.clone(), genome(5000, 1)] {
            for (_, key) in CanonicalKmers::new(&seq, bait.k) {
                assert_eq!(bait.get(key), exact.get(key));
            }
        }
        assert!(CanonicalKmers::new(&shared, bait.k).all(|(_, key)| bait.get(key) == u8::MAX));
        let added = genome(1000, 101);
        bait.add_sequence(0, &added);
        assert!(CanonicalKmers::new(&added, bait.k).all(|(_, key)| bait.get(key) & 1 != 0));
        bait.prepare_lookup();
        bait.extend_from_reads(1, [added.as_slice(), added.as_slice(), added.as_slice()], 3);
        assert!(CanonicalKmers::new(&added, bait.k).all(|(_, key)| bait.get(key) & 3 == 3));
    }

    #[test]
    fn scanned_masks_are_the_lookups_in_order() {
        // 450,000 keys: an index past `PREFETCH_FROM_BYTES`, scanned with prefetching
        let mut bait = BaitSet::new(21);
        for t in 0..3 {
            let target = bait.target(&format!("t{t}")).unwrap();
            bait.add_sequence(target, &genome(150_000, t as u64 + 31));
        }
        let mut reads = vec![Vec::new(), b"ACGT".to_vec(), genome(25, 2), genome(30, 3)];
        reads.push([genome(3000, 4), genome(20_000, 31)[500..4500].to_vec()].concat());
        let mut mixed = genome(20_000, 32)[..2000].to_vec();
        mixed.extend_from_slice(b"NNNNacgtRY");
        mixed.extend(genome(500, 5));
        reads.push(mixed);
        for prepared in [false, true] {
            if prepared {
                bait.prepare_lookup();
                let index = bait.lookup.as_ref().unwrap();
                assert!(index.bloom.len() * 64 + index.keys.len() * 9 >= PREFETCH_FROM_BYTES);
            }
            for r in &reads {
                let expected: Vec<(usize, u8)> = CanonicalKmers::new(r, 21)
                    .map(|(pos, km)| (pos, bait.map.get(&km).copied().unwrap_or(0)))
                    .collect();
                let mut got = Vec::new();
                bait.scan_masks(r, |pos, mask| got.push((pos, mask)));
                assert_eq!(got, expected);
            }
        }
        // absent keys: a filter false positive still reads 0 from the table
        let absent = genome(200_000, 77);
        let mut wrong = 0;
        bait.scan_masks(&absent, |pos, mask| {
            let km = CanonicalKmers::new(&absent[pos..pos + 21], 21).next().unwrap().1;
            wrong += usize::from(mask != bait.map.get(&km).copied().unwrap_or(0));
        });
        assert_eq!(wrong, 0);
    }

    #[test]
    fn any_key_is_any_lookup_at_every_position_of_the_hit() {
        let mut bait = BaitSet::new(21);
        let target = bait.target("t").unwrap();
        bait.add_sequence(target, &genome(450_000, 41));
        let present = genome(450_000, 41);
        let absent = genome(30_000, 78);
        for prepared in [false, true] {
            if prepared {
                bait.prepare_lookup();
                let index = bait.lookup.as_ref().unwrap();
                assert!(index.bloom.len() * 64 + index.keys.len() * 9 >= PREFETCH_FROM_BYTES);
            }
            let keys = |s: &[u8]| CanonicalKmers::new(s, 21).map(|(_, km)| km).collect::<Vec<_>>();
            assert!(!bait.any_key(std::iter::empty()));
            assert!(!bait.any_key(keys(&absent).into_iter()));
            // one present k-mer at the start, inside the prefetch ring, and at the very end
            for at in [0, 3, 9, 15, 16, 17, 200, 29_978] {
                let mut s = absent[..at].to_vec();
                s.extend_from_slice(&present[1000..1021]);
                s.extend_from_slice(&absent[at..at + 50.min(30_000 - at)]);
                let ks = keys(&s);
                assert_eq!(bait.any_key(ks.iter().copied()), ks.iter().any(|&km| bait.get(km) != 0), "at {at}");
                assert!(bait.any_key(ks.into_iter()), "at {at}");
            }
        }
    }

    #[test]
    fn errors_stay_out_of_a_target_deeper_than_the_histogram() {
        // 5,000 copies of a 2 kb genome; one error recurs in 20 of them
        let g = genome(2000, 11);
        let mut bad = g.clone();
        bad[1000] = if bad[1000] == b'A' { b'C' } else { b'A' };
        let mut reads: Vec<&[u8]> = vec![g.as_slice(); 5000];
        reads.extend(std::iter::repeat(bad.as_slice()).take(20));
        let solid: FxHashSet<u64> = solid_kmers(21, reads, 3, 0.1).into_iter().collect();
        let error: Vec<u64> = CanonicalKmers::new(&bad[980..1021], 21)
            .map(|(_, k)| k)
            .collect();
        assert!(error.iter().all(|k| !solid.contains(k)), "the error's k-mers");
        assert!(CanonicalKmers::new(&g, 21).all(|(_, k)| solid.contains(&k)));
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
