//! 2-bit canonical k-mers over DNA byte slices.

/// Largest k that fits a 2-bit packed k-mer in a `u64`.
pub const MAX_K: usize = 31;

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

/// Iterator over `(start, canonical_kmer)`; windows containing a non-ACGT base are skipped.
///
/// The canonical form is the smaller of the forward and reverse-complement encodings,
/// so a read and its reverse complement produce the same k-mer set.
pub struct CanonicalKmers<'a> {
    seq: &'a [u8],
    k: usize,
    mask: u64,
    shift: u32,
    fwd: u64,
    rev: u64,
    valid: usize,
    pos: usize,
}

impl<'a> CanonicalKmers<'a> {
    pub fn new(seq: &'a [u8], k: usize) -> Self {
        assert!((1..=MAX_K).contains(&k), "k must be in 1..={MAX_K}");
        Self {
            seq,
            k,
            mask: (1u64 << (2 * k)) - 1,
            shift: (2 * (k - 1)) as u32,
            fwd: 0,
            rev: 0,
            valid: 0,
            pos: 0,
        }
    }
}

impl Iterator for CanonicalKmers<'_> {
    type Item = (usize, u64);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        while self.pos < self.seq.len() {
            let base = self.seq[self.pos];
            self.pos += 1;
            match encode(base) {
                Some(code) => {
                    self.fwd = ((self.fwd << 2) | code) & self.mask;
                    self.rev = (self.rev >> 2) | ((3 - code) << self.shift);
                    self.valid += 1;
                    if self.valid >= self.k {
                        return Some((self.pos - self.k, self.fwd.min(self.rev)));
                    }
                }
                None => {
                    self.valid = 0;
                    self.fwd = 0;
                    self.rev = 0;
                }
            }
        }
        None
    }
}

/// The canonical k-mers of a 2-bit packed sequence (base i at bits 2(i % 4) of byte i / 4, A C G
/// T = 0 1 2 3, as the recruitment cache stores ACGT-only reads), in order: the keys
/// `CanonicalKmers` gives for the unpacked sequence, without decoding a base or testing for N.
pub struct PackedKmers<'a> {
    packed: &'a [u8],
    len: usize,
    k: usize,
    mask: u64,
    shift: u32,
    fwd: u64,
    rev: u64,
    pos: usize,
}

impl<'a> PackedKmers<'a> {
    pub fn new(packed: &'a [u8], len: usize, k: usize) -> Self {
        assert!((1..=MAX_K).contains(&k), "k must be in 1..={MAX_K}");
        assert!(packed.len() * 4 >= len, "packed sequence shorter than its length");
        Self {
            packed,
            len,
            k,
            mask: (1u64 << (2 * k)) - 1,
            shift: (2 * (k - 1)) as u32,
            fwd: 0,
            rev: 0,
            pos: 0,
        }
    }

    #[inline]
    fn push(&mut self) {
        let code = u64::from(self.packed[self.pos >> 2] >> (2 * (self.pos & 3))) & 3;
        self.fwd = ((self.fwd << 2) | code) & self.mask;
        self.rev = (self.rev >> 2) | ((3 - code) << self.shift);
        self.pos += 1;
    }
}

impl Iterator for PackedKmers<'_> {
    type Item = u64;

    #[inline]
    fn next(&mut self) -> Option<u64> {
        if self.len < self.k {
            return None;
        }
        while self.pos + 1 < self.k {
            self.push();
        }
        if self.pos >= self.len {
            return None;
        }
        self.push();
        Some(self.fwd.min(self.rev))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packed_kmers_are_the_canonical_kmers_of_the_unpacked_sequence() {
        for len in 0..70usize {
            let seq: Vec<u8> = (0..len).map(|i| b"ACGT"[(i * 7 + i / 3 + len) % 4]).collect();
            let mut packed = vec![0u8; len.div_ceil(4)];
            for (i, &b) in seq.iter().enumerate() {
                packed[i / 4] |= (encode(b).unwrap() as u8) << (2 * (i % 4));
            }
            for k in [1, 5, 21, 31] {
                let want: Vec<u64> = CanonicalKmers::new(&seq, k).map(|(_, km)| km).collect();
                let got: Vec<u64> = PackedKmers::new(&packed, len, k).collect();
                assert_eq!(got, want, "len {len} k {k}");
            }
        }
    }

    fn revcomp(s: &[u8]) -> Vec<u8> {
        s.iter()
            .rev()
            .map(|b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                x => *x,
            })
            .collect()
    }

    #[test]
    fn reverse_complement_gives_same_kmer_set() {
        let s = b"ACGTTGCAAGGCTTACGATCGATCGGA";
        let mut a: Vec<u64> = CanonicalKmers::new(s, 7).map(|(_, k)| k).collect();
        let rc = revcomp(s);
        let mut b: Vec<u64> = CanonicalKmers::new(&rc, 7).map(|(_, k)| k).collect();
        a.sort_unstable();
        b.sort_unstable();
        assert_eq!(a, b);
    }

    #[test]
    fn windows_with_n_are_skipped_and_positions_are_starts() {
        let s = b"ACGTNACGTA";
        let got: Vec<usize> = CanonicalKmers::new(s, 4).map(|(p, _)| p).collect();
        assert_eq!(got, vec![0, 5, 6]);
    }

    #[test]
    fn lowercase_equals_uppercase() {
        let a: Vec<u64> = CanonicalKmers::new(b"acgtacgt", 5)
            .map(|(_, k)| k)
            .collect();
        let b: Vec<u64> = CanonicalKmers::new(b"ACGTACGT", 5)
            .map(|(_, k)| k)
            .collect();
        assert_eq!(a, b);
    }
}
