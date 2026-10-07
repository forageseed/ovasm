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

#[cfg(test)]
mod tests {
    use super::*;

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
