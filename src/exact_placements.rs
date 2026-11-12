//! Exact bounded occurrence scanning for the seed-source audit.
use std::collections::VecDeque;
use std::io::{self, BufRead, Read, Write};
use std::str::FromStr;

const MAGIC: &[u8] = b"OVXPOS1\n";
const NONE: u32 = u32::MAX;
const LIMIT: u8 = 2;
const BUFFER: usize = 64 * 1024;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Scope {
    Panel,
    Records,
}

impl Scope {
    fn label(self) -> &'static str {
        match self {
            Self::Panel => "panel",
            Self::Records => "records",
        }
    }
}

impl FromStr for Scope {
    type Err = io::Error;

    fn from_str(value: &str) -> io::Result<Self> {
        match value {
            "panel" => Ok(Self::Panel),
            "records" => Ok(Self::Records),
            _ => Err(invalid("exact placement scope must be panel or records")),
        }
    }
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

fn base_index(base: u8) -> Option<usize> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

struct Node {
    next: [u32; 4],
    suffix: u32,
    output: u32,
    pattern: u32,
}

impl Node {
    fn new() -> Self {
        Self {
            next: [0; 4],
            suffix: 0,
            output: 0,
            pattern: NONE,
        }
    }
}

struct Pattern {
    length: u32,
    last_record: u64,
    hits: u8,
}

struct Matcher {
    trie: Vec<Node>,
    patterns: Vec<Pattern>,
    unsaturated: usize,
}

impl Matcher {
    fn new(mut words: impl BufRead) -> io::Result<Self> {
        let mut trie = vec![Node::new()];
        let mut patterns = Vec::new();
        let mut word = Vec::new();
        while words.read_until(b'\n', &mut word)? != 0 {
            if word.last() == Some(&b'\n') {
                word.pop();
                if word.last() == Some(&b'\r') {
                    word.pop();
                }
            }
            let length = u32::try_from(word.len())
                .ok()
                .filter(|&n| n > 0 && n < NONE)
                .ok_or_else(|| invalid("empty or too long exact placement pattern"))?;
            let mut node = 0;
            for &base in &word {
                let code = base_index(base)
                    .ok_or_else(|| invalid("exact placement patterns require uppercase ACGT"))?;
                if trie[node].next[code] == 0 {
                    let child = u32::try_from(trie.len())
                        .map_err(|_| invalid("too many exact placement states"))?;
                    trie.try_reserve(1).map_err(io::Error::other)?;
                    trie[node].next[code] = child;
                    trie.push(Node::new());
                }
                node = trie[node].next[code] as usize;
            }
            if trie[node].pattern != NONE {
                return Err(invalid("duplicate literal exact placement pattern"));
            }
            trie[node].pattern = u32::try_from(patterns.len())
                .ok()
                .filter(|&i| i < NONE)
                .ok_or_else(|| invalid("too many exact placement patterns"))?;
            patterns.try_reserve(1).map_err(io::Error::other)?;
            patterns.push(Pattern {
                length,
                last_record: u64::MAX,
                hits: 0,
            });
            word.clear();
        }
        if patterns.is_empty() {
            return Err(invalid("no exact placement patterns"));
        }

        // Failure transitions also retain the first shorter pattern ending at a state.
        let mut queue = VecDeque::new();
        queue
            .try_reserve(trie.len().saturating_sub(1))
            .map_err(io::Error::other)?;
        queue.extend(trie[0].next.into_iter().filter(|&node| node != 0));
        while let Some(index) = queue.pop_front() {
            let index = index as usize;
            let suffix = trie[index].suffix as usize;
            trie[index].output = if trie[suffix].pattern == NONE {
                trie[suffix].output
            } else {
                suffix as u32
            };
            for code in 0..4 {
                let child = trie[index].next[code];
                let fallback = trie[suffix].next[code];
                if child == 0 {
                    trie[index].next[code] = fallback;
                } else {
                    trie[child as usize].suffix = fallback;
                    queue.push_back(child);
                }
            }
        }
        Ok(Self {
            trie,
            unsaturated: patterns.len(),
            patterns,
        })
    }

    fn scan_chunk(
        &mut self,
        scope: Scope,
        record: u64,
        offset: u64,
        sequence: &[u8],
        state: &mut usize,
        output: &mut impl Write,
    ) -> io::Result<()> {
        if scope == Scope::Panel && self.unsaturated == 0 {
            return Ok(());
        }
        for (index, &base) in sequence.iter().enumerate() {
            *state = base_index(base).map_or(0, |code| self.trie[*state].next[code] as usize);
            let mut node = *state;
            while node != 0 {
                let pattern_id = self.trie[node].pattern;
                node = self.trie[node].output as usize;
                if pattern_id == NONE {
                    continue;
                }
                let pattern = &mut self.patterns[pattern_id as usize];
                if scope == Scope::Records && pattern.last_record != record {
                    pattern.last_record = record;
                    pattern.hits = 0;
                }
                if pattern.hits == LIMIT {
                    continue;
                }
                let start = offset + index as u64 + 1 - u64::from(pattern.length);
                writeln!(output, "H\t{record}\t{pattern_id}\t{start}")?;
                pattern.hits += 1;
                if scope == Scope::Panel && pattern.hits == LIMIT {
                    self.unsaturated -= 1;
                }
            }
        }
        Ok(())
    }
}

fn read_required(input: &mut impl Read, bytes: &mut [u8], detail: &str) -> io::Result<()> {
    input.read_exact(bytes).map_err(|error| {
        if error.kind() == io::ErrorKind::UnexpectedEof {
            invalid(detail)
        } else {
            error
        }
    })
}

/// Read length-prefixed literal sequence records. A terminal summary is emitted only
/// after consuming the explicit end frame and checking that no bytes follow it.
pub(crate) fn run(
    words: impl BufRead,
    mut input: impl Read,
    mut output: impl Write,
    scope: Scope,
) -> io::Result<()> {
    let mut matcher = Matcher::new(words)?;
    let mut magic = [0; MAGIC.len()];
    read_required(&mut input, &mut magic, "truncated exact placement protocol header")?;
    if magic != MAGIC {
        return Err(invalid("invalid exact placement protocol header"));
    }
    writeln!(output, "OVXPOS1\t{}", scope.label())?;
    let mut buffer = [0; BUFFER];
    let mut record = 0u64;
    loop {
        let mut frame = [0; 4];
        read_required(&mut input, &mut frame, "truncated exact placement record frame")?;
        let length = u32::from_le_bytes(frame);
        if length == NONE {
            let mut trailing = [0; 1];
            match input.read_exact(&mut trailing) {
                Ok(()) => return Err(invalid("trailing exact placement frame data")),
                Err(error) if error.kind() == io::ErrorKind::UnexpectedEof => break,
                Err(error) => return Err(error),
            }
        }
        let mut state = 0;
        let mut offset = 0;
        let length = u64::from(length);
        while offset < length {
            let size = (length - offset).min(BUFFER as u64) as usize;
            read_required(
                &mut input,
                &mut buffer[..size],
                "truncated exact placement sequence",
            )?;
            matcher.scan_chunk(scope, record, offset, &buffer[..size], &mut state, &mut output)?;
            offset += size as u64;
        }
        record = record
            .checked_add(1)
            .ok_or_else(|| invalid("too many exact placement records"))?;
    }
    writeln!(output, "E\t{record}\t{}", matcher.patterns.len())?;
    output.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn framed(records: &[&[u8]]) -> Vec<u8> {
        let mut input = MAGIC.to_vec();
        for record in records {
            input.extend_from_slice(&(record.len() as u32).to_le_bytes());
            input.extend_from_slice(record);
        }
        input.extend_from_slice(&NONE.to_le_bytes());
        input
    }

    fn search(words: &[u8], records: &[&[u8]], scope: Scope) -> Vec<(usize, usize, usize)> {
        let mut output = Vec::new();
        run(words, framed(records).as_slice(), &mut output, scope).unwrap();
        let text = String::from_utf8(output).unwrap();
        let mut lines = text.lines();
        assert_eq!(lines.next(), Some(format!("OVXPOS1\t{}", scope.label()).as_str()));
        assert_eq!(
            lines.next_back(),
            Some(format!("E\t{}\t{}", records.len(), words.lines().count()).as_str())
        );
        lines
            .map(|line| {
                let mut fields = line.split('\t');
                assert_eq!(fields.next(), Some("H"));
                let hit = (
                    fields.next().unwrap().parse().unwrap(),
                    fields.next().unwrap().parse().unwrap(),
                    fields.next().unwrap().parse().unwrap(),
                );
                assert_eq!(fields.next(), None);
                hit
            })
            .collect()
    }

    #[test]
    fn returns_bounded_overlapping_and_nested_literals_in_both_scopes() {
        let words = b"AAA\r\nAAAA\nACGA\nTCGT\nACGT\n";
        let records = [b"AAAAAA".as_slice(), b"TCGTNACGAACGA", b"ACGTNACGT", b"AAAA"];
        let by_record = search(words, &records, Scope::Records);
        assert!(by_record.contains(&(0, 0, 0)));
        assert!(by_record.contains(&(0, 0, 1)));
        assert!(by_record.contains(&(0, 1, 1)));
        assert!(by_record.contains(&(1, 2, 9)));
        assert!(by_record.contains(&(1, 3, 0)));
        assert!(by_record.contains(&(2, 4, 0)));
        assert!(by_record.contains(&(2, 4, 5)));
        assert_eq!(by_record.iter().filter(|&&(r, p, _)| r == 0 && p == 0).count(), 2);
        assert!(by_record.contains(&(3, 0, 0)));
        assert!(!search(words, &records, Scope::Panel).iter().any(|hit| hit.0 == 3));
    }

    #[test]
    fn matches_across_buffer_boundaries_and_resets_at_records_and_n() {
        let mut long = vec![b'N'; BUFFER - 2];
        long.extend_from_slice(b"ACGT");
        let records = [b"CA".as_slice(), b"AC", b"AANC", &long];
        let hits = search(b"CAAC\nAACC\nACGT\n", &records, Scope::Records);
        assert_eq!(hits, [(3, 2, BUFFER - 2)]);
    }

    #[test]
    fn saturated_panel_requires_complete_valid_framing() {
        let mut input = framed(&[b"AAAAA", b"AAAAA"]);
        input.truncate(input.len() - 5);
        let mut output = Vec::new();
        assert!(run(b"AA\n".as_slice(), input.as_slice(), &mut output, Scope::Panel).is_err());
        assert_eq!(output, b"OVXPOS1\tpanel\nH\t0\t0\t0\nH\t0\t0\t1\n");

        let mut input = framed(&[b"AAAAA"]);
        input.push(0);
        output.clear();
        assert!(run(b"AA\n".as_slice(), input.as_slice(), &mut output, Scope::Panel).is_err());
        assert!(!output.split(|&b| b == b'\n').any(|line| line.starts_with(b"E\t")));
    }

    #[test]
    fn rejects_invalid_protocol_header_and_patterns() {
        let input = framed(&[]);
        for words in [b"".as_slice(), b"\n", b"acgt\n", b"ACGN\n", b"ACGT\nACGT\n"] {
            let mut output = Vec::new();
            assert!(run(words, input.as_slice(), &mut output, Scope::Panel).is_err());
            assert!(output.is_empty());
        }
        let mut output = Vec::new();
        assert!(run(b"A\n".as_slice(), &b"wrong-v\n"[..], &mut output, Scope::Panel).is_err());
        assert!(output.is_empty());
    }

    #[test]
    fn matches_direct_substring_oracle_on_random_literal_and_record_batches() {
        let mut state = 324205u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state >> 32) as usize
        };
        let mut words = vec![b"A".to_vec(), b"AA".to_vec(), b"AAA".to_vec()];
        for _ in 0..70 {
            let len = 2 + next() % 12;
            let word: Vec<_> = (0..len).map(|_| b"ACGT"[next() % 4]).collect();
            if !words.contains(&word) {
                words.push(word);
            }
        }
        let records: Vec<Vec<u8>> = (0..95)
            .map(|_| (0..next() % 500).map(|_| b"ACGTN"[next() % 5]).collect())
            .collect();
        let refs: Vec<_> = records.iter().map(Vec::as_slice).collect();
        let encoded: Vec<_> = words.iter().flat_map(|w| w.iter().copied().chain([b'\n'])).collect();
        for scope in [Scope::Panel, Scope::Records] {
            let mut expected = Vec::new();
            for (pattern, word) in words.iter().enumerate() {
                let mut panel_count = 0;
                for (record, sequence) in records.iter().enumerate() {
                    for (start, _) in sequence
                        .windows(word.len())
                        .enumerate()
                        .filter(|(_, window)| *window == word)
                        .take(2)
                    {
                        if scope == Scope::Panel && panel_count >= 2 {
                            break;
                        }
                        expected.push((record, pattern, start));
                        panel_count += 1;
                    }
                }
            }
            let mut observed = search(&encoded, &refs, scope);
            expected.sort_unstable();
            observed.sort_unstable();
            assert_eq!(observed, expected);
        }
    }
}
