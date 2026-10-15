//! Whole original-record evidence of a supplied nuclear alignment. Matching shared DNA
//! alone is insufficient: the same continuous alignment must carry target-absent context
//! on BOTH sides of the junction. Missing evidence stays unknown, never nuclear by default.
use anyhow::{bail, Context, Result};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};
use std::{fs, path::Path};

use crate::{
    evidence::{self, EvidenceParams, Graph},
    kmer::CanonicalKmers,
    recruit::Rec,
};

#[derive(Deserialize)]
struct Document {
    schema: String,
    records: Vec<Record>,
}

/// Reference DNA is the exact [reference_start, reference_end) aligned substring,
/// forward in reference coordinates. CIGAR is extended (= X I D), target-forward.
/// query_sequence binds the evidence to a complete original input, not a local hit.
#[derive(Deserialize)]
struct Record {
    read_id: String,
    query_sequence: String,
    query_start: usize,
    query_end: usize,
    strand: String,
    reference_name: String,
    reference_start: usize,
    reference_sequence: String,
    cigar: String,
    #[serde(default)]
    junctions: Vec<(i64, i64)>,
}

#[derive(Serialize, Clone)]
pub struct ContextSummary {
    pub schema: &'static str,
    pub supplied_complete_records: usize,
    pub qualified_continuous_alignments: usize,
    pub minimum_unique_flank: usize,
    pub minimum_identity_and_query_coverage: f64,
    pub physical_origin_certified: bool,
    pub records: Vec<ContextRecord>,
}

#[derive(Serialize, Clone)]
pub struct ContextRecord {
    pub read_id: String,
    pub reference_name: String,
    pub qualified_continuous_alignment: bool,
    pub target_absent_query_intervals: Vec<(usize, usize)>,
    pub supplied_junctions: Vec<(i64, i64)>,
    pub supplied_junction_nuclear_flanked: Vec<bool>,
}

fn flanked(intervals: &[(usize, usize)], q0: i64, q1: i64) -> bool {
    let (lo, hi) = (q0.min(q1), q0.max(q1));
    intervals.iter().any(|&(_, b)| b as i64 <= lo) && intervals.iter().any(|&(a, _)| a as i64 >= hi)
}

struct Checked {
    sequence: Vec<u8>,
    intervals: Vec<(usize, usize)>,
}
pub struct ReadContext {
    reads: FxHashMap<String, Checked>,
    pub summary: ContextSummary,
}

fn target_words(g: &Graph, k: usize) -> FxHashSet<u64> {
    let mut words: FxHashSet<_> = g
        .seqs
        .iter()
        .flat_map(|s| CanonicalKmers::new(s, k).map(|(_, w)| w))
        .collect();
    // A word straddling an existing edge is target DNA too.
    for (i, &(a, b)) in g.links.iter().enumerate() {
        let x = g.oriented_seq(a);
        let y = g.oriented_seq(b);
        let ov = g.overlaps[i];
        let mut seq = x[x.len().saturating_sub(k - 1)..].to_vec();
        seq.extend_from_slice(&y[ov..y.len().min(ov + k - 1)]);
        words.extend(CanonicalKmers::new(&seq, k).map(|(_, w)| w));
    }
    words
}

fn check(r: &Record, words: &FxHashSet<u64>, p: &EvidenceParams) -> Result<(Checked, bool)> {
    if r.read_id.is_empty()
        || r.read_id.chars().any(char::is_whitespace)
        || r.reference_name.is_empty()
        || !matches!(r.strand.as_str(), "+" | "-")
    {
        bail!("invalid original-record or nuclear-reference identity");
    }
    let raw = r.query_sequence.as_bytes().to_ascii_uppercase();
    if !(r.query_start < r.query_end && r.query_end <= raw.len()) {
        bail!("nuclear context query coordinates exceed complete original read");
    }
    let query = if r.strand == "+" {
        raw.clone()
    } else {
        evidence::revcomp(&raw)
    };
    let reference = r.reference_sequence.as_bytes().to_ascii_uppercase();
    r.reference_start
        .checked_add(reference.len())
        .context("nuclear reference coordinate overflow")?;
    let mut shared = vec![false; reference.len()];
    for (pos, word) in CanonicalKmers::new(&reference, p.k) {
        if words.contains(&word) {
            shared[pos..pos + p.k].fill(true);
        }
    }
    let start = if r.strand == "+" {
        r.query_start
    } else {
        raw.len() - r.query_end
    };
    let end = if r.strand == "+" {
        r.query_end
    } else {
        raw.len() - r.query_start
    };
    let (mut q, mut t, mut n, mut matches, mut aligned, mut columns) =
        (start, 0usize, 0usize, 0usize, 0usize, 0usize);
    let mut unique = vec![false; raw.len()];
    let mut operations = 0;
    for c in r.cigar.bytes() {
        if c.is_ascii_digit() {
            n = n
                .checked_mul(10)
                .and_then(|x| x.checked_add((c - b'0') as usize))
                .context("CIGAR overflow")?;
            continue;
        }
        if n == 0 {
            bail!("empty or zero-length nuclear CIGAR operation");
        }
        let qn = q.checked_add(n).context("CIGAR overflow")?;
        let tn = t.checked_add(n).context("CIGAR overflow")?;
        match c {
            b'=' | b'X' => {
                if qn > end || tn > reference.len() {
                    bail!("nuclear CIGAR exceeds original query/reference");
                }
                for i in 0..n {
                    let equal = query[q + i] == reference[t + i];
                    if equal != (c == b'=') {
                        bail!("extended nuclear CIGAR disagrees with actual query/reference DNA");
                    }
                    if equal {
                        matches += 1;
                    }
                    // Both matching and mismatching aligned bases have checked coordinates;
                    // overall error tolerance is the existing evidence spacing tolerance.
                    unique[q + i] = b"ACGT".contains(&query[q + i])
                        && b"ACGT".contains(&reference[t + i])
                        && !shared[t + i];
                }
                q = qn;
                t = tn;
                aligned += n;
            }
            b'I' => {
                if qn > end {
                    bail!("nuclear insertion exceeds original query");
                }
                q = qn;
            }
            b'D' => {
                if tn > reference.len() {
                    bail!("nuclear deletion exceeds reference");
                }
                t = tn;
            }
            _ => {
                bail!("nuclear context requires explicit = X I D CIGAR, not clipping, M or splice")
            }
        }
        columns += n;
        n = 0;
        operations += 1;
    }
    if n != 0 || operations == 0 || q != end || t != reference.len() {
        bail!("incomplete nuclear CIGAR");
    }
    let qualified = aligned as f64 / raw.len() as f64 >= 1.0 - p.spacing_tol
        && matches as f64 / columns.max(1) as f64 >= 1.0 - p.spacing_tol;
    let mut intervals = Vec::new();
    let mut a = 0;
    while a < unique.len() {
        if !unique[a] {
            a += 1;
            continue;
        }
        let mut b = a + 1;
        while b < unique.len() && unique[b] {
            b += 1;
        }
        if qualified && b - a >= p.min_anchor {
            intervals.push(if r.strand == "+" {
                (a, b)
            } else {
                (raw.len() - b, raw.len() - a)
            });
        }
        a = b;
    }
    intervals.sort_unstable();
    Ok((
        Checked {
            sequence: raw,
            intervals,
        },
        qualified,
    ))
}

impl ReadContext {
    pub fn load(path: &Path, g: &Graph, p: &EvidenceParams) -> Result<Self> {
        let doc: Document = serde_json::from_slice(&fs::read(path)?)?;
        if doc.schema != "ovasm.whole-read-context.v1" {
            bail!("unsupported complete-read context schema");
        }
        let words = target_words(g, p.k);
        let mut reads = FxHashMap::default();
        let mut count = 0;
        let mut rows = Vec::new();
        for record in doc.records {
            let (checked, qualified) = check(&record, &words, p)
                .with_context(|| format!("nuclear context for {}", record.read_id))?;
            if qualified {
                count += 1;
            }
            for &(a, b) in &record.junctions {
                if a.min(b) < 0 || a.max(b) > checked.sequence.len() as i64 {
                    bail!("source junction exceeds original query");
                }
            }
            rows.push(ContextRecord {
                read_id: record.read_id.clone(),
                reference_name: record.reference_name.clone(),
                qualified_continuous_alignment: qualified,
                target_absent_query_intervals: checked.intervals.clone(),
                supplied_junction_nuclear_flanked: record
                    .junctions
                    .iter()
                    .map(|&(a, b)| flanked(&checked.intervals, a, b))
                    .collect(),
                supplied_junctions: record.junctions,
            });
            if reads.insert(record.read_id, checked).is_some() {
                bail!("duplicate original read in nuclear context");
            }
        }
        Ok(Self {
            summary: ContextSummary {
                schema: "ovasm.whole-read-context.v1",
                supplied_complete_records: reads.len(),
                qualified_continuous_alignments: count,
                minimum_unique_flank: p.min_anchor,
                minimum_identity_and_query_coverage: 1.0 - p.spacing_tol,
                physical_origin_certified: false,
                records: rows,
            },
            reads,
        })
    }
    pub fn nuclear_flanked(&self, r: &Rec, q0: i64, q1: i64) -> Result<bool> {
        let id = std::str::from_utf8(
            r.id.split(|b| b.is_ascii_whitespace())
                .next()
                .unwrap_or(&[]),
        )?;
        let Some(context) = self.reads.get(id) else {
            return Ok(false);
        };
        if context.sequence != r.seq.to_ascii_uppercase() {
            bail!("nuclear evidence does not bind the complete original record {id}");
        }
        Ok(flanked(&context.intervals, q0, q1))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn record(strand: &str) -> Record {
        let s = "ACGT".repeat(300);
        Record {
            read_id: "r".into(),
            query_sequence: s.clone(),
            query_start: 0,
            query_end: s.len(),
            strand: strand.into(),
            reference_name: "nucleus".into(),
            reference_start: 100,
            reference_sequence: s,
            cigar: "1200=".into(),
            junctions: vec![],
        }
    }
    #[test]
    fn complete_alignment_is_bound_to_actual_dna() {
        let mut r = record("+");
        let p = crate::junctions::hifi_evidence();
        assert!(check(&r, &FxHashSet::default(), &p).unwrap().1);
        r.reference_sequence.replace_range(600..601, "T");
        assert!(check(&r, &FxHashSet::default(), &p).is_err());
    }
    #[test]
    fn shared_sequence_is_not_nuclear_specific_context() {
        let r = record("+");
        let p = crate::junctions::hifi_evidence();
        let words = CanonicalKmers::new(r.reference_sequence.as_bytes(), p.k)
            .map(|(_, w)| w)
            .collect();
        assert!(check(&r, &words, &p).unwrap().0.intervals.is_empty());
    }
    #[test]
    fn ambiguous_bases_are_not_nuclear_specific_anchors() {
        let mut r = record("+");
        r.query_sequence = "N".repeat(1200);
        r.reference_sequence = r.query_sequence.clone();
        let (checked, _) = check(
            &r,
            &FxHashSet::default(),
            &crate::junctions::hifi_evidence(),
        )
        .unwrap();
        assert!(checked.intervals.is_empty());
    }
    #[test]
    fn local_match_cannot_certify_whole_record() {
        let mut r = record("-");
        r.query_end = 600;
        r.reference_sequence.truncate(600);
        r.cigar = "600=".into();
        let (checked, qualified) = check(
            &r,
            &FxHashSet::default(),
            &crate::junctions::hifi_evidence(),
        )
        .unwrap();
        assert!(!qualified);
        assert!(checked.intervals.is_empty());
    }
}
