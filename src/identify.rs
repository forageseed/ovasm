//! Closest reference genomes of an assembly (`ovasm identify`).
//!
//! An assembled organelle genome is compared with every reference genome of one or more sets
//! (e.g. the mitochondrial and plastid SeedDB) by k-mer containment: the fraction of the
//! assembly's distinct canonical k-mers present in the reference. Containment C converts to an
//! identity estimate C^(1/k) (Mash Screen; Ondov et al. 2019), which holds where the assembly
//! aligns to the reference and only ranks references beyond that. The report answers two
//! questions the assembly graph cannot: which organelle set the result resembles (a
//! "mitochondrion" closer to the plastid references is a target error) and which species it
//! resembles (a near-identical match to a reference of another species points to contamination
//! or a mislabelled sample).

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use rustc_hash::FxHashSet;
use serde::Serialize;

use crate::kmer::CanonicalKmers;

#[derive(Debug, Clone, Serialize)]
pub struct ReferenceHit {
    pub set: String,
    pub id: String,
    pub length: u64,
    /// Distinct query k-mers present in this reference.
    pub shared: usize,
    pub containment: f64,
    pub estimated_identity: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct IdentifyReport {
    pub tool: &'static str,
    pub version: &'static str,
    pub query: String,
    pub k: usize,
    pub query_length: u64,
    pub query_kmers: usize,
    /// Every reference, highest containment first.
    pub references: Vec<ReferenceHit>,
}

/// C^(1/k): the per-base identity at which a k-mer survives with probability C.
pub fn identity_from_containment(containment: f64, k: usize) -> f64 {
    if containment <= 0.0 {
        0.0
    } else {
        containment.powf(1.0 / k as f64)
    }
}

fn kmers_of(path: &Path, k: usize) -> Result<(FxHashSet<u64>, u64)> {
    let mut reader = needletail::parse_fastx_file(path)
        .with_context(|| format!("cannot read {}", path.display()))?;
    let (mut set, mut length) = (FxHashSet::default(), 0u64);
    while let Some(record) = reader.next() {
        let record = record.with_context(|| format!("malformed record in {}", path.display()))?;
        let seq = record.seq();
        length += seq.len() as u64;
        set.extend(CanonicalKmers::new(&seq, k).map(|(_, km)| km));
    }
    Ok((set, length))
}

pub fn run(
    query: &Path,
    references: &[(String, PathBuf)],
    k: usize,
    out: &Path,
) -> Result<IdentifyReport> {
    let (wanted, query_length) = kmers_of(query, k)?;
    let mut hits = Vec::new();
    for (set, path) in references {
        let mut reader = needletail::parse_fastx_file(path)
            .with_context(|| format!("cannot read {}", path.display()))?;
        while let Some(record) = reader.next() {
            let record =
                record.with_context(|| format!("malformed record in {}", path.display()))?;
            let seq = record.seq();
            let shared: FxHashSet<u64> = CanonicalKmers::new(&seq, k)
                .map(|(_, km)| km)
                .filter(|km| wanted.contains(km))
                .collect();
            let containment = if wanted.is_empty() {
                0.0
            } else {
                shared.len() as f64 / wanted.len() as f64
            };
            hits.push(ReferenceHit {
                set: set.clone(),
                id: String::from_utf8_lossy(record.id())
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string(),
                length: seq.len() as u64,
                shared: shared.len(),
                containment,
                estimated_identity: identity_from_containment(containment, k),
            });
        }
    }
    hits.sort_by(|a, b| {
        b.containment
            .total_cmp(&a.containment)
            .then_with(|| a.id.cmp(&b.id))
    });
    let report = IdentifyReport {
        tool: "ovasm identify",
        version: env!("CARGO_PKG_VERSION"),
        query: query.display().to_string(),
        k,
        query_length,
        query_kmers: wanted.len(),
        references: hits,
    };
    std::fs::write(out, serde_json::to_string_pretty(&report)?)?;
    Ok(report)
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

    /// `g` with a substitution every `every` bases.
    fn diverged(g: &[u8], every: usize) -> Vec<u8> {
        let mut d = g.to_vec();
        for i in (every / 2..d.len()).step_by(every) {
            d[i] = if d[i] == b'A' { b'C' } else { b'A' };
        }
        d
    }

    #[test]
    fn ranks_the_closest_reference_and_estimates_its_identity() {
        let dir = std::env::temp_dir().join(format!("ovasm-identify-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mito = genome(20_000, 1);
        let plastid = genome(20_000, 2);
        let fa = |name: &str, recs: &[(&str, &[u8])]| {
            let p = dir.join(name);
            let text: String = recs
                .iter()
                .map(|(id, s)| format!(">{id} desc\n{}\n", std::str::from_utf8(s).unwrap()))
                .collect();
            std::fs::write(&p, text).unwrap();
            p
        };
        // query: the plastid at ~99% identity (a substitution every 100 bp)
        let query = fa("q.fa", &[("q", &diverged(&plastid, 100))]);
        let mt = fa("mt.fa", &[("mtA", &mito), ("mtB", &genome(20_000, 3))]);
        let pt = fa(
            "pt.fa",
            &[("ptNear", &plastid), ("ptFar", &diverged(&plastid, 40))],
        );
        let refs = vec![
            ("mitochondrion".to_string(), mt),
            ("plastid".to_string(), pt),
        ];
        let r = run(&query, &refs, 21, &dir.join("id.json")).unwrap();
        assert_eq!(r.references[0].id, "ptNear");
        assert_eq!(r.references[0].set, "plastid");
        assert_eq!(r.references[1].id, "ptFar");
        let id = r.references[0].estimated_identity;
        assert!((id - 0.99).abs() < 0.005, "estimated identity {id}");
        assert!(r.references[2..]
            .iter()
            .all(|h| h.set == "mitochondrion" && h.shared == 0));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn identity_inverts_kmer_survival() {
        assert_eq!(identity_from_containment(0.0, 21), 0.0);
        assert!((identity_from_containment(0.95f64.powi(21), 21) - 0.95).abs() < 1e-12);
    }
}
