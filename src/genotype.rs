//! Genotyping an individual against the species' reference graph: which allele of every known
//! variant site its reads carry, whether two do (heteroplasmy), and which reference segments it
//! has at all.
//!
//! A site is a simple bubble of the graph: two oriented segments A and B with two or more ways
//! between them, each through one segment or straight from A to B (that allele lacks the
//! sequence). An allele is spelled with the last bases of A and the first of B around it, and
//! its depth is the median count, in the individual's reads, of the `K`-mers of that spelling
//! that occur nowhere else in the graph. A segment's presence is the median count of its own
//! such `K`-mers. Each allele also lists the reference runs whose walks take it.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::Result;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde_json::{json, Value};

use crate::evidence::{flip, revcomp, Graph, Oriented};
use crate::kmer::CanonicalKmers;

/// Word length of the alleles' spelling and of counting.
pub const K: usize = 31;
/// The minor allele of a site is called when it has this share of the site's depth.
const MINOR_SHARE: f64 = 0.1;
/// A site is called only when its alleles together reach this share of the organelle's depth.
const MIN_SITE_SHARE: f64 = 0.2;
/// A segment is present from this share of the organelle's depth, absent below `ABSENT_SHARE`.
const PRESENT_SHARE: f64 = 0.25;
const ABSENT_SHARE: f64 = 0.02;
/// Segments with fewer own `K`-mers do not set the organelle's depth.
const DEPTH_MIN_KMERS: usize = 10;

/// One way through a site: a segment, or None straight from the left flank to the right.
type Allele = Option<Oriented>;

struct Site {
    left: Oriented,
    right: Oriented,
    alleles: Vec<Allele>,
}

fn oseq(g: &Graph, o: Oriented) -> Vec<u8> {
    if o.1 {
        g.seqs[o.0].clone()
    } else {
        revcomp(&g.seqs[o.0])
    }
}

/// The simple bubbles of a blunt graph, each once (as the lesser of it and its reverse).
fn sites(g: &Graph) -> Vec<Site> {
    let mut succ: BTreeMap<Oriented, BTreeSet<Oriented>> = BTreeMap::new();
    for &(a, b) in &g.links {
        succ.entry(a).or_default().insert(b);
        succ.entry(flip(b)).or_default().insert(flip(a));
    }
    let empty = BTreeSet::new();
    let mut out = Vec::new();
    for (&a, next) in &succ {
        let mut by_right: BTreeMap<Oriented, BTreeSet<Allele>> = BTreeMap::new();
        for &x in next {
            if x.0 != a.0 {
                by_right.entry(x).or_default().insert(None);
            }
            for &b in succ.get(&x).unwrap_or(&empty) {
                if x.0 != a.0 && b.0 != a.0 && b.0 != x.0 {
                    by_right.entry(b).or_default().insert(Some(x));
                }
            }
        }
        for (b, alleles) in by_right {
            if alleles.len() < 2 || alleles.iter().all(Option::is_none) || (a, b) > (flip(b), flip(a)) {
                continue;
            }
            out.push(Site {
                left: a,
                right: b,
                alleles: alleles.into_iter().collect(),
            });
        }
    }
    out
}

/// An allele spelled with up to `K - 1` bases of its flanks on either side.
fn spell(g: &Graph, s: &Site, allele: Allele) -> Vec<u8> {
    let (l, r) = (oseq(g, s.left), oseq(g, s.right));
    let mut seq = l[l.len().saturating_sub(K - 1)..].to_vec();
    if let Some(x) = allele {
        seq.extend(oseq(g, x));
    }
    seq.extend_from_slice(&r[..r.len().min(K - 1)]);
    seq
}

fn kmers(seq: &[u8]) -> Vec<u64> {
    CanonicalKmers::new(seq, K).map(|(_, x)| x).collect()
}

fn median(mut v: Vec<u32>) -> Option<u32> {
    v.sort_unstable();
    (!v.is_empty()).then(|| v[v.len() / 2])
}

/// Genotype the reads of `reads` against the reference graph `g`, whose walks of each run are
/// `walks` (name, steps). Writes `genotypes.tsv` and `presence.tsv` into `out`.
pub fn run(g: &Graph, walks: &[(String, Vec<Oriented>)], reads: &[PathBuf], out: &Path) -> Result<Value> {
    let sites = sites(g);
    // every place a K-mer occurs in the graph: segment interiors and junctions, each once
    let mut places: FxHashMap<u64, u32> = FxHashMap::default();
    for s in &g.seqs {
        for x in kmers(s) {
            *places.entry(x).or_insert(0) += 1;
        }
    }
    for &(a, b) in &g.links {
        let (sa, sb) = (oseq(g, a), oseq(g, b));
        let junction = [&sa[sa.len().saturating_sub(K - 1)..], &sb[..sb.len().min(K - 1)]].concat();
        for x in kmers(&junction) {
            *places.entry(x).or_insert(0) += 1;
        }
    }
    let own = |seq: &[u8]| -> Vec<u64> {
        let mut v: Vec<u64> = kmers(seq).into_iter().filter(|x| places.get(x) == Some(&1)).collect();
        v.sort_unstable();
        v.dedup();
        v
    };
    let node_kmers: Vec<Vec<u64>> = g.seqs.par_iter().map(|s| own(s)).collect();
    let allele_kmers: Vec<Vec<Vec<u64>>> = sites
        .par_iter()
        .map(|s| s.alleles.iter().map(|&a| own(&spell(g, s, a))).collect())
        .collect();
    // count those K-mers in the reads
    let mut counts: FxHashMap<u64, u32> = FxHashMap::default();
    for x in node_kmers.iter().flatten().chain(allele_kmers.iter().flatten().flatten()) {
        counts.insert(*x, 0);
    }
    let (batches, reader) = crate::recruit::read_batches(reads);
    let mut reads_seen = 0u64;
    for batch in batches {
        reads_seen += batch.len() as u64;
        let hits: Vec<Vec<u64>> = batch
            .par_iter()
            .map(|r| kmers(&r.seq).into_iter().filter(|x| counts.contains_key(x)).collect())
            .collect();
        for x in hits.into_iter().flatten() {
            *counts.get_mut(&x).unwrap() += 1;
        }
    }
    reader
        .join()
        .map_err(|_| anyhow::anyhow!("genotyping reader panicked"))??;
    let depth_of = |ks: &[u64]| median(ks.iter().map(|x| counts[x]).collect());
    let node_depth: Vec<Option<u32>> = node_kmers.iter().map(|ks| depth_of(ks)).collect();
    let organelle = median(
        node_depth
            .iter()
            .zip(&node_kmers)
            .filter(|(d, ks)| ks.len() >= DEPTH_MIN_KMERS && d.is_some_and(|d| d > 0))
            .map(|(d, _)| d.unwrap())
            .collect(),
    )
    .unwrap_or(0) as f64;
    // which reference runs take each allele
    let run_of = |name: &str| name.split('#').next().unwrap_or(name).to_string();
    let mut taken: BTreeMap<(Oriented, Allele, Oriented), BTreeSet<String>> = BTreeMap::new();
    for (name, steps) in walks {
        let run = run_of(name);
        for w in steps.windows(2) {
            for (a, x, b) in [(w[0], None, w[1]), (flip(w[1]), None, flip(w[0]))] {
                taken.entry((a, x, b)).or_default().insert(run.clone());
            }
        }
        for w in steps.windows(3) {
            for (a, x, b) in [(w[0], Some(w[1]), w[2]), (flip(w[2]), Some(flip(w[1])), flip(w[0]))] {
                taken.entry((a, x, b)).or_default().insert(run.clone());
            }
        }
    }
    let label = |a: Allele| a.map_or("-".to_string(), |o| g.label(o));
    let mut tsv = String::from(
        "site\tleft\tright\tallele\tallele_bp\treference_runs\tkmers\tdepth\tshare\tcall\n",
    );
    let (mut called, mut no_call, mut ungenotyped, mut heteroplasmic, mut absent_from_runs) =
        (0usize, 0usize, 0usize, 0usize, 0usize);
    let mut minor_list = Vec::new();
    for (i, (s, ks)) in sites.iter().zip(&allele_kmers).enumerate() {
        let depths: Vec<Option<u32>> = ks.iter().map(|k| depth_of(k)).collect();
        if depths.iter().any(Option::is_none) {
            ungenotyped += 1;
        }
        let total: u32 = depths.iter().flatten().sum();
        let enough = organelle > 0.0 && total as f64 >= MIN_SITE_SHARE * organelle && depths.iter().all(Option::is_some);
        let major = depths.iter().enumerate().max_by_key(|(_, d)| d.unwrap_or(0)).map(|(j, _)| j);
        let mut minors = 0;
        for (j, &allele) in s.alleles.iter().enumerate() {
            let share = depths[j].map(|d| d as f64 / total.max(1) as f64);
            let call = if !enough {
                "no_call"
            } else if Some(j) == major {
                "major"
            } else if share.is_some_and(|f| f >= MINOR_SHARE) {
                minors += 1;
                "minor"
            } else {
                "absent"
            };
            let runs = taken.get(&(s.left, allele, s.right)).cloned().unwrap_or_default();
            if call == "major" && runs.is_empty() {
                absent_from_runs += 1;
            }
            writeln!(
                tsv,
                "{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{call}",
                g.label(s.left),
                g.label(s.right),
                label(allele),
                allele.map_or(0, |o| g.seqs[o.0].len()),
                if runs.is_empty() { "-".to_string() } else { runs.into_iter().collect::<Vec<_>>().join(",") },
                ks[j].len(),
                depths[j].map_or("-".to_string(), |d| d.to_string()),
                share.map_or("-".to_string(), |f| format!("{f:.3}")),
            )?;
        }
        if enough {
            called += 1;
            if minors > 0 {
                heteroplasmic += 1;
                minor_list.push(i);
            }
        } else {
            no_call += 1;
        }
    }
    fs::write(out.join("genotypes.tsv"), tsv)?;
    let mut presence = String::from("segment\tlength\tkmers\tdepth\trelative\tcall\n");
    let (mut present, mut absent, mut low, mut unknown) = (0usize, 0usize, 0usize, 0usize);
    for (i, (d, ks)) in node_depth.iter().zip(&node_kmers).enumerate() {
        let rel = d.map(|d| d as f64 / organelle.max(1.0));
        let call = match rel {
            None => {
                unknown += 1;
                "no_own_kmers"
            }
            Some(r) if r >= PRESENT_SHARE => {
                present += 1;
                "present"
            }
            Some(r) if r < ABSENT_SHARE => {
                absent += 1;
                "absent"
            }
            Some(_) => {
                low += 1;
                "low"
            }
        };
        writeln!(
            presence,
            "{}\t{}\t{}\t{}\t{}\t{call}",
            g.names[i],
            g.seqs[i].len(),
            ks.len(),
            d.map_or("-".to_string(), |d| d.to_string()),
            rel.map_or("-".to_string(), |r| format!("{r:.3}")),
        )?;
    }
    fs::write(out.join("presence.tsv"), presence)?;
    Ok(json!({"k": K, "reads": reads_seen, "organelle_kmer_depth": organelle,
        "sites": sites.len(), "called": called, "no_call": no_call,
        "sites_with_an_allele_without_own_kmers": ungenotyped,
        "heteroplasmic_sites": heteroplasmic, "heteroplasmic": minor_list,
        "major_alleles_in_no_reference_run": absent_from_runs,
        "segments": {"present": present, "low": low, "absent": absent, "no_own_kmers": unknown},
        "genotypes": out.join("genotypes.tsv"), "presence": out.join("presence.tsv"),
        "basis": "median read count of each allele's K-mers that occur once in the reference graph"}))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random(n: usize, seed: u64) -> String {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = x;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                b"ACGT"[((z ^ (z >> 31)) >> 62) as usize] as char
            })
            .collect()
    }

    #[test]
    fn a_snp_an_indel_and_a_missing_segment_are_genotyped() {
        // L [A|C] M [ins|-] N, plus a segment X only one run has after N
        let (l, m, n, x) = (random(300, 1), random(300, 2), random(300, 3), random(200, 4));
        let ins = random(40, 5);
        let gfa = format!(
            "S\tL\t{l}\nS\ta\tA\nS\tc\tC\nS\tM\t{m}\nS\ti\t{ins}\nS\tN\t{n}\nS\tX\t{x}\n\
             L\tL\t+\ta\t+\t0M\nL\tL\t+\tc\t+\t0M\nL\ta\t+\tM\t+\t0M\nL\tc\t+\tM\t+\t0M\n\
             L\tM\t+\ti\t+\t0M\nL\ti\t+\tN\t+\t0M\nL\tM\t+\tN\t+\t0M\nL\tN\t+\tX\t+\t0M\n"
        );
        let dir = std::env::temp_dir().join(format!("ovasm-genotype-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let gp = dir.join("g.gfa");
        fs::write(&gp, &gfa).unwrap();
        let g = crate::evidence::parse_gfa(&gp).unwrap();
        let id = |n: &str| g.names.iter().position(|x| x == n).unwrap();
        let walks = vec![
            ("one#1#mt".to_string(), vec![(id("L"), true), (id("a"), true), (id("M"), true), (id("i"), true), (id("N"), true), (id("X"), true)]),
            ("two#1#mt".to_string(), vec![(id("L"), true), (id("c"), true), (id("M"), true), (id("N"), true)]),
        ];
        // the individual: C at the SNP, the insertion, no X; 80% C and 20% A (heteroplasmy)
        let ind = format!("{l}C{m}{ins}{n}");
        let minor = format!("{l}A{m}{ins}{n}");
        let mut fq = String::new();
        for (genome, copies) in [(&ind, 8), (&minor, 2)] {
            for _ in 0..copies {
                for a in (0..genome.len() - 100).step_by(5) {
                    let r = &genome[a..a + 100];
                    writeln!(fq, "@r\n{r}\n+\n{}", "I".repeat(100)).unwrap();
                }
            }
        }
        let reads = dir.join("r.fq");
        fs::write(&reads, fq).unwrap();
        let s = run(&g, &walks, &[reads], &dir).unwrap();
        assert_eq!(s["sites"], 2, "{s}");
        assert_eq!(s["called"], 2, "{s}");
        assert_eq!(s["heteroplasmic_sites"], 1, "{s}");
        let tsv = fs::read_to_string(dir.join("genotypes.tsv")).unwrap();
        let call = |allele: &str| {
            tsv.lines()
                .find(|l| l.split('\t').nth(3) == Some(allele))
                .map(|l| l.rsplit('\t').next().unwrap().to_string())
                .unwrap()
        };
        assert_eq!(call("c+"), "major");
        assert_eq!(call("a+"), "minor");
        assert_eq!(call("i+"), "major");
        assert_eq!(call("-"), "absent");
        let presence = fs::read_to_string(dir.join("presence.tsv")).unwrap();
        assert!(presence.lines().any(|l| l.starts_with("X\t") && l.ends_with("\tabsent")), "{presence}");
        assert!(presence.lines().any(|l| l.starts_with("M\t") && l.ends_with("\tpresent")), "{presence}");
        assert!(tsv.lines().any(|l| l.split('\t').nth(3) == Some("c+") && l.split('\t').nth(5) == Some("two")));
        fs::remove_dir_all(dir).unwrap();
    }
}
