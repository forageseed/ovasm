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
    // each site: the major allele when called, and the runs taking each allele (`mosaic`)
    let mut observed: Vec<Option<usize>> = Vec::new();
    let mut allele_runs: Vec<Vec<BTreeSet<String>>> = Vec::new();
    let mut site_calls: Vec<Vec<(&'static str, Option<f64>)>> = Vec::new();
    for (i, (s, ks)) in sites.iter().zip(&allele_kmers).enumerate() {
        let depths: Vec<Option<u32>> = ks.iter().map(|k| depth_of(k)).collect();
        if depths.iter().any(Option::is_none) {
            ungenotyped += 1;
        }
        let total: u32 = depths.iter().flatten().sum();
        let enough = organelle > 0.0 && total as f64 >= MIN_SITE_SHARE * organelle && depths.iter().all(Option::is_some);
        let major = depths.iter().enumerate().max_by_key(|(_, d)| d.unwrap_or(0)).map(|(j, _)| j);
        let mut minors = 0;
        observed.push(major.filter(|_| enough));
        allele_runs.push(Vec::new());
        site_calls.push(Vec::new());
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
            allele_runs[i].push(runs.clone());
            site_calls[i].push((call, share));
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
    let mosaic_summary = write_mosaic(g, &sites, walks, &observed, &allele_runs, &site_calls, out)?;
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
        "major_alleles_in_no_reference_run": absent_from_runs, "mosaic": mosaic_summary,
        "segments": {"present": present, "low": low, "absent": absent, "no_own_kmers": unknown},
        "genotypes": out.join("genotypes.tsv"), "presence": out.join("presence.tsv"),
        "basis": "median read count of each allele's K-mers that occur once in the reference graph"}))
}

/// A run's allele differs from the reads' call with this probability (`mosaic`).
const MISMATCH: f64 = 0.01;
/// The individual's closest run changes between neighbouring sites with this probability: an
/// organelle genome is inherited whole, so a mosaic of runs is rare.
const SWITCH: f64 = 0.001;

/// The Li-Stephens mosaic of an individual over its reference runs: at each site of a chain the
/// runs are states, their alleles (`hap[t][r]`, empty where a run's walks do not pass the site)
/// are matched with the call (`obs[t]`) at probability 1 - `MISMATCH`, and the state changes by
/// `SWITCH` between sites. Returns, per site, the runs' posterior without the site's own call
/// (to predict it), and the Viterbi run.
fn mosaic(hap: &[Vec<Vec<usize>>], obs: &[Option<usize>], runs: usize) -> (Vec<Vec<f64>>, Vec<usize>) {
    let n = hap.len();
    if n == 0 || runs == 0 {
        return (Vec::new(), Vec::new());
    }
    let r = runs as f64;
    let e: Vec<Vec<f64>> = (0..n)
        .map(|t| {
            (0..runs)
                .map(|s| match obs[t] {
                    Some(a) if !hap[t][s].is_empty() => {
                        if hap[t][s].contains(&a) {
                            1.0 - MISMATCH
                        } else {
                            MISMATCH
                        }
                    }
                    _ => 1.0,
                })
                .collect()
        })
        .collect();
    let norm = |v: &mut Vec<f64>| {
        let s: f64 = v.iter().sum();
        v.iter_mut().for_each(|x| *x /= s);
    };
    let mut alpha = vec![vec![0.0; runs]; n];
    for s in 0..runs {
        alpha[0][s] = e[0][s] / r;
    }
    norm(&mut alpha[0]);
    for t in 1..n {
        for s in 0..runs {
            alpha[t][s] = e[t][s] * ((1.0 - SWITCH) * alpha[t - 1][s] + SWITCH / r);
        }
        norm(&mut alpha[t]);
    }
    let mut beta = vec![vec![1.0; runs]; n];
    for t in (0..n - 1).rev() {
        let w: Vec<f64> = (0..runs).map(|s| e[t + 1][s] * beta[t + 1][s]).collect();
        let sw: f64 = w.iter().sum();
        for s in 0..runs {
            beta[t][s] = (1.0 - SWITCH) * w[s] + SWITCH * sw / r;
        }
        norm(&mut beta[t]);
    }
    let without: Vec<Vec<f64>> = (0..n)
        .map(|t| {
            let mut p: Vec<f64> = (0..runs).map(|s| alpha[t][s] / e[t][s] * beta[t][s]).collect();
            norm(&mut p);
            p
        })
        .collect();
    // Viterbi, in logs
    let stay = ((1.0 - SWITCH) + SWITCH / r).ln();
    let leave = (SWITCH / r).ln();
    let mut v: Vec<f64> = (0..runs).map(|s| (1.0 / r).ln() + e[0][s].ln()).collect();
    let mut back = vec![vec![0usize; runs]; n];
    for t in 1..n {
        let mut next = vec![0.0; runs];
        for s in 0..runs {
            let (from, best) = (0..runs)
                .map(|q| (q, v[q] + if q == s { stay } else { leave }))
                .max_by(|a, b| a.1.total_cmp(&b.1))
                .unwrap();
            back[t][s] = from;
            next[s] = best + e[t][s].ln();
        }
        v = next;
    }
    let mut path = vec![0usize; n];
    path[n - 1] = (0..runs).max_by(|&a, &b| v[a].total_cmp(&v[b])).unwrap();
    for t in (1..n).rev() {
        path[t - 1] = back[t][path[t]];
    }
    (without, path)
}

/// The allele the runs predict at a site, from their weights, and its probability.
fn predict(hap: &[Vec<usize>], weights: &[f64], alleles: usize) -> Option<(usize, f64)> {
    let mut p = vec![0.0; alleles];
    let mut informative = 0.0;
    for (s, w) in weights.iter().enumerate() {
        if hap[s].is_empty() {
            continue;
        }
        informative += w;
        for &a in &hap[s] {
            p[a] += w / hap[s].len() as f64;
        }
    }
    if informative <= 0.0 {
        return None;
    }
    p.iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(a, &x)| (a, x / informative))
}

/// The individual as a mosaic of the reference runs (`mosaic`), its sites in the order of the
/// run walk that passes most of them (sites off it are predicted from the runs' overall
/// weights): which run each stretch is most like, every site predicted without its own call
/// (agreement with the called ones measures how far the runs predict this individual), and the
/// uncalled sites imputed. Writes `mosaic.tsv`.
fn write_mosaic(
    g: &Graph,
    sites: &[Site],
    walks: &[(String, Vec<Oriented>)],
    observed: &[Option<usize>],
    allele_runs: &[Vec<BTreeSet<String>>],
    calls: &[Vec<(&'static str, Option<f64>)>],
    out: &Path,
) -> Result<Value> {
    let runs: Vec<String> = allele_runs
        .iter()
        .flatten()
        .flatten()
        .cloned()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let run_index: BTreeMap<&str, usize> = runs.iter().enumerate().map(|(i, r)| (r.as_str(), i)).collect();
    // hap[site][run]: the alleles the run takes there
    let hap: Vec<Vec<Vec<usize>>> = allele_runs
        .iter()
        .map(|alleles| {
            let mut h = vec![Vec::new(); runs.len()];
            for (a, rs) in alleles.iter().enumerate() {
                for r in rs {
                    h[run_index[r.as_str()]].push(a);
                }
            }
            h
        })
        .collect();
    // the backbone: the walk passing most sites, and each site's place on it
    // (site, allele, read reversed?) of each way through a site
    let mut keys: BTreeMap<(Oriented, Allele, Oriented), (usize, usize, bool)> = BTreeMap::new();
    for (i, s) in sites.iter().enumerate() {
        for (j, &a) in s.alleles.iter().enumerate() {
            keys.insert((s.left, a, s.right), (i, j, false));
            keys.insert((flip(s.right), a.map(flip), flip(s.left)), (i, j, true));
        }
    }
    // each site a walk passes: (its step before the site, the allele it takes, reversed?)
    let places = |steps: &[Oriented]| -> BTreeMap<usize, (usize, usize, bool)> {
        let mut at = BTreeMap::new();
        for (p, w) in steps.windows(2).enumerate() {
            if let Some(&(i, j, rev)) = keys.get(&(w[0], None, w[1])) {
                at.entry(i).or_insert((p, j, rev));
            }
        }
        for (p, w) in steps.windows(3).enumerate() {
            if let Some(&(i, j, rev)) = keys.get(&(w[0], Some(w[1]), w[2])) {
                at.entry(i).or_insert((p, j, rev));
            }
        }
        at
    };
    // the first of the walks passing most sites (the VCF's coordinates do not change with the
    // order equal walks come in)
    let (_, backbone, at) = walks
        .iter()
        .enumerate()
        .map(|(k, (name, steps))| (k, name.clone(), places(steps)))
        .max_by_key(|(k, _, at)| (at.len(), std::cmp::Reverse(*k)))
        .unwrap_or_default();
    let mut chain: Vec<usize> = at.keys().copied().collect();
    chain.sort_by_key(|i| at[i].0);
    let chain_hap: Vec<Vec<Vec<usize>>> = chain.iter().map(|&i| hap[i].clone()).collect();
    let chain_obs: Vec<Option<usize>> = chain.iter().map(|&i| observed[i]).collect();
    let (without, path) = mosaic(&chain_hap, &chain_obs, runs.len());
    let mut weights = vec![0.0; runs.len()];
    for p in &without {
        for (s, x) in p.iter().enumerate() {
            weights[s] += x / without.len() as f64;
        }
    }
    let on_chain: BTreeMap<usize, usize> = chain.iter().enumerate().map(|(t, &i)| (i, t)).collect();
    let label = |s: &Site, a: usize| s.alleles[a].map_or("-".to_string(), |o| g.label(o));
    let mut tsv = String::from(
        "site\tbackbone_step\tcalled\tpredicted\tprobability\tclosest_run\timputed\tagrees\n",
    );
    let (mut scored, mut agree, mut imputed, mut imputed_p) = (0usize, 0usize, 0usize, 0.0);
    // each site's allele for the VCF: (allele, imputed?, probability when imputed)
    let mut chosen: Vec<Option<(usize, bool, f64)>> = Vec::new();
    for (i, s) in sites.iter().enumerate() {
        let t = on_chain.get(&i).copied();
        let w = t.map_or(&weights, |t| &without[t]);
        let pred = predict(&hap[i], w, s.alleles.len());
        let agrees = match (observed[i], pred) {
            (Some(o), Some((a, _))) => {
                scored += 1;
                agree += usize::from(o == a);
                if o == a { "yes" } else { "no" }
            }
            _ => "-",
        };
        let is_imputed = observed[i].is_none() && pred.is_some();
        if is_imputed {
            imputed += 1;
            imputed_p += pred.unwrap().1;
        }
        chosen.push(match (observed[i], pred) {
            (Some(o), _) => Some((o, false, 1.0)),
            (None, Some((a, p))) => Some((a, true, p)),
            _ => None,
        });
        writeln!(
            tsv,
            "{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{agrees}",
            t.map_or("-".to_string(), |t| at[&chain[t]].0.to_string()),
            observed[i].map_or("-".to_string(), |a| label(s, a)),
            pred.map_or("-".to_string(), |(a, _)| label(s, a)),
            pred.map_or("-".to_string(), |(_, p)| format!("{p:.3}")),
            t.map_or("-".to_string(), |t| runs[path[t]].clone()),
            if is_imputed { "yes" } else { "no" },
        )?;
    }
    fs::write(out.join("mosaic.tsv"), tsv)?;
    let records = write_vcf(g, sites, walks, &backbone, &at, &chosen, calls, out)?;
    // the Viterbi stretches along the backbone
    let mut stretches = Vec::new();
    let mut start = 0;
    for t in 1..=path.len() {
        if t == path.len() || path[t] != path[start] {
            stretches.push(json!({"run": runs[path[start]], "first_step": at[&chain[start]].0,
                "last_step": at[&chain[t - 1]].0, "sites": t - start}));
            start = t;
        }
    }
    let mut ranked: Vec<(String, f64)> = runs.iter().cloned().zip(weights.iter().copied()).collect();
    ranked.sort_by(|a, b| b.1.total_cmp(&a.1));
    Ok(json!({"runs": ranked.iter().map(|(r, w)| json!({"run": r, "weight": w})).collect::<Vec<_>>(),
        "backbone": backbone, "sites_on_backbone": chain.len(), "stretches": stretches,
        "called_sites_predicted_without_their_call": scored, "of_them_agreeing": agree,
        "agreement": if scored > 0 { agree as f64 / scored as f64 } else { 0.0 },
        "imputed_sites": imputed, "imputed_mean_probability": if imputed > 0 { imputed_p / imputed as f64 } else { 0.0 },
        "mosaic": out.join("mosaic.tsv"), "vcf": out.join("genotypes.vcf"), "vcf_records": records,
        "basis": "Li-Stephens mosaic over the reference runs (mismatch 0.01, switch 0.001 per site); imputed alleles are the runs', not the reads'"}))
}

/// The genotypes as a VCF on the backbone walk (`write_mosaic`), haploid: every site where the
/// individual's allele (called, else imputed) or a minor allele differs from the walk's. Alleles
/// are spelled in the walk's direction after its base before the site (POS); INFO gives the
/// alleles' read shares (AF), IMP and PP for an imputed allele, and the site's index (SITE).
#[allow(clippy::too_many_arguments)]
fn write_vcf(
    g: &Graph,
    sites: &[Site],
    walks: &[(String, Vec<Oriented>)],
    backbone: &str,
    at: &BTreeMap<usize, (usize, usize, bool)>,
    chosen: &[Option<(usize, bool, f64)>],
    calls: &[Vec<(&'static str, Option<f64>)>],
    out: &Path,
) -> Result<usize> {
    let Some((_, steps)) = walks.iter().find(|(n, _)| n == backbone) else {
        fs::write(out.join("genotypes.vcf"), "##fileformat=VCFv4.2\n")?;
        return Ok(0);
    };
    // where each step ends on the walk
    let mut end = Vec::with_capacity(steps.len());
    let mut len = 0usize;
    for &o in steps {
        len += g.seqs[o.0].len();
        end.push(len);
    }
    let mut records: Vec<(usize, String)> = Vec::new();
    for (&i, &(p, b, rev)) in at {
        let Some((a, imp, prob)) = chosen[i] else {
            continue;
        };
        let s = &sites[i];
        let alts: Vec<usize> = (0..s.alleles.len())
            .filter(|&j| j != b && (j == a || calls[i][j].0 == "minor" || calls[i][j].0 == "major"))
            .collect();
        if alts.is_empty() {
            continue;
        }
        let pad = *oseq(g, steps[p]).last().unwrap();
        let spelled = |j: usize| -> String {
            let x = if rev { s.alleles[j].map(flip) } else { s.alleles[j] };
            let mut v = vec![pad];
            if let Some(x) = x {
                v.extend(oseq(g, x));
            }
            String::from_utf8_lossy(&v).to_ascii_uppercase()
        };
        let gt = if a == b { 0 } else { 1 + alts.iter().position(|&j| j == a).unwrap() };
        let af: Vec<String> = alts
            .iter()
            .map(|&j| calls[i][j].1.map_or(".".to_string(), |f| format!("{f:.3}")))
            .collect();
        let mut info = format!("SITE={i};AF={}", af.join(","));
        if imp {
            write!(info, ";IMP;PP={prob:.3}")?;
        }
        records.push((
            end[p],
            format!(
                "{backbone}\t{}\t.\t{}\t{}\t.\tPASS\t{info}\tGT\t{gt}",
                end[p],
                spelled(b),
                alts.iter().map(|&j| spelled(j)).collect::<Vec<_>>().join(",")
            ),
        ));
    }
    records.sort();
    let mut vcf = String::from("##fileformat=VCFv4.2\n");
    writeln!(vcf, "##contig=<ID={backbone},length={len}>")?;
    vcf.push_str("##INFO=<ID=SITE,Number=1,Type=Integer,Description=\"Site of genotypes.tsv\">\n");
    vcf.push_str("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Read share of each ALT allele\">\n");
    vcf.push_str("##INFO=<ID=IMP,Number=0,Type=Flag,Description=\"Allele imputed from the reference runs, not called from reads\">\n");
    vcf.push_str("##INFO=<ID=PP,Number=1,Type=Float,Description=\"Probability of the imputed allele\">\n");
    vcf.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Haploid genotype: the individual's allele\">\n");
    vcf.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n");
    for (_, r) in &records {
        writeln!(vcf, "{r}")?;
    }
    fs::write(out.join("genotypes.vcf"), vcf)?;
    Ok(records.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_mosaic_follows_the_closest_run_and_predicts_a_missing_call() {
        // three runs over ten sites of two alleles; the individual is run 1 for the first half
        // and run 2 after, and site 3 is not called
        let runs_alleles = [
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 1, 0, 1, 0, 0, 1, 1],
            [0, 1, 1, 0, 1, 0, 1, 1, 0, 1],
        ];
        let hap: Vec<Vec<Vec<usize>>> = (0..10).map(|t| (0..3).map(|r| vec![runs_alleles[r][t]]).collect()).collect();
        let truth: Vec<usize> = (0..10).map(|t| if t < 5 { runs_alleles[1][t] } else { runs_alleles[2][t] }).collect();
        let obs: Vec<Option<usize>> = (0..10).map(|t| (t != 3).then_some(truth[t])).collect();
        let (without, path) = mosaic(&hap, &obs, 3);
        assert_eq!(&path[..3], &[1, 1, 1]);
        assert_eq!(&path[7..], &[2, 2, 2]);
        let (a, p) = predict(&hap[3], &without[3], 2).unwrap();
        assert_eq!(a, truth[3]);
        assert!(p > 0.9, "{p}");
    }

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
        // on run one's walk (it passes both sites): the SNP differs (C for A); the insertion is one's
        let vcf = fs::read_to_string(dir.join("genotypes.vcf")).unwrap();
        let records: Vec<&str> = vcf.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(records.len(), 1, "{vcf}");
        let f: Vec<&str> = records[0].split('\t').collect();
        assert_eq!((f[0], f[1], f[9]), ("one#1#mt", "300", "1"), "{vcf}");
        let last = l.as_bytes()[299] as char;
        assert_eq!((f[3], f[4]), (format!("{last}A").as_str(), format!("{last}C").as_str()), "{vcf}");
        fs::remove_dir_all(dir).unwrap();
    }
}
