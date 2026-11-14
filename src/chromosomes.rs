//! How many chromosomes an organelle's molecules are, and which molecules are conformations of
//! the same one (what recombination across a repeat turns into each other).
//!
//! `linearize` writes the fewest molecules that explain the graph, and that count is no count of
//! chromosomes: the Trifolium striatum mitochondrion is two circles (242,524 and 41,849 bp) that
//! two repeats with reads for both of their pairings (134:140 and 147:161) turn into each other,
//! one chromosome; Nipponbare's three molecules are one chromosome too. The calls:
//!
//! 1. sequences joined by a read-supported link are one chromosome;
//! 2. the molecules of one chromosome are its conformations; each has a share, the reads of the
//!    pairing it takes at the least supported repeat it crosses;
//! 3. molecules are different chromosomes when no read-supported link joins them, or when every
//!    read across a repeat they share takes the pairing that keeps them apart, enough reads for
//!    the other pairing's share to be below `MAX_MISSING_SHARE` at 95%;
//! 4. a repeat they share that no read crosses, or too few, leaves a range of counts;
//! 5. a component far shallower than the main one (`SHALLOW`) without sequence of the organelle's
//!    seeds is possibly nuclear or contamination, and not counted, when its depth is the whole
//!    library's (`--depths`);
//! 6. any other component (none is joined to the main one) at `SOURCE_FOLD` times the main one's
//!    depth or more, or at that share of it or less, is of a source still to be found
//!    (`source_pending`: another mitochondrial type, a high-copy molecule of its own, or
//!    contamination): not counted, only in the range's maximum. In the kept reads' depth this also
//!    takes a component an incomplete recruitment left shallow: Glycine max's 59.5 kb component is
//!    at 37% of the main one's depth in the kept reads, and at 1,253x against 536x in the whole
//!    library; its 62,661 bp circle (442x) is counted.
//!
//! The molecules and the graph do not change; their FASTA headers and GFA paths are annotated.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use rustc_hash::FxHashSet;
use serde::Serialize;
use serde_json::Value;

use crate::kmer::CanonicalKmers;

/// Fewest reads on a link for it to join sequences (`solve linearize --min-reads`).
const MIN_LINK_READS: u64 = 2;
/// A repeat whose reads all take one pairing keeps apart the molecules that pairing separates
/// when the missing pairing's share is below this at 95%: 1 - 0.05^(1/n) for n reads, so from 59
/// reads on (Nipponbare's u29 and u30, 139 and 187 reads; its u10, 10 reads, would not do).
pub const MAX_MISSING_SHARE: f64 = 0.05;
/// A component under this share of the main component's depth (rule 5).
pub const SHALLOW: f64 = 0.5;
/// A component unjoined to the main one at this many times its depth, or this share of it, is of
/// a source still to be found (rule 6).
pub const SOURCE_FOLD: f64 = 2.0;
/// Seed sequence a component must hold, in bp of reads-free seed k-mer blocks, to count as the
/// organelle's own (rule 5): a gene or two of the land-plant seeds.
const MIN_SEED_BP: usize = 300;
const SEED_K: usize = 21;
/// A short-read run does not deliver a molecule whose every component is free of the main one,
/// holds under `MIN_SEED_BP` of seed sequence, closes no circle and is shorter than this many
/// reads (`withhold`): ZM4's read-pair recruitment left 175 pieces of 60-683 bp (20.9 kb, 166 of
/// them on the alfalfa nuclear genome, 155 at five copies or more) in two components unjoined to
/// its mitochondrion, at 28-141x against its 123x, and five molecules of 226-583 bp of them were
/// delivered beside its 306,690 bp one. Depth does not tell them: a nuclear family in the
/// organelle's depth band stands at its depth. A small circle of a multipartite mitochondrion
/// without seed genes is kept (it closes), as is any free component of a few reads' length.
pub const WITHHOLD_READS: usize = 5;

/// One-sided 95% upper bound of a share seen 0 times in `n`.
pub fn missing_share_upper(n: u64) -> f64 {
    if n == 0 {
        1.0
    } else {
        1.0 - 0.05f64.powf(1.0 / n as f64)
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct RepeatCall {
    pub segment: String,
    pub length: usize,
    /// Reads for (entry0 -> exit0 + entry1 -> exit1) and (entry0 -> exit1 + entry1 -> exit0);
    /// none when no read crosses it.
    pub pairings: Option<(u64, u64)>,
    /// interconvertible, apart, too_few_reads or not_crossed
    pub call: &'static str,
    /// 95% upper bound of the missing pairing's share, when one pairing has no reads.
    pub missing_share_upper: Option<f64>,
    /// The molecules (FASTA order, from 1) that cross it.
    pub molecules: Vec<usize>,
}

#[derive(Debug, Clone, Serialize)]
pub struct Conformation {
    /// Its molecule in molecules.fasta (from 1) and its name there.
    pub molecule: usize,
    pub name: String,
    pub length: usize,
    pub circular: bool,
    /// Reads of the pairing it takes at the least supported repeat it crosses, over the reads of
    /// both pairings there; 1 when it crosses none with reads.
    pub share: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct Chromosome {
    pub index: usize,
    pub conformations: Vec<Conformation>,
    pub segments: usize,
    pub bp: usize,
    pub depth: f64,
    pub seed_bp: Option<usize>,
    /// Its segments in the assembly graph (`withhold`).
    #[serde(skip)]
    pub segment_names: Vec<String>,
    /// counted, or source_pending (rule 6)
    pub status: &'static str,
    /// Why it is source_pending.
    pub reason: Option<String>,
    pub repeats: Vec<RepeatCall>,
}

#[derive(Debug, Clone, Serialize)]
pub struct Excluded {
    pub segments: Vec<String>,
    pub bp: usize,
    pub depth: f64,
    pub seed_bp: Option<usize>,
    /// possibly_nuclear_or_contamination (rule 5), or unplaced (no molecule)
    pub reason: &'static str,
}

/// A molecule not delivered (`withhold`), written to molecules.withheld.fasta.
#[derive(Debug, Clone, Serialize)]
pub struct Withheld {
    /// Its place in molecules.fasta before (from 1), its name and length.
    pub molecule: usize,
    pub name: String,
    pub length: usize,
    /// The components holding it (chromosome indices of chromosomes.json) and their segments
    /// in the assembly graph (assembly.gfa; their sequence in components.withheld.fasta).
    pub components: Vec<usize>,
    pub segments: Vec<String>,
    /// The four criteria, each as it stands: joined to the main component (no), seed sequence
    /// against `min_seed_bp`, a circle among its molecules (no), bases against `max_bp`.
    pub joined_to_main: bool,
    pub seed_bp: Option<usize>,
    pub min_seed_bp: usize,
    pub circular: bool,
    pub components_bp: usize,
    pub max_bp: usize,
    pub reason: &'static str,
}

#[derive(Debug, Clone, Serialize)]
pub struct ChromosomeReport {
    /// An estimate for the represented graph, not a certified physical chromosome count.
    pub count: Option<usize>,
    /// Bounds for represented graph groups, not bounds on unobserved biological molecules.
    pub range: (usize, usize),
    /// Disconnected components contain linear paths that may be assembly fragments.
    pub fragmented_components: bool,
    /// Components of a source still to be found (rule 6), in `range`'s maximum only.
    pub source_pending: usize,
    pub reasons: Vec<String>,
    /// "library" (`--depths`, whole-library k-mer depth) or "kept_reads" (assembly dp)
    pub depth_source: &'static str,
    pub max_missing_share: f64,
    pub shallow: f64,
    pub source_fold: f64,
    pub chromosomes: Vec<Chromosome>,
    pub excluded: Vec<Excluded>,
    /// Molecules not delivered (`withhold`); none unless a run asks.
    pub withheld: Vec<Withheld>,
}

struct Molecule {
    name: String,
    path: Vec<(String, bool)>,
    length: usize,
    circular: bool,
}

fn oriented(s: &str) -> (String, bool) {
    let forward = !s.ends_with('-');
    (s.trim_end_matches(['+', '-']).to_string(), forward)
}

fn flip(o: &(String, bool)) -> (String, bool) {
    (o.0.clone(), !o.1)
}

/// Which pairing (0 or 1) a molecule takes at `repeat` (in the orientation evidence reports it)
/// at each place its path crosses it.
fn pairings_taken(m: &Molecule, repeat: &(String, bool), entries: &[(String, bool)], exits: &[(String, bool)]) -> Vec<usize> {
    let n = m.path.len();
    let mut out = Vec::new();
    for i in 0..n {
        if m.path[i].0 != repeat.0 || (!m.circular && (i == 0 || i + 1 == n)) {
            continue;
        }
        let (prev, next) = (&m.path[(i + n - 1) % n], &m.path[(i + 1) % n]);
        let (entry, exit) = if m.path[i].1 == repeat.1 {
            (prev.clone(), next.clone())
        } else {
            (flip(next), flip(prev))
        };
        let e = entries.iter().position(|x| *x == entry);
        let x = exits.iter().position(|y| *y == exit);
        if let (Some(e), Some(x)) = (e, x) {
            out.push(usize::from(e != x));
        }
    }
    out
}

/// Bp of `seq` in runs of k-mers found in `seeds`.
fn seed_bp(seq: &[u8], seeds: &FxHashSet<u64>) -> usize {
    let (mut bp, mut run_end) = (0usize, 0usize);
    for (pos, kmer) in CanonicalKmers::new(seq, SEED_K) {
        if seeds.contains(&kmer) {
            let end = pos + SEED_K;
            bp += end - run_end.max(pos);
            run_end = end;
        }
    }
    bp
}

/// The calls of one organelle's output: its graph, evidence and linearization (FASTA order of
/// molecules), optional whole-library depths per segment and optional seed k-mers.
pub fn call(
    gfa: &Path,
    evidence: &Value,
    linearization: &Value,
    depths: Option<&BTreeMap<String, f64>>,
    seeds: Option<&FxHashSet<u64>>,
) -> Result<ChromosomeReport> {
    crate::repeat_depth::validate_report_graph(gfa, evidence)?;
    let text = fs::read_to_string(gfa).with_context(|| format!("cannot read {}", gfa.display()))?;
    let mut length = BTreeMap::new();
    let mut dp = BTreeMap::new();
    let mut seq_seed = BTreeMap::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.first() != Some(&"S") || f.len() < 3 {
            continue;
        }
        length.insert(f[1].to_string(), f[2].len());
        if let Some(d) = f[3..].iter().find_map(|t| t.strip_prefix("dp:f:")?.parse::<f64>().ok()) {
            dp.insert(f[1].to_string(), d);
        }
        if let Some(s) = seeds {
            seq_seed.insert(f[1].to_string(), seed_bp(f[2].as_bytes(), s));
        }
    }
    let depth_of = |s: &str| match depths {
        Some(d) => d.get(s).copied(),
        None => dp.get(s).copied(),
    };
    // rule 1: components of read-supported links
    let names: Vec<String> = length.keys().cloned().collect();
    let index: BTreeMap<&str, usize> = names.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let mut parent: Vec<usize> = (0..names.len()).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    for l in evidence["links"].as_array().into_iter().flatten() {
        if l["reads"].as_u64().unwrap_or(0) < MIN_LINK_READS {
            continue;
        }
        let a = oriented(l["from"].as_str().unwrap_or_default()).0;
        let b = oriented(l["to"].as_str().unwrap_or_default()).0;
        if let (Some(&a), Some(&b)) = (index.get(a.as_str()), index.get(b.as_str())) {
            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
            parent[ra] = rb;
        }
    }
    let mut members: BTreeMap<usize, Vec<String>> = BTreeMap::new();
    for (i, s) in names.iter().enumerate() {
        members.entry(find(&mut parent, i)).or_default().push(s.clone());
    }
    let mut molecules = Vec::new();
    for c in linearization["components"].as_array().into_iter().flatten() {
        for m in c["best"]["molecules"].as_array().into_iter().flatten() {
            molecules.push(Molecule {
                name: format!("linear.{}", molecules.len() + 1),
                path: m["path"].as_str().unwrap_or_default().split_whitespace().map(oriented).collect(),
                length: m["length"].as_u64().unwrap_or(0) as usize,
                circular: m["circular"].as_bool().unwrap_or(false),
            });
        }
    }
    let repeats: BTreeMap<String, &Value> = evidence["repeats"]
        .as_array()
        .into_iter()
        .flatten()
        .map(|r| (oriented(r["repeat"].as_str().unwrap_or_default()).0, r))
        .collect();
    let repeat_segments: BTreeSet<String> = evidence["repeat_segments"]
        .as_array()
        .into_iter()
        .flatten()
        .filter_map(|s| s.as_str().map(|s| oriented(s).0))
        .chain(repeats.keys().cloned())
        .collect();
    // the components, largest first; the main one is the largest with molecules
    let mut comps: Vec<Vec<String>> = members.into_values().collect();
    comps.sort_by_key(|c| std::cmp::Reverse(c.iter().map(|s| length[s]).sum::<usize>()));
    let weighted = |c: &[String]| {
        let (mut sum, mut bp) = (0.0, 0usize);
        for s in c {
            if let Some(d) = depth_of(s) {
                sum += d * length[s] as f64;
                bp += length[s];
            }
        }
        (bp > 0).then(|| sum / bp as f64)
    };
    let uses = |m: &Molecule, c: &BTreeSet<&str>| m.path.iter().any(|o| c.contains(o.0.as_str()));
    let main_depth = comps
        .iter()
        .find(|c| {
            let set: BTreeSet<&str> = c.iter().map(String::as_str).collect();
            molecules.iter().any(|m| uses(m, &set))
        })
        .and_then(|c| weighted(c));
    let mut report = ChromosomeReport {
        count: None,
        range: (0, 0),
        fragmented_components: false,
        source_pending: 0,
        reasons: Vec::new(),
        depth_source: if depths.is_some() { "library" } else { "kept_reads" },
        max_missing_share: MAX_MISSING_SHARE,
        shallow: SHALLOW,
        source_fold: SOURCE_FOLD,
        chromosomes: Vec::new(),
        excluded: Vec::new(),
        withheld: Vec::new(),
    };
    let (mut lo, mut hi) = (0usize, 0usize);
    let mut main_seen = false;
    for c in &comps {
        let set: BTreeSet<&str> = c.iter().map(String::as_str).collect();
        let bp: usize = c.iter().map(|s| length[s]).sum();
        let depth = weighted(c).unwrap_or(0.0);
        let seed: Option<usize> = seeds.map(|_| c.iter().map(|s| seq_seed[s]).sum());
        let mine: Vec<usize> = (0..molecules.len()).filter(|&i| uses(&molecules[i], &set)).collect();
        if mine.is_empty() {
            report.excluded.push(Excluded {
                segments: c.clone(),
                bp,
                depth,
                seed_bp: seed,
                reason: "unplaced",
            });
            continue;
        }
        let is_main = !main_seen;
        main_seen = true;
        let shallow = main_depth.is_some_and(|m| depth < SHALLOW * m);
        let no_seed = seed.is_some_and(|b| b < MIN_SEED_BP);
        if shallow && no_seed && depths.is_some() {
            report.excluded.push(Excluded {
                segments: c.clone(),
                bp,
                depth,
                seed_bp: seed,
                reason: "possibly_nuclear_or_contamination",
            });
            continue;
        }
        // rules 2-4 over the repeats its molecules share: molecules sharing a repeat with reads for
        // both pairings are one chromosome; one with too few reads, or none, joins them only at the
        // range's minimum; one read for one pairing alone, enough of them, keeps them apart. A
        // molecule sharing no repeat with the others is one chromosome with them (rule 1).
        let mut calls = Vec::new();
        let mut share = vec![1.0f64; molecules.len()];
        let slot: BTreeMap<usize, usize> = mine.iter().enumerate().map(|(j, &i)| (i, j)).collect();
        let mut joined: Vec<usize> = (0..mine.len()).collect();
        let mut joined_or_open = joined.clone();
        let mut involved = vec![false; mine.len()];
        let mut open = 0usize;
        for s in c.iter().filter(|s| repeat_segments.contains(*s)) {
            let crossing: Vec<usize> = mine
                .iter()
                .copied()
                .filter(|&i| molecules[i].path.iter().any(|o| &o.0 == s))
                .collect();
            let shared = crossing.len() > 1;
            let unite = |sets: &mut [&mut Vec<usize>]| {
                for w in crossing.windows(2) {
                    let (a, b) = (slot[&w[0]], slot[&w[1]]);
                    for p in sets.iter_mut() {
                        let (ra, rb) = (find(p, a), find(p, b));
                        p[ra] = rb;
                    }
                }
            };
            if shared {
                for i in &crossing {
                    involved[slot[i]] = true;
                }
            }
            let Some(r) = repeats.get(s) else {
                if shared {
                    open += 1;
                    unite(&mut [&mut joined_or_open]);
                }
                calls.push(RepeatCall {
                    segment: s.clone(),
                    length: length[s],
                    pairings: None,
                    call: "not_crossed",
                    missing_share_upper: None,
                    molecules: crossing.iter().map(|i| i + 1).collect(),
                });
                continue;
            };
            let p = r["pairing_reads"].as_array().map_or((0, 0), |v| {
                (v[0].as_u64().unwrap_or(0), v.get(1).and_then(Value::as_u64).unwrap_or(0))
            });
            let n = p.0 + p.1;
            let upper = (p.0.min(p.1) == 0).then(|| missing_share_upper(n));
            let call = match upper {
                None => "interconvertible",
                Some(_) if n == 0 => "not_crossed",
                Some(u) if u < MAX_MISSING_SHARE => "apart",
                Some(_) => "too_few_reads",
            };
            if shared {
                match call {
                    "apart" => {}
                    "interconvertible" => unite(&mut [&mut joined, &mut joined_or_open]),
                    _ => {
                        open += 1;
                        unite(&mut [&mut joined_or_open]);
                    }
                }
            }
            let rep = oriented(r["repeat"].as_str().unwrap_or_default());
            let side = |key: &str| -> Vec<(String, bool)> {
                r[key].as_array().into_iter().flatten().filter_map(|x| x.as_str().map(oriented)).collect()
            };
            let (entries, exits) = (side("entries"), side("exits"));
            for &i in &crossing {
                for k in pairings_taken(&molecules[i], &rep, &entries, &exits) {
                    let reads = if k == 0 { p.0 } else { p.1 };
                    if n > 0 {
                        share[i] = share[i].min(reads as f64 / n as f64);
                    }
                }
            }
            calls.push(RepeatCall {
                segment: s.clone(),
                length: length[s],
                pairings: Some(p),
                call,
                missing_share_upper: upper,
                molecules: crossing.iter().map(|i| i + 1).collect(),
            });
        }
        // groups of the molecules that share repeats (the others join any of them): with the
        // open repeats joining, the fewest chromosomes; without, the most
        let groups = |p: &mut Vec<usize>| {
            let roots: BTreeSet<usize> =
                (0..mine.len()).filter(|&j| involved[j]).map(|j| find(p, j)).collect();
            roots.len().max(1)
        };
        let at_least = groups(&mut joined_or_open);
        let at_most = groups(&mut joined);
        // rule 6: how far its depth is from the main one's, either way
        let fold = main_depth
            .filter(|&m| m > 0.0 && depth > 0.0)
            .map(|m| (depth / m).max(m / depth));
        let (status, reason) = if !is_main && fold.is_some_and(|f| f >= SOURCE_FOLD) {
            let main = main_depth.unwrap_or(0.0);
            let why = format!(
                "a component of {bp} bp unjoined to the main one is at {depth:.1} against its {main:.1} \
                 ({} depth, {:.2} times): its source is to be found (another mitochondrial type, a \
                 high-copy molecule of its own, or contamination)",
                report.depth_source,
                depth / main,
            );
            report.reasons.push(why.clone());
            report.source_pending += 1;
            ("source_pending", Some(why))
        } else {
            lo += at_least;
            ("counted", None)
        };
        hi += at_most;
        if at_most > at_least {
            report.reasons.push(format!(
                "{} repeat(s) shared by molecules {:?} without enough reads across them",
                open,
                mine.iter().map(|i| i + 1).collect::<Vec<_>>()
            ));
        }
        report.chromosomes.push(Chromosome {
            index: report.chromosomes.len() + 1,
            conformations: mine
                .iter()
                .map(|&i| Conformation {
                    molecule: i + 1,
                    name: molecules[i].name.clone(),
                    length: molecules[i].length,
                    circular: molecules[i].circular,
                    share: share[i],
                })
                .collect(),
            segments: c.len(),
            segment_names: c.clone(),
            bp,
            depth,
            seed_bp: seed,
            status,
            reason,
            repeats: calls,
        });
    }
    // A disconnected linear path may be a piece of a chromosome shared with another component.
    // No edge is not evidence of biological separation: do not count these as established
    // chromosomes merely because incomplete assembly split them. Preserve components and paths.
    let fragmented = report.chromosomes.len() > 1
        && report.chromosomes.iter().any(|c| c.conformations.iter().any(|m| !m.circular));
    if fragmented {
        report.fragmented_components = true;
        report.reasons.push("Disconnected graph components include linear paths; they may be assembly fragments. The range describes represented graph groups, not biological chromosome-count bounds.".into());
        for c in &mut report.chromosomes {
            if c.conformations.iter().any(|m| !m.circular) {
                let why = "unconnected linear assembly path; biological separation is unresolved";
                c.reason = Some(match c.reason.take() {
                    Some(previous) => format!("{previous}; {why}"),
                    None => why.into(),
                });
                c.status = "fragment_unresolved";
            }
        }
        report.range = (usize::from(lo > 0), hi);
        report.count = None;
    } else {
        report.range = (lo, hi);
        report.count = (lo == hi).then_some(lo);
    }
    Ok(report)
}

/// The chromosome and conformation of each molecule (FASTA order, from 1).
fn labels(report: &ChromosomeReport) -> BTreeMap<usize, (usize, usize, f64, &'static str)> {
    let mut out = BTreeMap::new();
    for c in &report.chromosomes {
        for (j, conf) in c.conformations.iter().enumerate() {
            out.insert(conf.molecule, (c.index, j + 1, conf.share, c.status));
        }
    }
    out
}

/// The molecules a short-read run does not deliver (`WITHHOLD_READS`): every component holding
/// one is free of the main component (the first), holds seed sequence and under `MIN_SEED_BP` of
/// it, has no circular molecule and is under `max_bp`. Without seeds nothing is withheld.
pub fn withhold(report: &mut ChromosomeReport, max_bp: usize) {
    let mut of: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (ci, c) in report.chromosomes.iter().enumerate() {
        for m in &c.conformations {
            of.entry(m.molecule).or_default().push(ci);
        }
    }
    let fragment = |ci: usize| {
        let c = &report.chromosomes[ci];
        ci > 0
            && c.seed_bp.is_some_and(|b| b < MIN_SEED_BP)
            && c.conformations.iter().all(|m| !m.circular)
            && c.bp < max_bp
    };
    let mut withheld = Vec::new();
    for (&m, cs) in &of {
        if !cs.iter().all(|&ci| fragment(ci)) {
            continue;
        }
        let conf = report.chromosomes[cs[0]].conformations.iter().find(|c| c.molecule == m).unwrap();
        withheld.push(Withheld {
            molecule: m,
            name: conf.name.clone(),
            length: conf.length,
            components: cs.iter().map(|&ci| report.chromosomes[ci].index).collect(),
            segments: cs.iter().flat_map(|&ci| report.chromosomes[ci].segment_names.iter().cloned()).collect(),
            joined_to_main: false,
            seed_bp: cs.iter().map(|&ci| report.chromosomes[ci].seed_bp).sum(),
            min_seed_bp: MIN_SEED_BP,
            circular: false,
            components_bp: cs.iter().map(|&ci| report.chromosomes[ci].bp).sum(),
            max_bp,
            reason: "free_fragment_without_seed",
        });
    }
    if !withheld.is_empty() {
        report.reasons.push(format!(
            "{} molecule(s) of {} bp in all not delivered: free components without seed sequence, \
             no circle and under {max_bp} bp (molecules.withheld.fasta)",
            withheld.len(),
            withheld.iter().map(|w| w.length).sum::<usize>()
        ));
    }
    report.withheld = withheld;
}

/// Take the withheld molecules' segments (`withhold`) out of the organelle graph at `path`
/// (organelle.gfa, its segments named `<sample>.<segment>[.<piece>]`), with their links and the
/// paths made of them alone; a segment a kept path walks stays. The assembly graph keeps them.
/// The number of segments taken out.
pub fn drop_withheld_gfa(path: &Path, report: &ChromosomeReport) -> Result<usize> {
    let gone: BTreeSet<&str> = report.withheld.iter().flat_map(|w| w.segments.iter().map(String::as_str)).collect();
    if gone.is_empty() {
        return Ok(0);
    }
    let text = fs::read_to_string(path)?;
    let sample = text
        .lines()
        .find(|l| l.starts_with("H\t"))
        .and_then(|h| h.split('\t').find_map(|f| f.strip_prefix("sn:Z:")))
        .unwrap_or_default()
        .to_string();
    let base = |name: &str| -> String {
        let own = name.strip_prefix(&format!("{sample}.")).unwrap_or(name);
        own.split('.').next().unwrap_or(own).to_string()
    };
    let steps = |walk: &str| -> Vec<String> {
        walk.split(',').map(|s| s.trim_end_matches(['+', '-']).to_string()).collect()
    };
    // segments a path that is not all withheld walks
    let mut walked: BTreeSet<String> = BTreeSet::new();
    for line in text.lines().filter(|l| l.starts_with("P\t")) {
        let segs = steps(line.split('\t').nth(2).unwrap_or_default());
        if !segs.iter().all(|s| gone.contains(base(s).as_str())) {
            walked.extend(segs);
        }
    }
    let drop = |name: &str| gone.contains(base(name).as_str()) && !walked.contains(name);
    let (mut out, mut n) = (String::with_capacity(text.len()), 0);
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let keep = match f.first().copied() {
            Some("S") => {
                let d = drop(f.get(1).copied().unwrap_or_default());
                n += d as usize;
                !d
            }
            Some("L") => !drop(f.get(1).copied().unwrap_or_default()) && !drop(f.get(3).copied().unwrap_or_default()),
            Some("P") => !steps(f.get(2).copied().unwrap_or_default()).iter().all(|s| drop(s)),
            _ => true,
        };
        if keep {
            out.push_str(line);
            out.push('\n');
        }
    }
    fs::write(path, out)?;
    Ok(n)
}

/// The withheld molecules' segments of the assembly graph `gfa`, with their component and molecule.
pub fn write_withheld_segments(gfa: &Path, to: &Path, report: &ChromosomeReport) -> Result<()> {
    let mut of: BTreeMap<&str, (usize, usize)> = BTreeMap::new();
    for w in &report.withheld {
        for s in &w.segments {
            let c = report.chromosomes.iter().find(|c| c.segment_names.contains(s)).map_or(0, |c| c.index);
            of.entry(s.as_str()).or_insert((c, w.molecule));
        }
    }
    let mut out = String::new();
    for line in fs::read_to_string(gfa)?.lines().filter(|l| l.starts_with("S\t")) {
        let f: Vec<&str> = line.split('\t').collect();
        if let (Some(name), Some(seq)) = (f.get(1), f.get(2)) {
            if let Some((c, m)) = of.get(name) {
                out.push_str(&format!(">{name} component={c} molecule={m} length={}\n{seq}\n", seq.len()));
            }
        }
    }
    fs::write(to, out)?;
    Ok(())
}

/// Move the withheld molecules (`withhold`) from `path` to `to`, in their order.
pub fn split_withheld(path: &Path, to: &Path, report: &ChromosomeReport) -> Result<()> {
    let gone: BTreeSet<usize> = report.withheld.iter().map(|w| w.molecule).collect();
    let text = fs::read_to_string(path)?;
    let (mut kept, mut moved) = (String::with_capacity(text.len()), String::new());
    let mut n = 0;
    for line in text.lines() {
        if line.starts_with('>') {
            n += 1;
        }
        let out = if gone.contains(&n) { &mut moved } else { &mut kept };
        out.push_str(line);
        out.push('\n');
    }
    fs::write(path, kept)?;
    fs::write(to, moved)?;
    Ok(())
}

/// Append `chromosome=`, `conformation=` and `share=` to the molecules' FASTA headers.
pub fn annotate_fasta(path: &Path, report: &ChromosomeReport) -> Result<()> {
    let text = fs::read_to_string(path)?;
    let labels = labels(report);
    let mut out = String::with_capacity(text.len() + 64 * labels.len());
    let mut n = 0;
    for line in text.lines() {
        if line.starts_with('>') {
            // Replace our own old annotations when re-running this command, while preserving
            // other header fields and conservative provenance from earlier stages.
            let fragment = labels.get(&(n + 1)).is_some_and(|label| label.3 == "fragment_unresolved");
            let clean = line.split_whitespace().filter(|field| {
                !(fragment && *field == "reads_decide=false")
                    && !["chromosome=", "conformation=", "share=", "chromosome_status=", "graph_component_group="]
                        .iter().any(|prefix| field.starts_with(prefix))
            }).collect::<Vec<_>>().join(" ");
            out.push_str(&clean);
            n += 1;
            if let Some((c, k, share, status)) = labels.get(&n) {
                if *status == "fragment_unresolved" {
                    out.push_str(&format!(" graph_component_group={c} chromosome_status=fragment_unresolved"));
                    out.push_str(" reads_decide=false");
                } else {
                    out.push_str(&format!(" chromosome={c} conformation={k} share={share:.3}"));
                    if *status != "counted" {
                        out.push_str(&format!(" chromosome_status={status}"));
                    }
                }
            }
        } else {
            out.push_str(line);
        }
        out.push('\n');
    }
    fs::write(path, out)?;
    Ok(())
}

/// Tag the paths of the representative solution (`<sample>#1#<prefix>.<n>`, molecule n) in
/// organelle.gfa with `CH:i:` (chromosome) and `CF:i:` (conformation); names stay as they are.
pub fn annotate_gfa(path: &Path, report: &ChromosomeReport) -> Result<()> {
    let text = fs::read_to_string(path)?;
    let labels = labels(report);
    let mut out = String::with_capacity(text.len());
    for line in text.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        let mut new_label = None;
        if fields.first() == Some(&"P") {
            let name = fields.get(1).copied().unwrap_or_default();
            let mut parts = name.rsplitn(2, '#');
            let last = parts.next().unwrap_or_default();
            let solution = parts.next().and_then(|rest| rest.rsplit('#').next());
            let n = last.rsplit('.').next().and_then(|x| x.parse::<usize>().ok());
            if let (Some("1"), Some(n)) = (solution, n) {
                new_label = labels.get(&n);
            }
        }
        if let Some((c, k, _, status)) = new_label {
            let clean: Vec<&str> = fields.iter().copied().enumerate().filter_map(|(i, f)| {
                (i < 4 || (!["CH:", "CF:", "FG:"].iter().any(|prefix| f.starts_with(prefix))
                    && !(*status == "fragment_unresolved" && f == "up:Z:fragment_unresolved"))).then_some(f)
            }).collect();
            out.push_str(&clean.join("\t"));
            if *status == "fragment_unresolved" {
                out.push_str(&format!("\tFG:i:{c}"));
                if !clean.iter().skip(4).any(|f| f.starts_with("up:")) {
                    out.push_str("\tup:Z:fragment_unresolved");
                }
            } else {
                out.push_str(&format!("\tCH:i:{c}\tCF:i:{k}"));
            }
        } else {
            out.push_str(line);
        }
        out.push('\n');
    }
    fs::write(path, out)?;
    Ok(())
}

/// Seed k-mers of the organelle's seed files.
pub fn seed_kmers(paths: &[std::path::PathBuf]) -> Result<FxHashSet<u64>> {
    let mut set = FxHashSet::default();
    for p in paths {
        let mut reader = needletail::parse_fastx_file(p)
            .with_context(|| format!("cannot read seeds {}", p.display()))?;
        while let Some(rec) = reader.next() {
            let rec = rec?;
            for (_, k) in CanonicalKmers::new(&rec.seq(), SEED_K) {
                set.insert(k);
            }
        }
    }
    Ok(set)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn gfa(dir: &Path, segs: &[(&str, usize, f64)]) -> std::path::PathBuf {
        let mut text = String::from("H\tVN:Z:1.0\n");
        for (i, (name, len, dp)) in segs.iter().enumerate() {
            let seq: String = (0..*len).map(|j| b"ACGT"[(j * 7 + i * 13 + j / 5) % 4] as char).collect();
            text.push_str(&format!("S\t{name}\t{seq}\tdp:f:{dp}\n"));
        }
        let p = dir.join("assembly.gfa");
        fs::write(&p, text).unwrap();
        p
    }

    fn dir(tag: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("ovasm-chrom-{tag}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&d);
        fs::create_dir_all(&d).unwrap();
        d
    }

    fn link(a: &str, b: &str, reads: u64) -> Value {
        json!({"from": a, "to": b, "reads": reads})
    }

    #[test]
    fn missing_share_bound_needs_59_reads_for_5_percent() {
        assert!(missing_share_upper(59) < MAX_MISSING_SHARE);
        assert!(missing_share_upper(58) >= MAX_MISSING_SHARE);
        assert!((missing_share_upper(10) - 0.259).abs() < 0.001);
    }

    #[test]
    fn two_circles_a_repeat_turns_into_each_other_are_one_chromosome() {
        // Trifolium striatum's mitochondrion: two circles across repeats u4 and u5, both pairings
        // read (134:140, 147:161)
        let d = dir("trifolium");
        let g = gfa(&d, &[("u0", 900, 150.0), ("u1", 700, 120.0), ("u2", 600, 130.0), ("u3", 600, 130.0), ("u4", 300, 220.0), ("u5", 300, 230.0)]);
        let ev = json!({
            "links": [link("u0+", "u4-", 100), link("u4-", "u3+", 100), link("u3+", "u5-", 100), link("u5-", "u0+", 100),
                      link("u1+", "u5+", 90), link("u5+", "u2-", 90), link("u2-", "u4+", 90), link("u4+", "u1+", 90)],
            "repeats": [
                {"repeat": "u4+", "entries": ["u2-", "u3-"], "exits": ["u0-", "u1+"], "pairing_reads": [134, 140]},
                {"repeat": "u5+", "entries": ["u1+", "u3+"], "exits": ["u0+", "u2-"], "pairing_reads": [147, 161]}],
            "repeat_segments": ["u4", "u5"]});
        let lin = json!({"components": [{"best": {"molecules": [
            {"path": "u0+ u4- u3+ u5-", "circular": true, "length": 242524},
            {"path": "u1+ u5+ u2- u4+", "circular": true, "length": 41849}]}}]});
        let r = call(&g, &ev, &lin, None, None).unwrap();
        assert_eq!(r.count, Some(1));
        assert_eq!(r.chromosomes[0].conformations.len(), 2);
        assert!(r.chromosomes[0].repeats.iter().all(|c| c.call == "interconvertible"));
        fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn unlinked_components_and_a_repeat_read_one_way_only_are_separate() {
        let d = dir("apart");
        let g = gfa(&d, &[("a", 900, 100.0), ("b", 900, 100.0), ("r", 300, 200.0), ("c", 800, 95.0)]);
        // a and b share repeat r, all 120 reads across it close each on itself; c is unlinked
        let base = json!({
            "links": [link("a+", "r+", 60), link("r+", "a+", 60), link("b+", "r+", 60), link("r+", "b+", 60), link("c+", "c+", 40)],
            "repeat_segments": ["r"]});
        let mut ev = base.clone();
        ev["repeats"] = json!([{"repeat": "r+", "entries": ["a+", "b+"], "exits": ["a+", "b+"], "pairing_reads": [120, 0]}]);
        let lin = json!({"components": [
            {"best": {"molecules": [{"path": "a+ r+", "circular": true, "length": 1200}, {"path": "b+ r+", "circular": true, "length": 1200}]}},
            {"best": {"molecules": [{"path": "c+", "circular": true, "length": 800}]}}]});
        let r = call(&g, &ev, &lin, None, None).unwrap();
        assert_eq!(r.count, Some(3));
        // ten reads only: a range
        ev["repeats"][0]["pairing_reads"] = json!([10, 0]);
        let r = call(&g, &ev, &lin, None, None).unwrap();
        assert_eq!((r.count, r.range), (None, (2, 3)));
        // no read across it: a range as well
        let r = call(&g, &base, &lin, None, None).unwrap();
        assert_eq!((r.count, r.range), (None, (2, 3)));
        fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn an_uncrossed_repeat_between_molecules_already_one_chromosome_changes_nothing() {
        // Nipponbare: molecules 1 and 3 share u31 (158:187) and u0, which no read crosses;
        // molecule 2 shares no repeat with them and is one chromosome with them by its links
        let d = dir("nipponbare");
        let g = gfa(&d, &[("a", 900, 70.0), ("b", 900, 70.0), ("c", 900, 70.0), ("r", 300, 140.0), ("q", 300, 140.0)]);
        let ev = json!({
            "links": [link("a+", "r+", 50), link("r+", "c+", 50), link("c+", "q+", 50), link("q+", "a+", 50), link("b+", "a+", 30)],
            "repeats": [{"repeat": "r+", "entries": ["a+", "q+"], "exits": ["c+", "q+"], "pairing_reads": [158, 187]}],
            "repeat_segments": ["r", "q"]});
        let lin = json!({"components": [{"best": {"molecules": [
            {"path": "a+ r+ q+", "circular": true, "length": 1500},
            {"path": "b+", "circular": false, "length": 900},
            {"path": "c+ q+ r+", "circular": true, "length": 1500}]}}]});
        let r = call(&g, &ev, &lin, None, None).unwrap();
        assert_eq!(r.count, Some(1));
        assert_eq!(r.chromosomes[0].conformations.len(), 3);
        fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn disconnected_linear_fragments_are_not_counted_as_chromosomes() {
        let d = dir("fragmented");
        let g = gfa(&d, &[("a", 1200, 100.0), ("b", 1000, 100.0), ("c", 900, 100.0)]);
        let ev = json!({"links": [], "repeats": [], "repeat_segments": []});
        let lin = json!({"components": [
            {"best": {"molecules": [{"path": "a+", "circular": false, "length": 1200}]}},
            {"best": {"molecules": [{"path": "b+", "circular": false, "length": 1000}]}},
            {"best": {"molecules": [{"path": "c+", "circular": false, "length": 900}]}}]});
        let r = call(&g, &ev, &lin, None, None).unwrap();
        assert_eq!(r.count, None);
        assert_eq!(r.range, (1, 3));
        assert!(r.fragmented_components);
        assert!(r.chromosomes.iter().all(|c| c.status == "fragment_unresolved"));
        assert!(r.reasons.iter().any(|s| s.contains("not biological chromosome-count bounds")));
        let fa = d.join("molecules.fa");
        fs::write(&fa, ">one chromosome=7 conformation=2 share=0.500\nACGT\n>two\nACGT\n>three\nACGT\n").unwrap();
        annotate_fasta(&fa, &r).unwrap();
        annotate_fasta(&fa, &r).unwrap();
        let text = fs::read_to_string(fa).unwrap();
        assert!(text.contains("graph_component_group=1 chromosome_status=fragment_unresolved reads_decide=false"));
        assert!(!text.contains(" chromosome="));
        let paths = d.join("paths.gfa");
        fs::write(&paths, "P\tsample#1#mt.1\ta+\t*\tCH:i:7\tCF:i:2\n").unwrap();
        annotate_gfa(&paths, &r).unwrap();
        annotate_gfa(&paths, &r).unwrap();
        let text = fs::read_to_string(paths).unwrap();
        assert!(text.contains("FG:i:1\tup:Z:fragment_unresolved"));
        assert!(!text.contains("CH:i:"));
        assert_eq!(text.matches("up:Z:fragment_unresolved").count(), 1);
        fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn a_component_far_from_the_main_ones_depth_is_excluded_or_its_source_pending() {
        let d = dir("shallow");
        let g = gfa(&d, &[("m", 2000, 130.0), ("x", 800, 40.0)]);
        let ev = json!({"links": [link("m+", "m+", 50), link("x+", "x+", 20)], "repeats": [], "repeat_segments": []});
        let lin = json!({"components": [
            {"best": {"molecules": [{"path": "m+", "circular": true, "length": 2000}]}},
            {"best": {"molecules": [{"path": "x+", "circular": true, "length": 800}]}}]});
        let seeds = FxHashSet::default();
        // the kept reads' depth (31% of the main one's): its source is pending, a range
        let r = call(&g, &ev, &lin, None, Some(&seeds)).unwrap();
        assert_eq!((r.count, r.range, r.source_pending), (None, (1, 2), 1));
        assert_eq!(r.chromosomes[1].status, "source_pending");
        assert!(r.chromosomes[1].reason.as_deref().is_some_and(|w| w.contains("kept_reads")));
        assert_eq!((r.chromosomes[0].status, r.chromosomes[0].reason.as_deref()), ("counted", None));
        // the whole library's, no seed sequence: not counted
        let lib: BTreeMap<String, f64> = [("m".into(), 480.0), ("x".into(), 60.0)].into();
        let r = call(&g, &ev, &lin, Some(&lib), Some(&seeds)).unwrap();
        assert_eq!(r.count, Some(1));
        assert_eq!(r.excluded[0].reason, "possibly_nuclear_or_contamination");
        // as Glycine max's 59.5 kb component in the library (1,253x against 536x): pending
        let lib: BTreeMap<String, f64> = [("m".into(), 536.0), ("x".into(), 1253.0)].into();
        let r = call(&g, &ev, &lin, Some(&lib), Some(&seeds)).unwrap();
        assert_eq!((r.count, r.range, r.source_pending), (None, (1, 2), 1));
        assert_eq!(r.chromosomes[1].status, "source_pending");
        // as its 62,661 bp circle (442x against 536x): within twice, counted
        let lib: BTreeMap<String, f64> = [("m".into(), 536.0), ("x".into(), 442.0)].into();
        let r = call(&g, &ev, &lin, Some(&lib), Some(&seeds)).unwrap();
        assert_eq!((r.count, r.source_pending), (Some(2), 0));
        fs::remove_dir_all(&d).ok();
    }

    #[test]
    fn a_free_fragment_without_seed_is_withheld_and_a_free_circle_is_not() {
        let d = dir("withhold");
        let g = gfa(&d, &[("m", 2000, 130.0), ("x", 300, 120.0), ("y", 5000, 125.0)]);
        let ev = json!({"links": [link("m+", "m+", 50), link("y+", "y+", 30)], "repeats": [], "repeat_segments": []});
        let lin = json!({"components": [
            {"best": {"molecules": [{"path": "m+", "circular": true, "length": 2000}]}},
            {"best": {"molecules": [{"path": "x+", "circular": false, "length": 300}]}},
            {"best": {"molecules": [{"path": "y+", "circular": true, "length": 5000}]}}]});
        let seeds = FxHashSet::default();
        let mut r = call(&g, &ev, &lin, None, Some(&seeds)).unwrap();
        // nothing unless asked
        assert!(r.withheld.is_empty());
        withhold(&mut r, 750);
        assert_eq!(r.withheld.iter().map(|w| w.molecule).collect::<Vec<_>>(), vec![2]);
        assert_eq!(r.withheld[0].reason, "free_fragment_without_seed");
        // a fragment longer than the bound is delivered
        let mut r = call(&g, &ev, &lin, None, Some(&seeds)).unwrap();
        withhold(&mut r, 300);
        assert!(r.withheld.is_empty());
        // without seeds nothing is judged
        let mut r = call(&g, &ev, &lin, None, None).unwrap();
        withhold(&mut r, 750);
        assert!(r.withheld.is_empty());
        let fa = d.join("molecules.fasta");
        fs::write(&fa, ">linear.1
AC
>linear.2
GT
>linear.3
TT
").unwrap();
        let mut r = call(&g, &ev, &lin, None, Some(&seeds)).unwrap();
        withhold(&mut r, 750);
        split_withheld(&fa, &d.join("molecules.withheld.fasta"), &r).unwrap();
        assert_eq!(fs::read_to_string(&fa).unwrap(), ">linear.1
AC
>linear.3
TT
");
        assert_eq!(fs::read_to_string(d.join("molecules.withheld.fasta")).unwrap(), ">linear.2
GT
");
        // the organelle graph loses the fragment, its link and its path; the free circle stays
        let og = d.join("organelle.gfa");
        let gfa_text = [
            "H\tVN:Z:1.0\tsn:Z:sample",
            "S\tsample.m\tACGT", "S\tsample.x.1\tACGT", "S\tsample.y\tACGT",
            "L\tsample.m\t+\tsample.m\t+\t0M", "L\tsample.x.1\t+\tsample.x.1\t+\t0M", "L\tsample.y\t+\tsample.y\t+\t0M",
            "P\tsample#1#mt.1\tsample.m+\t*", "P\tsample#2#mt.1\tsample.x.1+\t*", "P\tsample#3#mt.1\tsample.y+\t*",
        ]
        .join("\n");
        fs::write(&og, gfa_text + "\n").unwrap();
        assert_eq!(drop_withheld_gfa(&og, &r).unwrap(), 1);
        let text = fs::read_to_string(&og).unwrap();
        assert!(!text.contains("sample.x"), "{text}");
        assert!(text.contains("S\tsample.y") && text.contains("S\tsample.m") && text.contains("sample#3#mt.1"), "{text}");
        let asm = d.join("assembly.gfa");
        fs::write(&asm, "S\tm\tACGT\nS\tx\tGGCC\nS\ty\tTTAA\n").unwrap();
        write_withheld_segments(&asm, &d.join("components.withheld.fasta"), &r).unwrap();
        assert_eq!(
            fs::read_to_string(d.join("components.withheld.fasta")).unwrap(),
            format!(">x component={} molecule=2 length=4
GGCC
", r.withheld[0].components[0])
        );
        fs::remove_dir_all(&d).ok();
    }
}
