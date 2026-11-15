//! Which configuration of a graph a species panel supports, where the target's reads leave it
//! open (`ovasm pan reconstruct`, `ovasm solve pan-configuration`).
//!
//! `linearize` joins the single-copy segments (anchors) by bridges: walks from one anchor end
//! through repeats to another. Every anchor has two ends and a closed configuration takes one
//! bridge at each, so it is a perfect matching of the anchor ends. Short reads do not span a
//! repeat longer than their fragments; every pairing of its copies is then a bridge the graph
//! allows and the reads rank them alike (an Arabidopsis mitochondrion: 48 anchors, 24 repeats,
//! over 300,000 covers, all at 48 ln 1/2).
//!
//! A panel run (another individual's assembly in the reference graph) holds a bridge when one
//! of its walks holds the first anchor's last `FLANK` own bases followed, at the distance the
//! repeats make, by the last anchor's first `FLANK` own bases ("own": not shared with the
//! neighbour through a link overlap, which every entry of a repeat has alike), in a stretch the
//! run's reads decided (`PanelWalk::cuts`).
//!
//! The model has no constant fitted to results:
//! * Organelles are inherited whole, so the target's junctions are those of its closest
//!   relative unless they changed since. Which run that is follows from the target's alleles at
//!   the graph's variant sites: with a mismatch rate unknown, below one half and integrated
//!   under Jeffreys' prior, d mismatches and m matches have the probability `ln_marginal`
//!   gives (the mean of sin(t)^(2d) cos(t)^(2m) over t in (0, pi/4)); a site the run does not
//!   pass says nothing (1 / alleles). The runs' posterior weights are these, normalized
//!   (`relative_weights`).
//! * At an anchor end with bridges 1..K, `votes` is the weighted share of the runs that hold
//!   each (a run holding several counts for each in equal parts) and `frequency` the share
//!   among all runs alike, with half a count added per bridge (Krichevsky-Trofimov).
//! * A junction that changed since the closest relative is drawn from the species' frequency
//!   (parent-independent mutation): P(bridge) = rho * votes + (1 - rho) * frequency. `rho` is
//!   the panel's own: every run is predicted from the others in turn, and rho is the fixed
//!   point of rho = (sum of the responsibilities rho a / (rho a + (1 - rho) b) + 1/2) / (n + 1)
//!   over the junctions the runs decide (a: the votes for what the run holds, b: its
//!   frequency).
//! * The configuration taken maximizes the expected number of correctly joined anchor ends,
//!   the sum of P over its bridges at both of their ends: a maximum-weight perfect matching.
//!   Expectation is linear, so no independence between junctions is assumed. The matching is
//!   exact: the end graph falls into small components (a repeat's entries and exits) whose
//!   perfect matchings are enumerated.
//!
//! The target's reads come first throughout: a matching is ranked by what `linearize` ranks a
//! cover by from reads (bridges that reads measured before the graph's pairings, shares of
//! reads, copy numbers from depth) and the panel orders only what those leave tied.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;

use crate::evidence::{flip, revcomp, Graph, Oriented};

/// Own bases of an anchor either side of a junction that place a bridge on a panel walk.
pub const FLANK: usize = 31;

/// One walk of a panel run.
pub struct PanelWalk {
    /// PanSN sample of the walk: walks of one run are configurations of one individual.
    pub run: String,
    pub sequence: Vec<u8>,
    pub circular: bool,
    /// Offsets (bases from the walk's start; 0 for a circular walk's closing junction) where the
    /// run's own graph branches: the run's reads did not decide what follows what across them,
    /// so a bridge that spans one is not held. Empty for a walk its reads decided throughout.
    pub cuts: Vec<usize>,
}

/// The target's call at one variant site of the reference graph (`genotype`): per allele the
/// runs whose walks take it and whether it is the target's major allele.
#[derive(Clone, Debug, Default)]
pub struct Site {
    pub alleles: Vec<(BTreeSet<String>, bool)>,
}

/// `genotypes.tsv` as `genotype::run` writes it (columns site, reference_runs, call).
pub fn read_genotypes(path: &Path) -> Result<Vec<Site>> {
    let text = fs::read_to_string(path).with_context(|| format!("cannot read {}", path.display()))?;
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    let col = |name: &str| header.iter().position(|h| *h == name).with_context(|| format!("{}: no column {name}", path.display()));
    let (site, runs, call) = (col("site")?, col("reference_runs")?, col("call")?);
    let mut order: Vec<String> = Vec::new();
    let mut by: BTreeMap<String, Site> = BTreeMap::new();
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() <= site.max(runs).max(call) {
            continue;
        }
        if !by.contains_key(f[site]) {
            order.push(f[site].to_string());
        }
        by.entry(f[site].to_string()).or_default().alleles.push((
            f[runs].split(',').filter(|r| !r.is_empty()).map(str::to_string).collect(),
            f[call] == "major",
        ));
    }
    Ok(order.iter().map(|s| by[s].clone()).collect())
}

/// ln of the probability of `mismatches` and `matches` at sites that each differ at an unknown
/// rate below one half (a relative agrees more often than not; without the bound a run that
/// differs at every site would count as much as one that differs at none), integrated under
/// Jeffreys' prior. With the rate sin^2(t) that prior is uniform in t, so this is the mean of
/// sin(t)^(2 d) cos(t)^(2 m) over t in (0, pi/4): 1 for no site, 1/2 + 1/pi for one match,
/// 1/2 - 1/pi for one mismatch. Midpoint rule, fine enough for the peak's width.
fn ln_marginal(matches: usize, mismatches: usize) -> f64 {
    let n = 4096usize.max((200.0 * ((matches + mismatches) as f64).sqrt()) as usize).min(1 << 20);
    let h = std::f64::consts::FRAC_PI_4 / n as f64;
    let (d, m) = (2.0 * mismatches as f64, 2.0 * matches as f64);
    let terms: Vec<f64> = (0..n)
        .map(|i| {
            let t = (i as f64 + 0.5) * h;
            d * t.sin().ln() + m * t.cos().ln()
        })
        .collect();
    let top = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    top + (terms.iter().map(|v| (v - top).exp()).sum::<f64>() / n as f64).ln()
}

/// Per site and allele, which runs take it; `target`: the allele compared against at each site.
struct Alleles {
    takes: Vec<Vec<Vec<bool>>>,
    /// `ln_marginal` by (matches, mismatches): the panel's leave-one-out asks for the same few.
    marginals: std::cell::RefCell<FxHashMap<(usize, usize), f64>>,
}

impl Alleles {
    fn new(sites: &[Site], runs: &[String]) -> Self {
        let index: FxHashMap<&str, usize> = runs.iter().enumerate().map(|(i, r)| (r.as_str(), i)).collect();
        let takes: Vec<Vec<Vec<bool>>> = sites
            .iter()
            .map(|s| {
                s.alleles
                    .iter()
                    .map(|(holders, _)| {
                        let mut v = vec![false; runs.len()];
                        for h in holders {
                            if let Some(&i) = index.get(h.as_str()) {
                                v[i] = true;
                            }
                        }
                        v
                    })
                    .collect::<Vec<Vec<bool>>>()
            })
            .collect();
        Self { takes, marginals: Default::default() }
    }

    /// The allele a run takes at each site, where it takes exactly one.
    fn of_run(&self, run: usize) -> Vec<Option<usize>> {
        self.takes
            .iter()
            .map(|site| {
                let mut has = site.iter().enumerate().filter(|(_, v)| v[run]).map(|(a, _)| a);
                match (has.next(), has.next()) {
                    (Some(a), None) => Some(a),
                    _ => None,
                }
            })
            .collect()
    }

    /// Posterior weight of each run being the closest relative of an individual with alleles
    /// `target`; `skip` (the individual itself, in the panel's leave-one-out) gets none.
    fn relative_weights(&self, target: &[Option<usize>], skip: Option<usize>, runs: usize) -> Vec<f64> {
        let mut ln = vec![f64::NEG_INFINITY; runs];
        for (r, l) in ln.iter_mut().enumerate() {
            if Some(r) == skip {
                continue;
            }
            let (mut matches, mut mismatches, mut other) = (0usize, 0usize, 0.0f64);
            for (site, t) in self.takes.iter().zip(target) {
                let Some(t) = *t else { continue };
                let n = site.len().max(2) as f64;
                if site[t][r] {
                    matches += 1;
                } else if site.iter().any(|v| v[r]) {
                    mismatches += 1;
                    other -= (n - 1.0).ln();
                } else {
                    other -= n.ln();
                }
            }
            *l = *self.marginals.borrow_mut().entry((matches, mismatches)).or_insert_with(|| ln_marginal(matches, mismatches)) + other;
        }
        let top = ln.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if !top.is_finite() {
            return vec![0.0; runs];
        }
        let mut w: Vec<f64> = ln.iter().map(|l| (l - top).exp()).collect();
        let sum: f64 = w.iter().sum();
        w.iter_mut().for_each(|x| *x /= sum);
        w
    }
}

fn encode(kmer: &[u8]) -> Option<u64> {
    let mut x = 0u64;
    for &b in kmer {
        x = (x << 2)
            | match b {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => return None,
            };
    }
    Some(x)
}

/// A bridge on a sequence: the first anchor's last own bases, the last anchor's first own
/// bases, and the bases between them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Signature {
    a: u64,
    b: u64,
    distance: usize,
}

fn signature(g: &Graph, overlap: &FxHashMap<(Oriented, Oriented), usize>, walk: &[Oriented]) -> Option<Signature> {
    let n = walk.len();
    if n < 2 {
        return None;
    }
    let ov = |a: Oriented, b: Oriented| overlap.get(&(a, b)).copied().unwrap_or(0);
    let (sa, sb) = (g.oriented_seq(walk[0]), g.oriented_seq(walk[n - 1]));
    let mut total = sa.len();
    for w in walk.windows(2) {
        total += g.seqs[w[1].0].len().checked_sub(ov(w[0], w[1]))?;
    }
    let (first, last) = (ov(walk[0], walk[1]), ov(walk[n - 2], walk[n - 1]));
    let a_end = sa.len().checked_sub(first)?;
    let own_b = sb.len().checked_sub(last)?;
    let b_start = total.checked_sub(own_b)?;
    if a_end < FLANK || own_b < FLANK || b_start < a_end {
        return None;
    }
    Some(Signature { a: encode(&sa[a_end - FLANK..a_end])?, b: encode(&sb[last..last + FLANK])?, distance: b_start - a_end })
}

/// Where each probe lies on `seq` (ascending).
fn probe_positions(seq: &[u8], probes: &FxHashSet<u64>) -> FxHashMap<u64, Vec<usize>> {
    let mask: u64 = (1u64 << (2 * FLANK)) - 1;
    let (mut x, mut valid) = (0u64, 0usize);
    let mut out: FxHashMap<u64, Vec<usize>> = FxHashMap::default();
    for (i, &b) in seq.iter().enumerate() {
        let c = match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => {
                valid = 0;
                x = 0;
                continue;
            }
        };
        x = ((x << 2) | c) & mask;
        valid += 1;
        if valid >= FLANK && probes.contains(&x) {
            out.entry(x).or_default().push(i + 1 - FLANK);
        }
    }
    out
}

/// Distance tolerance of a bridge: other individuals' copies of a repeat differ a little in
/// length (as `reference_prior`: within 5%, at least 20 bases).
fn tolerance(distance: usize) -> usize {
    20.max(distance / 20)
}

/// Which bridges one run holds, in any of its walks.
fn observe(walks: &[&PanelWalk], sigs: &[Option<Signature>], probes: &FxHashSet<u64>, longest: usize) -> Vec<bool> {
    let mut held = vec![false; sigs.len()];
    for w in walks {
        let n = w.sequence.len();
        if n < 2 * FLANK {
            continue;
        }
        for reverse in [false, true] {
            let mut s = if reverse { revcomp(&w.sequence) } else { w.sequence.to_ascii_uppercase() };
            let mut cuts: Vec<usize> = w.cuts.iter().map(|&c| if reverse { n - c.min(n) } else { c.min(n) }).collect();
            if w.circular {
                let wrap: Vec<u8> = s.iter().cycle().take(n.min(longest + 2 * FLANK + tolerance(longest) + 1)).copied().collect();
                s.extend(wrap);
                // a cut at the closing junction is at 0 and at n alike
                let again: Vec<usize> = cuts.iter().flat_map(|&c| [c % n, c % n + n, c % n + 2 * n]).collect();
                cuts = again;
            }
            cuts.sort_unstable();
            cuts.dedup();
            let pos = probe_positions(&s, probes);
            for (i, sig) in sigs.iter().enumerate() {
                let Some(sig) = sig else { continue };
                if held[i] {
                    continue;
                }
                let (Some(pa), Some(pb)) = (pos.get(&sig.a), pos.get(&sig.b)) else { continue };
                let tol = tolerance(sig.distance);
                'starts: for &p in pa.iter().take_while(|&&p| p < n) {
                    let want = p + FLANK + sig.distance;
                    let lo = want.saturating_sub(tol);
                    let from = pb.partition_point(|&q| q < lo);
                    for &q in pb[from..].iter().take_while(|&&q| q <= want + tol) {
                        // decided by the run's reads: its graph does not branch inside the junction
                        let next_cut = cuts.partition_point(|&c| c <= p);
                        if cuts.get(next_cut).map_or(true, |&c| c >= q + FLANK) {
                            held[i] = true;
                            break 'starts;
                        }
                    }
                }
            }
        }
    }
    held
}

/// The panel's support of every bridge.
#[derive(Clone, Debug, Serialize)]
pub struct Posterior {
    pub runs: Vec<String>,
    /// Posterior weight of each run being the target's closest relative.
    pub relative_weights: Vec<f64>,
    /// Variant sites of the reference graph the target has a major allele at.
    pub sites: usize,
    /// Probability that a junction is as the closest relative has it, from the panel's own
    /// leave-one-out, and the junctions that estimate rests on.
    pub rho: f64,
    pub rho_junctions: usize,
    /// Per bridge (in the order given): probability that the molecule leaves its first anchor
    /// that way, the relatives' weighted share, the species' frequency, runs holding it and
    /// runs holding any bridge at that anchor end.
    pub probability: Vec<f64>,
    pub votes: Vec<f64>,
    pub frequency: Vec<f64>,
    pub runs_holding: Vec<usize>,
    pub runs_informative: Vec<usize>,
}

/// Votes and frequency of every bridge from the runs' observations; `skip`: a run left out.
fn shares(ends: &[Vec<usize>], held: &[Vec<bool>], weights: &[f64], skip: Option<usize>, bridges: usize) -> (Vec<f64>, Vec<f64>, Vec<usize>) {
    let (mut votes, mut frequency, mut informative) = (vec![0.0; bridges], vec![0.0; bridges], vec![0usize; bridges]);
    for at in ends {
        let (mut count, mut weighted) = (vec![0.0f64; at.len()], vec![0.0f64; at.len()]);
        let (mut n, mut wsum) = (0usize, 0.0f64);
        for (r, h) in held.iter().enumerate() {
            if Some(r) == skip {
                continue;
            }
            let k = at.iter().filter(|&&i| h[i]).count();
            if k == 0 {
                continue;
            }
            n += 1;
            wsum += weights[r];
            for (j, &i) in at.iter().enumerate() {
                if h[i] {
                    count[j] += 1.0 / k as f64;
                    weighted[j] += weights[r] / k as f64;
                }
            }
        }
        for (j, &i) in at.iter().enumerate() {
            frequency[i] = (count[j] + 0.5) / (n as f64 + at.len() as f64 / 2.0);
            votes[i] = if wsum > 0.0 { weighted[j] / wsum } else { frequency[i] };
            informative[i] = n;
        }
    }
    (votes, frequency, informative)
}

/// The panel's posterior for the bridges `walks` (as `linearize` builds them: every bridge in
/// both directions, each leaving its first anchor) of graph `g`.
pub fn posterior(g: &Graph, walks: &[Vec<Oriented>], panel: &[PanelWalk], sites: &[Site]) -> Posterior {
    let mut overlap: FxHashMap<(Oriented, Oriented), usize> = FxHashMap::default();
    for (&(a, b), &o) in g.links.iter().zip(&g.overlaps) {
        overlap.entry((a, b)).or_insert(o);
        overlap.entry((flip(b), flip(a))).or_insert(o);
    }
    let sigs: Vec<Option<Signature>> = walks.iter().map(|w| signature(g, &overlap, w)).collect();
    let probes: FxHashSet<u64> = sigs.iter().flatten().flat_map(|s| [s.a, s.b]).collect();
    let longest = sigs.iter().flatten().map(|s| s.distance).max().unwrap_or(0);
    let runs: Vec<String> = panel.iter().map(|w| w.run.clone()).collect::<BTreeSet<_>>().into_iter().collect();
    let held: Vec<Vec<bool>> = runs
        .iter()
        .map(|r| {
            let own: Vec<&PanelWalk> = panel.iter().filter(|w| &w.run == r).collect();
            observe(&own, &sigs, &probes, longest)
        })
        .collect();
    // the bridges leaving each anchor end
    let mut by_end: BTreeMap<Oriented, Vec<usize>> = BTreeMap::new();
    for (i, w) in walks.iter().enumerate() {
        by_end.entry(w[0]).or_default().push(i);
    }
    let ends: Vec<Vec<usize>> = by_end.into_values().collect();
    let alleles = Alleles::new(sites, &runs);
    let target: Vec<Option<usize>> = sites
        .iter()
        .map(|s| {
            let mut major = s.alleles.iter().enumerate().filter(|(_, a)| a.1).map(|(i, _)| i);
            match (major.next(), major.next()) {
                (Some(a), None) => Some(a),
                _ => None,
            }
        })
        .collect();
    // rho from the panel: each run predicted from the others
    let mut pairs: Vec<(f64, f64)> = Vec::new();
    for j in 0..runs.len() {
        let w = alleles.relative_weights(&alleles.of_run(j), Some(j), runs.len());
        let (votes, frequency, informative) = shares(&ends, &held, &w, Some(j), walks.len());
        for at in &ends {
            let own: Vec<usize> = at.iter().copied().filter(|&i| held[j][i]).collect();
            if at.len() < 2 || own.is_empty() || own.len() == at.len() || informative[at[0]] == 0 {
                continue;
            }
            pairs.push((own.iter().map(|&i| votes[i]).sum(), own.iter().map(|&i| frequency[i]).sum()));
        }
    }
    let mut rho = 0.5f64;
    for _ in 0..1000 {
        let r: f64 = pairs
            .iter()
            .map(|&(a, b)| {
                let (x, y) = (rho * a, (1.0 - rho) * b);
                if x + y > 0.0 { x / (x + y) } else { 0.5 }
            })
            .sum();
        let next = (r + 0.5) / (pairs.len() as f64 + 1.0);
        let moved = (next - rho).abs();
        rho = next;
        if moved < 1e-12 {
            break;
        }
    }
    let weights = alleles.relative_weights(&target, None, runs.len());
    let (votes, frequency, informative) = shares(&ends, &held, &weights, None, walks.len());
    let probability: Vec<f64> = (0..walks.len()).map(|i| rho * votes[i] + (1.0 - rho) * frequency[i]).collect();
    // a bridge without a signature (an anchor shorter than the flank): the end's bridges alike
    let mut probability = probability;
    for at in &ends {
        if at.iter().any(|&i| sigs[i].is_none()) {
            for &i in at {
                probability[i] = 1.0 / at.len() as f64;
            }
        }
    }
    Posterior {
        relative_weights: weights,
        sites: target.iter().flatten().count(),
        rho,
        rho_junctions: pairs.len(),
        probability,
        votes,
        frequency,
        runs_holding: (0..walks.len()).map(|i| held.iter().filter(|h| h[i]).count()).collect(),
        runs_informative: informative,
        runs,
    }
}

/// A bridge between two anchor ends for `best_matching`.
#[derive(Clone, Debug)]
pub struct Edge {
    pub a: usize,
    pub b: usize,
    /// ln of the bridge's share of the reads at its ends (`linearize::exit_fractions`).
    pub read_ln: f64,
    /// A bridge reads did not measure: the graph's pairing through an uncrossed repeat, or weak.
    pub incomplete: bool,
    /// Probability of the bridge at its one end plus that at its other.
    pub prior: f64,
    /// Repeats the bridge lays down.
    pub passes: Vec<usize>,
}

/// The matching `best_matching` takes.
#[derive(Clone, Debug)]
pub struct Matching {
    /// Index of the edge taken at every end.
    pub edge_at: Vec<usize>,
    pub molecules: usize,
    /// Sum of the edges' priors: the expected number of correctly joined ends.
    pub expected: f64,
    /// Matchings the reads and the prior rank alike (saturating).
    pub equally_good: u64,
    pub examined: u64,
    pub uses_incomplete: bool,
    /// The target's single records the matching explains: those through two repeats or more,
    /// and all.
    pub records: (u64, u64),
    /// Components of the end graph given a partial matching (too many open ends or over budget).
    pub partial_components: usize,
}

fn bucket(x: f64) -> i128 {
    (x * 1e9).round() as i128
}

/// Molecules of a matching: anchor k has ends 2k and 2k + 1, joined through the anchor; an end
/// without a partner (`usize::MAX`) ends a linear molecule.
fn molecules(partner: &[usize]) -> usize {
    let mut seen = vec![false; partner.len()];
    let mut n = 0;
    // linear molecules, from each open end to the next
    for start in 0..partner.len() {
        if seen[start] || partner[start] != usize::MAX {
            continue;
        }
        n += 1;
        let mut x = start;
        loop {
            seen[x] = true;
            let y = x ^ 1;
            seen[y] = true;
            if partner[y] == usize::MAX {
                break;
            }
            x = partner[y];
        }
    }
    for start in 0..partner.len() {
        if seen[start] {
            continue;
        }
        n += 1;
        let mut x = start;
        loop {
            seen[x] = true;
            let y = partner[x];
            seen[y] = true;
            x = y ^ 1;
            if x == start {
                break;
            }
        }
    }
    n
}

/// The matching of the `ends` anchor ends (anchor k: ends 2k and 2k + 1) by `edges` that
/// `linearize::cover_order` ranks first on what reads measured (no incomplete bridge if one
/// such matching exists; then the target's single records it explains, `records` of the edges
/// taken: those through two repeats or more, then all, as `ReadFit::order`; then the shares of
/// reads less `copy_penalty` per repeat copy laid down beyond or short of `expected`), then
/// with the greatest prior, then with the fewest molecules. `coupled`: ends one record ties
/// together, so that the edges it needs are weighed in one component and it is counted once.
/// An end no edge reaches is where the graph stops and ends a linear molecule
/// (`edge_at`: `usize::MAX`). The ends that edges reach are all matched where that can be;
/// a component where it cannot (bridges through one repeat and the ends they join) leaves the
/// fewest open, as `cover_order` takes the fewest linear molecules first. A component that
/// would need more than `MAX_OPEN` such ends (a fragmented graph: a shallow or very low-depth
/// sample) or has more matchings than `budget` steps enumerate gets a partial matching of its
/// own (`partial`): the bridges the reads measured or the panel favours, the rest of its ends
/// left open; the other components are matched exactly as before.
#[allow(clippy::too_many_arguments)]
pub fn best_matching(
    ends: usize,
    edges: &[Edge],
    expected: &dyn Fn(usize) -> Option<usize>,
    copy_penalty: f64,
    coupled: &[Vec<usize>],
    records: &dyn Fn(&[usize]) -> (u64, u64),
    budget: u64,
) -> Option<Matching> {
    if ends == 0 || ends % 2 == 1 {
        return None;
    }
    let mut at: Vec<Vec<usize>> = vec![Vec::new(); ends];
    for (i, e) in edges.iter().enumerate() {
        if e.a == e.b || e.a >= ends || e.b >= ends {
            continue;
        }
        at[e.a].push(i);
        at[e.b].push(i);
    }
    if at.iter().all(Vec::is_empty) {
        return None;
    }
    // components of the end graph
    let mut root: Vec<usize> = (0..ends).collect();
    fn find(root: &mut [usize], mut x: usize) -> usize {
        while root[x] != x {
            root[x] = root[root[x]];
            x = root[x];
        }
        x
    }
    // bridges through one repeat belong together: its copies are counted over all of them
    let mut through: FxHashMap<usize, usize> = FxHashMap::default();
    for e in edges {
        if e.a != e.b && e.a < ends && e.b < ends {
            let (x, y) = (find(&mut root, e.a), find(&mut root, e.b));
            root[x.max(y)] = x.min(y);
            for &r in &e.passes {
                let first = *through.entry(r).or_insert(e.a);
                let (x, y) = (find(&mut root, e.a), find(&mut root, first));
                root[x.max(y)] = x.min(y);
            }
        }
    }
    // ends one of the target's records ties together are weighed together
    for group in coupled {
        for pair in group.windows(2) {
            if pair[0] < ends && pair[1] < ends {
                let (x, y) = (find(&mut root, pair[0]), find(&mut root, pair[1]));
                root[x.max(y)] = x.min(y);
            }
        }
    }
    let mut comps: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for x in 0..ends {
        // an end no edge reaches is in no component: the graph stops there
        if at[x].is_empty() {
            continue;
        }
        let r = find(&mut root, x);
        comps.entry(r).or_default().push(x);
    }
    // per component: the best matchings among those without an incomplete edge, and among all
    type Key = (u64, u64, i128, i128);
    struct Best {
        key: Option<Key>,
        ties: Vec<Vec<usize>>, // the edges of each
        count: u64,
    }
    const KEEP: usize = 64;
    /// Every perfect matching of `members`: the first free end takes each of its edges in turn.
    /// False when the budget ran out.
    #[allow(clippy::too_many_arguments)]
    fn enumerate(
        members: &[usize],
        at: &[Vec<usize>],
        edges: &[Edge],
        used: &mut [bool],
        taken: &mut Vec<usize>,
        steps: &mut u64,
        budget: u64,
        visit: &mut dyn FnMut(&[usize]),
    ) -> bool {
        *steps += 1;
        if *steps > budget {
            return false;
        }
        let Some(end) = members.iter().copied().find(|&x| !used[x]) else {
            visit(taken);
            return true;
        };
        for &i in &at[end] {
            let e = &edges[i];
            let other = if e.a == end { e.b } else { e.a };
            if used[other] {
                continue;
            }
            used[end] = true;
            used[other] = true;
            taken.push(i);
            let ok = enumerate(members, at, edges, used, taken, steps, budget, visit);
            taken.pop();
            used[end] = false;
            used[other] = false;
            if !ok {
                return false;
            }
        }
        true
    }
    /// Most ends of one component left open beyond those no edge reaches.
    const MAX_OPEN: usize = 2;
    /// The next choice of `pick.len()` of `n` in lexicographic order; false after the last.
    fn next_choice(pick: &mut [usize], n: usize) -> bool {
        let k = pick.len();
        let mut i = k;
        while i > 0 {
            i -= 1;
            if pick[i] + (k - i) < n {
                pick[i] += 1;
                for j in i + 1..k {
                    pick[j] = pick[j - 1] + 1;
                }
                return true;
            }
        }
        false
    }
    /// A partial matching of one component: its edges that the reads measured or that the
    /// prior favours (`prior` above 1: more than even odds at its two ends on average), measured
    /// first, then the likelier, then the greater read share; each taken when both its ends are
    /// still free. The component's other ends stay open.
    fn partial(members: &[usize], at: &[Vec<usize>], edges: &[Edge]) -> Vec<usize> {
        let mut candidates: Vec<usize> = members.iter().flat_map(|&x| at[x].iter().copied()).collect();
        candidates.sort_unstable();
        candidates.dedup();
        candidates.retain(|&i| !edges[i].incomplete || edges[i].prior > 1.0);
        candidates.sort_by(|&i, &j| {
            let (a, b) = (&edges[i], &edges[j]);
            a.incomplete
                .cmp(&b.incomplete)
                .then(b.prior.total_cmp(&a.prior))
                .then(b.read_ln.total_cmp(&a.read_ln))
                .then(i.cmp(&j))
        });
        let mut used: FxHashMap<usize, ()> = FxHashMap::default();
        let mut taken = Vec::new();
        for i in candidates {
            let e = &edges[i];
            if used.contains_key(&e.a) || used.contains_key(&e.b) {
                continue;
            }
            used.insert(e.a, ());
            used.insert(e.b, ());
            taken.push(i);
        }
        taken
    }
    let mut examined = 0u64;
    let mut partial_components = 0usize;
    let mut options: Vec<[Best; 2]> = Vec::new();
    for members in comps.values() {
        let mut best = [Best { key: None, ties: Vec::new(), count: 0 }, Best { key: None, ties: Vec::new(), count: 0 }];
        // each component its own budget, so that one large component does not withdraw the others
        let (mut steps, mut exhausted) = (0u64, false);
        // Every end joined if that can be; else the fewest left open, two more at a time: each
        // open end ends a linear molecule, and `cover_order` takes the fewest of those first. A
        // graph that lacks one link has such a component (a repeat with one entry too few).
        let mut open = members.len() % 2;
        loop {
            let mut pick: Vec<usize> = (0..open).collect();
            loop {
                let mut used = vec![false; ends];
                for &p in &pick {
                    used[members[p]] = true;
                }
                let mut visit = |taken: &[usize]| {
                    let mut passes: FxHashMap<usize, usize> = FxHashMap::default();
                    let (mut read, mut prior, mut incomplete) = (0.0f64, 0.0f64, false);
                    for &i in taken {
                        let e = &edges[i];
                        read += e.read_ln;
                        prior += e.prior;
                        incomplete |= e.incomplete;
                        for &r in &e.passes {
                            *passes.entry(r).or_default() += 1;
                        }
                    }
                    let deviation: usize = passes.iter().map(|(&r, &got)| expected(r).map_or(0, |want| got.abs_diff(want))).sum();
                    let (joint, threads) = records(taken);
                    let key = (joint, threads, bucket(read - copy_penalty * deviation as f64), bucket(prior));
                    for (class, b) in best.iter_mut().enumerate() {
                        if class == 0 && incomplete {
                            continue;
                        }
                        match b.key {
                            Some(k) if key < k => {}
                            Some(k) if key == k => {
                                b.count = b.count.saturating_add(1);
                                if b.ties.len() < KEEP {
                                    b.ties.push(taken.to_vec());
                                }
                            }
                            _ => {
                                b.key = Some(key);
                                b.count = 1;
                                b.ties = vec![taken.to_vec()];
                            }
                        }
                    }
                };
                if !enumerate(members, &at, edges, &mut used, &mut Vec::new(), &mut steps, budget, &mut visit) {
                    exhausted = true;
                    break;
                }
                if !next_choice(&mut pick, members.len()) {
                    break;
                }
            }
            if exhausted || best[1].key.is_some() {
                break;
            }
            open += 2;
            if open > MAX_OPEN.min(members.len()) {
                break;
            }
        }
        examined = examined.saturating_add(steps);
        if exhausted || best[1].key.is_none() {
            let taken = partial(members, &at, edges);
            partial_components += 1;
            // the partial matching stands for both classes: it is the only one this component has
            best = [
                Best { key: Some((0, 0, 0, 0)), ties: vec![taken.clone()], count: 1 },
                Best { key: Some((0, 0, 0, 0)), ties: vec![taken], count: 1 },
            ];
        }
        options.push(best);
    }
    // reads first, as `cover_order`: a cover of measured bridges alone, where one exists
    let class = if options.iter().all(|b| b[0].key.is_some()) { 0 } else { 1 };
    let lists: Vec<&Vec<Vec<usize>>> = options.iter().map(|b| &b[class].ties).collect();
    let equally_good = options.iter().fold(1u64, |n, b| n.saturating_mul(b[class].count));
    let build = |choice: &[usize]| -> (Vec<usize>, Vec<usize>) {
        let (mut partner, mut edge_at) = (vec![usize::MAX; ends], vec![usize::MAX; ends]);
        for (list, &c) in lists.iter().zip(choice) {
            for &i in &list[c] {
                let e = &edges[i];
                partner[e.a] = e.b;
                partner[e.b] = e.a;
                edge_at[e.a] = i;
                edge_at[e.b] = i;
            }
        }
        (partner, edge_at)
    };
    // among the ties the fewest molecules: every combination when few, else one component at a time
    let mut choice = vec![0usize; lists.len()];
    let combinations = lists.iter().fold(1u64, |n, l| n.saturating_mul(l.len() as u64));
    let mut fewest = molecules(&build(&choice).0);
    if combinations <= 4096 {
        let mut current = choice.clone();
        'all: loop {
            let n = molecules(&build(&current).0);
            if n < fewest {
                fewest = n;
                choice = current.clone();
            }
            let mut c = 0;
            loop {
                if c == current.len() {
                    break 'all;
                }
                current[c] += 1;
                if current[c] < lists[c].len() {
                    break;
                }
                current[c] = 0;
                c += 1;
            }
        }
    } else {
        let mut improved = true;
        while improved {
            improved = false;
            for c in 0..lists.len() {
                for option in 0..lists[c].len() {
                    let mut trial = choice.clone();
                    trial[c] = option;
                    let n = molecules(&build(&trial).0);
                    if n < fewest {
                        fewest = n;
                        choice = trial;
                        improved = true;
                    }
                }
            }
        }
    }
    let (_, edge_at) = build(&choice);
    let mut distinct: Vec<usize> = edge_at.iter().copied().filter(|&i| i != usize::MAX).collect();
    distinct.sort_unstable();
    distinct.dedup();
    let explained = lists.iter().zip(&choice).map(|(list, &c)| records(&list[c])).fold((0, 0), |a, b| (a.0 + b.0, a.1 + b.1));
    Some(Matching {
        records: explained,
        molecules: fewest,
        expected: distinct.iter().map(|&i| edges[i].prior).sum(),
        uses_incomplete: distinct.iter().any(|&i| edges[i].incomplete),
        edge_at,
        equally_good,
        examined,
        partial_components,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_component_with_too_many_open_ends_gets_a_partial_matching() {
        // one end (0) reaches five others through one repeat: at most one of them can be joined
        // to it, so four would stay open. The measured bridge is taken; the others stay open.
        let e = |b: usize, incomplete: bool, prior: f64| Edge { a: 0, b, read_ln: -1.0, incomplete, prior, passes: vec![7] };
        let edges = vec![e(1, true, 0.4), e(2, true, 0.4), e(3, false, 0.4), e(4, true, 0.4), e(5, true, 0.4)];
        let m = best_matching(6, &edges, &|_| None, 0.0, &[], &|_| (0, 0), 1_000).expect("a partial matching");
        assert_eq!(m.partial_components, 1);
        assert_eq!(m.edge_at[0], 2);
        assert_eq!(m.edge_at[3], 2);
        assert!([1, 2, 4, 5].iter().all(|&x| m.edge_at[x] == usize::MAX));
        // nothing measured and the panel favours none: every end stays open
        let edges = vec![e(1, true, 0.4), e(2, true, 0.4), e(3, true, 0.4), e(4, true, 0.4), e(5, true, 0.4)];
        let m = best_matching(6, &edges, &|_| None, 0.0, &[], &|_| (0, 0), 1_000).expect("a partial matching");
        assert!(m.edge_at.iter().all(|&x| x == usize::MAX));
    }

    fn seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize]
            })
            .collect()
    }

    /// Two anchors entering and two leaving one repeat: a, c -> r -> b, d.
    fn repeat_graph(overlap: usize) -> (Graph, Vec<Vec<Oriented>>) {
        let shared_in = seq(overlap, 90);
        let shared_out = seq(overlap, 91);
        let mut r = shared_in.clone();
        r.extend(seq(300, 5));
        r.extend(shared_out.iter());
        let with_tail = |mut s: Vec<u8>, tail: &[u8]| {
            s.extend(tail);
            s
        };
        let with_head = |head: &[u8], s: Vec<u8>| {
            let mut out = head.to_vec();
            out.extend(s);
            out
        };
        let g = Graph {
            names: ["a", "b", "c", "d", "r"].iter().map(|s| s.to_string()).collect(),
            seqs: vec![
                with_tail(seq(400, 1), &shared_in),
                with_head(&shared_out, seq(400, 2)),
                with_tail(seq(400, 3), &shared_in),
                with_head(&shared_out, seq(400, 4)),
                r,
            ],
            links: vec![((0, true), (4, true)), ((2, true), (4, true)), ((4, true), (1, true)), ((4, true), (3, true))],
            overlaps: vec![overlap; 4],
        };
        let fwd = |a: usize, b: usize| vec![(a, true), (4, true), (b, true)];
        let rev = |a: usize, b: usize| vec![(b, false), (4, false), (a, false)];
        let walks = vec![fwd(0, 1), fwd(0, 3), fwd(2, 1), fwd(2, 3), rev(0, 1), rev(0, 3), rev(2, 1), rev(2, 3)];
        (g, walks)
    }

    fn spell(g: &Graph, walk: &[Oriented]) -> Vec<u8> {
        crate::linearize::spell(g, walk, false).unwrap()
    }

    fn run(name: &str, g: &Graph, first: (usize, usize), second: (usize, usize), cuts: Vec<usize>) -> PanelWalk {
        let mut s = spell(g, &[(first.0, true), (4, true), (first.1, true)]);
        s.extend(seq(200, 77));
        s.extend(spell(g, &[(second.0, true), (4, true), (second.1, true)]));
        PanelWalk { run: name.into(), sequence: s, circular: false, cuts }
    }

    #[test]
    fn a_panel_run_holds_the_pairing_it_spells_whatever_the_link_overlap() {
        for overlap in [0, 60] {
            let (g, walks) = repeat_graph(overlap);
            let panel = vec![run("x", &g, (0, 1), (2, 3), Vec::new())];
            let p = posterior(&g, &walks, &panel, &[]);
            // a->b and c->d and their reverses are held; a->d and c->b are not
            assert_eq!(p.runs_holding, vec![1, 0, 0, 1, 1, 0, 0, 1], "overlap {overlap}");
            assert!(p.probability[0] > p.probability[1] && p.probability[3] > p.probability[2]);
            assert!((p.probability[0] + p.probability[1] - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn a_junction_the_runs_own_graph_leaves_open_is_not_held() {
        let (g, walks) = repeat_graph(0);
        // the run's own graph branches either side of the repeat's first copy, not of its second
        let decided = run("x", &g, (0, 1), (2, 3), vec![400, 700]);
        let p = posterior(&g, &walks, &[decided], &[]);
        assert_eq!(p.runs_holding, vec![0, 0, 0, 1, 0, 0, 0, 1]);
        // the same cuts on a circular walk, and a cut at the closing junction
        let (g, walks) = repeat_graph(0);
        let mut circle = run("x", &g, (0, 1), (2, 3), vec![0]);
        circle.circular = true;
        let p = posterior(&g, &walks, &[circle], &[]);
        assert_eq!(p.runs_holding, vec![1, 0, 0, 1, 1, 0, 0, 1]);
    }

    #[test]
    fn the_closest_relative_outweighs_the_majority_and_rho_is_the_panels_own() {
        let (g, walks) = repeat_graph(0);
        // five runs pair a-b, c-d; two pair a-d, c-b; the target's alleles are those of the two
        let mut panel = Vec::new();
        for i in 0..5 {
            panel.push(run(&format!("major{i}"), &g, (0, 1), (2, 3), Vec::new()));
        }
        for i in 0..2 {
            panel.push(run(&format!("minor{i}"), &g, (0, 3), (2, 1), Vec::new()));
        }
        let majors: BTreeSet<String> = (0..5).map(|i| format!("major{i}")).collect();
        let minors: BTreeSet<String> = (0..2).map(|i| format!("minor{i}")).collect();
        let sites: Vec<Site> = (0..40).map(|_| Site { alleles: vec![(majors.clone(), false), (minors.clone(), true)] }).collect();
        let p = posterior(&g, &walks, &panel, &sites);
        assert_eq!(p.sites, 40);
        let minor_weight: f64 = p.runs.iter().zip(&p.relative_weights).filter(|(r, _)| r.starts_with("minor")).map(|(_, w)| w).sum();
        assert!(minor_weight > 0.999, "{minor_weight}");
        // within the panel every run's junctions are those of its closest relative
        assert!(p.rho > 0.9, "{}", p.rho);
        assert!(p.rho_junctions > 0);
        assert!(p.probability[1] > 0.8 && p.probability[1] > p.probability[0], "{:?}", p.probability);
        // without alleles the runs count alike and the majority leads
        let flat = posterior(&g, &walks, &panel, &[]);
        assert!(flat.probability[0] > flat.probability[1]);
        assert!(flat.relative_weights.iter().all(|w| (w - 1.0 / 7.0).abs() < 1e-12));
    }

    #[test]
    fn the_marginal_likelihood_is_jeffreys_over_rates_below_one_half() {
        let pi = std::f64::consts::PI;
        // no site: 1; one site: it matches with 1/2 + 1/pi and differs with 1/2 - 1/pi
        assert!(ln_marginal(0, 0).abs() < 1e-9);
        assert!((ln_marginal(1, 0).exp() - (0.5 + 1.0 / pi)).abs() < 1e-7);
        assert!((ln_marginal(0, 1).exp() - (0.5 - 1.0 / pi)).abs() < 1e-7);
        // one of each: (4 / pi) * integral of sin^2 cos^2 over (0, pi/4) = 1/8
        assert!((ln_marginal(1, 1).exp() - 0.125).abs() < 1e-7);
        // far below one half the bound hardly matters: one more mismatch at 0 of 10 costs
        // B(1/2, 21/2) / B(3/2, 19/2) = 19, a little more with the bound
        let cost = (ln_marginal(10, 0) - ln_marginal(9, 1)).exp();
        assert!(cost > 19.0 && cost < 19.2, "{cost}");
        // a run that differs everywhere is no relative
        assert!(ln_marginal(40, 0) - ln_marginal(0, 40) > 20.0);
        // many sites: the peak is narrow and the rule still holds it (a third mismatch in 343 sites costs 340.5 / 2.5)
        let step = (ln_marginal(341, 2) - ln_marginal(340, 3)).exp();
        assert!((step - 340.5 / 2.5).abs() / (340.5 / 2.5) < 1e-3, "{step}");
    }

    /// Every matching, perfect or not: the most ends joined, then reads, then the prior.
    fn brute_force(ends: usize, edges: &[Edge]) -> (usize, i128, i128) {
        fn go(x: usize, ends: usize, edges: &[Edge], used: &mut Vec<bool>, joined: usize, read: f64, prior: f64, best: &mut (usize, i128, i128)) {
            let Some(e) = (x..ends).find(|&e| !used[e]) else {
                let key = (joined, bucket(read), bucket(prior));
                if key > *best {
                    *best = key;
                }
                return;
            };
            // left open
            used[e] = true;
            go(e + 1, ends, edges, used, joined, read, prior, best);
            for edge in edges.iter().filter(|g| g.a != g.b && (g.a == e || g.b == e)) {
                let o = if edge.a == e { edge.b } else { edge.a };
                if used[o] {
                    continue;
                }
                used[o] = true;
                go(e + 1, ends, edges, used, joined + 2, read + edge.read_ln, prior + edge.prior, best);
                used[o] = false;
            }
            used[e] = false;
        }
        let mut best = (0, i128::MIN, i128::MIN);
        go(0, ends, edges, &mut vec![false; ends], 0, 0.0, 0.0, &mut best);
        best
    }

    #[test]
    fn the_matching_is_the_exhaustive_optimum_on_random_end_graphs() {
        let mut x = 12345u64;
        let mut next = |m: u64| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x % m
        };
        let mut solved = 0;
        for _ in 0..300 {
            let ends = 2 * (2 + next(4) as usize);
            let mut edges = Vec::new();
            for a in 0..ends {
                for b in a + 1..ends {
                    if next(3) == 0 {
                        edges.push(Edge { a, b, read_ln: -(next(3) as f64), incomplete: false, prior: next(1000) as f64 / 1000.0, passes: Vec::new() });
                    }
                }
            }
            // None only where one component would need more than two open ends (or no edge at all)
            let Some(m) = best_matching(ends, &edges, &|_| None, 2.0, &[], &|_| (0, 0), 1_000_000) else { continue };
            solved += 1;
            let mut taken: Vec<usize> = m.edge_at.iter().copied().filter(|&i| i != usize::MAX).collect();
            taken.sort_unstable();
            taken.dedup();
            let joined = m.edge_at.iter().filter(|&&i| i != usize::MAX).count();
            assert_eq!(taken.len() * 2, joined);
            let read: f64 = taken.iter().map(|&i| edges[i].read_ln).sum();
            assert_eq!((joined, bucket(read), bucket(m.expected)), brute_force(ends, &edges));
        }
        assert!(solved > 100, "{solved}");
    }

    #[test]
    fn reads_come_before_the_prior_and_ties_take_the_fewest_molecules() {
        // anchors 0 and 1 (ends 0,1 and 2,3): one circle through both (1-2, 3-0) or two circles (0-1, 2-3)
        let edge = |a, b, read_ln, incomplete, prior| Edge { a, b, read_ln, incomplete, prior, passes: Vec::new() };
        let tied = vec![edge(1, 2, 0.0, true, 1.0), edge(3, 0, 0.0, true, 1.0), edge(0, 1, 0.0, true, 1.0), edge(2, 3, 0.0, true, 1.0)];
        let m = best_matching(4, &tied, &|_| None, 2.0, &[], &|_| (0, 0), 1000).unwrap();
        assert_eq!((m.molecules, m.equally_good), (1, 2));
        // the prior prefers the two circles
        let mut prior = tied.clone();
        prior[2].prior = 1.6;
        prior[3].prior = 1.6;
        let m = best_matching(4, &prior, &|_| None, 2.0, &[], &|_| (0, 0), 1000).unwrap();
        assert_eq!((m.molecules, m.equally_good), (2, 1));
        assert!((m.expected - 3.2).abs() < 1e-12);
        // but bridges that reads measured come first
        let mut reads = prior.clone();
        reads[0].incomplete = false;
        reads[1].incomplete = false;
        let m = best_matching(4, &reads, &|_| None, 2.0, &[], &|_| (0, 0), 1000).unwrap();
        assert_eq!((m.molecules, m.uses_incomplete), (1, false));
        // and so does a copy number: the two circles lay repeat 9 down once where depth says twice
        let mut copies = prior.clone();
        copies[0].passes = vec![9];
        copies[1].passes = vec![9];
        copies[2].passes = vec![9];
        let m = best_matching(4, &copies, &|r| (r == 9).then_some(2), 2.0, &[], &|_| (0, 0), 1000).unwrap();
        assert_eq!(m.molecules, 1);
        // ends no bridge reaches end a linear molecule: 0 - (1-2) - 3
        let open = best_matching(4, &tied[..1], &|_| None, 2.0, &[], &|_| (0, 0), 1000).unwrap();
        assert_eq!((open.molecules, open.edge_at.clone()), (1, vec![usize::MAX, 0, 0, usize::MAX]));
        // three ends that bridges reach cannot all be joined: one stays open, the one the prior misses least
        let odd = vec![edge(1, 2, 0.0, true, 1.0), edge(1, 3, 0.0, true, 1.5)];
        let m = best_matching(4, &odd, &|_| None, 2.0, &[], &|_| (0, 0), 1000).unwrap();
        assert_eq!((m.molecules, m.equally_good, m.edge_at.clone()), (1, 1, vec![usize::MAX, 1, usize::MAX, 1]));
        // no bridge at all: no matching
        assert!(best_matching(4, &[], &|_| None, 2.0, &[], &|_| (0, 0), 1000).is_none());
        // too small a budget: that component falls back to its partial matching (here no bridge is
        // measured or favoured, so every end stays open) instead of withdrawing the matching
        let m = best_matching(4, &tied, &|_| None, 2.0, &[], &|_| (0, 0), 1).unwrap();
        assert_eq!(m.partial_components, 1);
        assert!(m.edge_at.iter().all(|&x| x == usize::MAX));
    }

    #[test]
    fn the_targets_own_records_come_before_the_prior_and_tie_ends_together() {
        // four anchors (ends 0-7), two choices that share nothing: ends 0-3 pair as (0-2, 1-3) or
        // (0-3, 1-2), ends 4-7 as (4-6, 5-7) or (4-7, 5-6); the prior wants the first of each
        let edge = |a, b, prior| Edge { a, b, read_ln: 0.0, incomplete: true, prior, passes: Vec::new() };
        let edges = vec![
            edge(0, 2, 1.0), edge(1, 3, 1.0), edge(0, 3, 0.5), edge(1, 2, 0.5),
            edge(4, 6, 1.0), edge(5, 7, 1.0), edge(4, 7, 0.5), edge(5, 6, 0.5),
        ];
        // one record of the target runs through the bridges 0-3 and 4-7: seven reads
        let across = |taken: &[usize]| if taken.contains(&2) && taken.contains(&6) { (0, 7) } else { (0, 0) };
        // weighed apart, neither half can show it holds the record, and the prior has its way
        let apart = best_matching(8, &edges, &|_| None, 2.0, &[], &across, 1000).unwrap();
        assert_eq!((apart.edge_at[0], apart.edge_at[4], apart.records), (0, 4, (0, 0)));
        // tied together by the record, the bridges it runs through are taken against the prior
        let tied = best_matching(8, &edges, &|_| None, 2.0, &[vec![0, 4]], &across, 1000).unwrap();
        assert_eq!((tied.edge_at[0], tied.edge_at[4], tied.records), (2, 6, (0, 7)));
        assert!((tied.expected - 2.0).abs() < 1e-12 && (apart.expected - 4.0).abs() < 1e-12);
        // a record through two repeats counts before any number of the others
        let joint = |taken: &[usize]| if taken.contains(&0) { (1, 1) } else if taken.contains(&2) { (0, 9) } else { (0, 0) };
        let m = best_matching(8, &edges, &|_| None, 2.0, &[], &joint, 1000).unwrap();
        assert_eq!((m.edge_at[0], m.records), (0, (1, 1)));
        // an end no edge reaches may be named by a record: it stays open and blocks nothing
        let open = best_matching(8, &edges[..5], &|_| None, 2.0, &[vec![0, 5]], &|_| (0, 0), 1000).unwrap();
        assert_eq!((open.edge_at[4], open.edge_at[5], open.edge_at[7]), (4, usize::MAX, usize::MAX));
    }

    #[test]
    fn genotype_tables_group_alleles_by_site() {
        let dir = std::env::temp_dir().join(format!("ovasm-pan-posterior-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let path = dir.join("genotypes.tsv");
        fs::write(&path, "site\tleft\tright\tallele\tallele_bp\treference_runs\tkmers\tdepth\tshare\tcall\n0\t1+\t3+\t2+\t1\tx,y\t2\t90\t1.000\tmajor\n0\t1+\t3+\t9+\t1\tz\t2\t0\t0.000\tabsent\n1\t3+\t5+\t4+\t1\t\t2\t0\t0.000\tno_call\n").unwrap();
        let sites = read_genotypes(&path).unwrap();
        assert_eq!(sites.len(), 2);
        assert_eq!(sites[0].alleles.len(), 2);
        assert!(sites[0].alleles[0].1 && !sites[0].alleles[1].1);
        assert!(sites[0].alleles[0].0.contains("y") && sites[1].alleles[0].0.is_empty());
        fs::remove_dir_all(&dir).unwrap();
    }
}
