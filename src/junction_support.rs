//! Read support of every join a reconstructed molecule makes through a repeat, and of every
//! repeat's copy number along it.
//!
//! A short-read molecule is a walk of the de novo graph. Where the walk leaves a repeat it
//! pairs the unique sequence before the repeat with the unique sequence after it; when the
//! reads do not decide those pairings the linearisation still lays out one walk (rice
//! SAMEA115396721: one 376 kb linear molecule from 18 pairings no read measured). The held-out
//! HiFi reads of that sample cross its two worst joins 0 and 1 times while every base of the
//! molecule has 117-2,361 HiFi reads, so the sequence is the sample's and the joins are not.
//!
//! Each join of unique step `a` (ending at x) through repeat steps to unique step `b` (starting
//! at y) is tested with the target's own reads: a fragment (one read, or a read pair placed in
//! the expected orientation) supports it when it holds a graph-unique k-mer of `a` ending at or
//! before x and one of `b` starting at or after y. How many such fragments to expect is not a
//! constant: the fragments lying wholly within long unique steps give the span distribution of
//! this library, and a join with g repeat bases is expected to be crossed by
//! sum_f max(0, s_f - g - K + 1) / U fragments (s_f anchor span, U positions they were drawn
//! from), scaled by the depth of the shallower flank against the molecule's single-copy depth.
//!   * `supported`: observed fragments are not improbably few for that expectation;
//!   * `unsupported`: so few would occur with probability under `ALPHA` under Poisson, *and*
//!     the count over its expectation lies below the lowest `ALPHA` quantile of the same ratio
//!     at control joins: points every `CONTROL_STEP` bp inside unique steps, where a stretch of
//!     `CONTROL_GAP` bp is certainly joined to what follows, counted the same way. Coverage is
//!     lumpier than Poisson (rice 721: controls fall to 0.43 of expectation at their 0.1%
//!     quantile, so a join with 105 of 142 expected crossings is supported, 0-3 of 97-192 is
//!     not). The reads could cross such a join and do not, so the molecule is cut there;
//!   * `undetermined_long_repeat`: no fragment crosses and the expectation is too low for that
//!     to mean anything (the repeat is longer than the library spans): left joined, reported.
//!
//! Copy numbers: a segment's k-mer depth over the single-copy depth (length-weighted median of
//! once-visited unique steps of at least `MIN_COPY_STEP` bp), with a 95% interval from the
//! spread of those single-copy steps; a walk visiting a segment fewer (or more) times than that
//! interval allows is reported as under- (over-) traversed, and one whose best estimate rounds
//! above (below) the visits while the interval still reaches them as `*_uncertain`; a segment
//! sharing most of its k-mers with the other organelle's reference takes that organelle's
//! reads as well and is reported as such instead. Copy numbers change no sequence.
//!
//! Constants that are not fitted to these data: `K`, `ALPHA`, `CONTROL_GAP`/`CONTROL_STEP`
//! (where controls are placed), `MIN_COPY_STEP`, `SHARED_WITH_OTHER`, `PLACE_TOLERANCE`.

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use rustc_hash::FxHashMap;
use serde::Serialize;

use crate::evidence::revcomp;
use crate::recruit::for_each_rec;

/// Word length of placing reads (as `close`).
pub const K: usize = 31;
/// Probability under which so few crossing fragments make a join unsupported.
pub const ALPHA: f64 = 1e-3;
/// Unique steps at least this long give the single-copy depth and its spread.
const MIN_COPY_STEP: usize = 1000;
/// Placements of one read further apart than this are two different places.
const PLACE_TOLERANCE: i64 = 20;
/// Repeat length of the control joins placed in unique sequence, and their spacing.
const CONTROL_GAP: usize = 100;
const CONTROL_STEP: usize = 250;
/// A segment sharing at least this share of its k-mers with the other organelle's reference
/// takes reads of that organelle too: its depth is not its copy number here.
const SHARED_WITH_OTHER: f64 = 0.5;
/// Search window around the expected start of a step in the molecule (segments of the graph
/// may overlap by up to a k-mer or so where a bubble was collapsed).
const STEP_SEARCH: usize = 3000;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Step {
    pub segment: String,
    pub forward: bool,
    pub start: usize,
    pub end: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct Join {
    pub molecule: String,
    pub left: String,
    pub right: String,
    pub repeat: Vec<String>,
    /// last base of the left unique step, first base of the right one (molecule coordinates;
    /// `y` may exceed the length for the join across the origin of a circle)
    pub x: usize,
    pub y: usize,
    pub repeat_bp: usize,
    pub crossing_fragments: u64,
    pub expected: f64,
    pub p_at_most_observed: f64,
    pub class: &'static str,
}

#[derive(Clone, Debug, Serialize)]
pub struct CopyCall {
    pub segment: String,
    pub length: usize,
    pub depth: f64,
    pub copy_estimate: f64,
    pub copy_ci95: (f64, f64),
    pub visits: usize,
    pub class: &'static str,
}

#[derive(Debug, Default, Serialize)]
pub struct Report {
    pub k: usize,
    pub alpha: f64,
    pub molecules_with_path: usize,
    pub molecules_without_path: usize,
    pub anchors: usize,
    pub placed_reads: u64,
    pub fragments: u64,
    pub library_fragments: usize,
    pub library_positions: f64,
    pub span_median: Option<i64>,
    pub span_max: Option<i64>,
    pub control_points: usize,
    pub control_ratio_floor: f64,
    pub single_copy_depth: Option<f64>,
    pub single_copy_log_sd: Option<f64>,
    pub joins: Vec<Join>,
    pub copies: Vec<CopyCall>,
    pub unsupported: usize,
    pub undetermined_long_repeat: usize,
    pub supported: usize,
    pub under_traversed: usize,
    pub over_traversed: usize,
    pub uncertain_copy_calls: usize,
}

/// One molecule of the output: name, sequence, circular, walk of graph segments.
pub struct Molecule {
    pub name: String,
    pub sequence: Vec<u8>,
    pub circular: bool,
    pub path: Option<Vec<(String, bool)>>,
}

/// Segment sequences and k-mer depths (`dp:f:`) of a GFA.
pub fn read_segments(gfa: &Path) -> Result<BTreeMap<String, (Vec<u8>, f64)>> {
    let text = std::fs::read_to_string(gfa).with_context(|| format!("cannot read {}", gfa.display()))?;
    let mut out = BTreeMap::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.first() == Some(&"S") && f.len() >= 3 {
            let dp = f[3..]
                .iter()
                .find_map(|t| t.strip_prefix("dp:f:").and_then(|v| v.parse::<f64>().ok()))
                .unwrap_or(0.0);
            out.insert(f[1].to_string(), (f[2].as_bytes().to_ascii_uppercase(), dp));
        }
    }
    Ok(out)
}

/// `path=u1+,u26-,...` of a FASTA header.
pub fn header_path(header: &str) -> Option<Vec<(String, bool)>> {
    let p = header.split_whitespace().find_map(|w| w.strip_prefix("path="))?;
    p.split(',')
        .map(|s| {
            let (name, o) = s.split_at(s.len().checked_sub(1)?);
            match o {
                "+" => Some((name.to_string(), true)),
                "-" => Some((name.to_string(), false)),
                _ => None,
            }
        })
        .collect()
}

/// Where each step of `path` lies in `seq`: searched near where the previous step ended,
/// since collapsed bubbles can leave neighbouring segments overlapping. None when a step is
/// not found (the walk does not spell the molecule).
pub fn locate_steps(seq: &[u8], path: &[(String, bool)], segs: &BTreeMap<String, (Vec<u8>, f64)>) -> Option<Vec<Step>> {
    let mut steps = Vec::with_capacity(path.len());
    let mut pos = 0usize;
    for (i, (name, fw)) in path.iter().enumerate() {
        let s = &segs.get(name)?.0;
        let s = if *fw { s.clone() } else { revcomp(s) };
        let probe = &s[..s.len().min(60)];
        let lo = pos.saturating_sub(STEP_SEARCH);
        let hi = (pos + STEP_SEARCH + probe.len()).min(seq.len());
        let found = if lo < hi { find(&seq[lo..hi], probe).map(|j| lo + j) } else { None };
        let start = match found {
            Some(j) => j,
            // the molecule's ends may stop inside their first or last step
            None if i == 0 || i + 1 == path.len() => pos.min(seq.len()),
            None => return None,
        };
        let end = (start + s.len()).min(seq.len());
        steps.push(Step { segment: name.clone(), forward: *fw, start, end });
        pos = end;
    }
    Some(steps)
}

fn find(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || hay.len() < needle.len() {
        return None;
    }
    hay.windows(needle.len()).position(|w| w == needle)
}

fn pack(w: &[u8]) -> Option<u64> {
    let mut v = 0u64;
    for &b in w {
        v = (v << 2)
            | match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return None,
            };
    }
    Some(v)
}

fn canonical(w: &[u8]) -> Option<u64> {
    let f = pack(w)?;
    let r = pack(&revcomp(w))?;
    Some(f.min(r))
}

fn name_of(id: &[u8]) -> Vec<u8> {
    let name = id.split(|&b| b == b' ' || b == b'\t').next().unwrap_or(id);
    name.strip_suffix(b"/1").or_else(|| name.strip_suffix(b"/2")).unwrap_or(name).to_vec()
}

/// One read placed on one molecule: forward?, implied leftmost/rightmost molecule position,
/// first and last anchor position (k-mer starts) agreeing with that placement.
#[derive(Clone, Copy, Debug)]
struct Placed {
    mol: u32,
    fw: bool,
    left: i64,
    right: i64,
    amin: i64,
    amax: i64,
}

/// P(X <= n) for X ~ Poisson(lambda).
pub fn poisson_cdf(n: u64, lambda: f64) -> f64 {
    if lambda <= 0.0 {
        return 1.0;
    }
    let mut term = (-lambda).exp();
    let mut sum = term;
    for i in 1..=n {
        term *= lambda / i as f64;
        sum += term;
    }
    sum.min(1.0)
}

fn weighted_median(mut v: Vec<(f64, f64)>) -> Option<f64> {
    v.retain(|x| x.0.is_finite() && x.1 > 0.0);
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.0.total_cmp(&b.0));
    let total: f64 = v.iter().map(|x| x.1).sum();
    let mut acc = 0.0;
    for (x, w) in &v {
        acc += w;
        if acc * 2.0 >= total {
            return Some(*x);
        }
    }
    v.last().map(|x| x.0)
}

/// Joins, their read support and copy numbers of `molecules` (walks of `segs`) from `reads`.
pub fn analyse(
    molecules: &[Molecule],
    segs: &BTreeMap<String, (Vec<u8>, f64)>,
    reads: &[PathBuf],
    other_organelle: &[Vec<u8>],
) -> Result<(Report, Vec<Vec<Step>>)> {
    let mut report = Report { k: K, alpha: ALPHA, ..Default::default() };
    // steps of every molecule with a walk that spells it
    let mut all_steps: Vec<Vec<Step>> = Vec::new();
    for m in molecules {
        let steps = m.path.as_ref().and_then(|p| locate_steps(&m.sequence, p, segs));
        if steps.is_some() {
            report.molecules_with_path += 1;
        } else {
            report.molecules_without_path += 1;
        }
        all_steps.push(steps.unwrap_or_default());
    }
    let mut visits: BTreeMap<&str, usize> = BTreeMap::new();
    for s in all_steps.iter().flatten() {
        *visits.entry(s.segment.as_str()).or_default() += 1;
    }
    // single-copy depth from long once-visited segments, and its spread
    let singles: Vec<(f64, f64)> = all_steps
        .iter()
        .flatten()
        .filter(|s| visits[s.segment.as_str()] == 1 && s.end - s.start >= MIN_COPY_STEP)
        .map(|s| (segs[&s.segment].1, (s.end - s.start) as f64))
        .collect();
    let base = weighted_median(singles.clone()).filter(|&b| b > 0.0);
    report.single_copy_depth = base;
    let log_sd = base.and_then(|b| {
        let dev: Vec<(f64, f64)> = singles.iter().map(|&(d, w)| ((d / b).ln().abs(), w)).collect();
        weighted_median(dev).map(|mad| (mad * 1.4826).max(0.05))
    });
    report.single_copy_log_sd = log_sd;
    let copy = |seg: &str| -> Option<f64> { base.map(|b| segs[seg].1 / b) };
    // a step is a repeat when the walks visit it more than once or its depth says two copies
    let repeat = |seg: &str| -> bool { visits[seg] > 1 || copy(seg).is_some_and(|c| c.round() >= 2.0) };

    // graph-unique k-mers
    let mut count: FxHashMap<u64, u32> = FxHashMap::default();
    for (seq, _) in segs.values() {
        if seq.len() < K {
            continue;
        }
        for w in seq.windows(K) {
            if let Some(c) = canonical(w) {
                *count.entry(c).or_default() += 1;
            }
        }
    }
    // anchors: graph-unique k-mers inside unique steps, either strand -> (mol, fw, pos)
    let mut anchors: FxHashMap<u64, Option<(u32, bool, u32)>> = FxHashMap::default();
    for (mi, (m, steps)) in molecules.iter().zip(&all_steps).enumerate() {
        for s in steps.iter().filter(|s| !repeat(&s.segment)) {
            if s.end < s.start + K {
                continue;
            }
            for p in s.start..=s.end - K {
                let w = &m.sequence[p..p + K];
                let Some(c) = canonical(w) else { continue };
                if count.get(&c) != Some(&1) {
                    continue;
                }
                for (word, fw) in [(pack(w), true), (pack(&revcomp(w)), false)] {
                    let Some(word) = word else { continue };
                    let v = Some((mi as u32, fw, p as u32));
                    anchors.entry(word).and_modify(|e| if *e != v { *e = None }).or_insert(v);
                }
            }
        }
    }
    report.anchors = anchors.values().filter(|v| v.is_some()).count() / 2;
    drop(count);

    // place reads
    let mut by_name: FxHashMap<Vec<u8>, Vec<Placed>> = FxHashMap::default();
    for path in reads {
        for_each_rec(path, |r| {
            let seq = r.seq.to_ascii_uppercase();
            let n = seq.len() as i64;
            if seq.len() < K {
                return;
            }
            let mut first: Option<Placed> = None;
            for (i, w) in seq.windows(K).enumerate() {
                let Some(word) = pack(w) else { continue };
                let Some(Some((mol, fw, at))) = anchors.get(&word).copied() else { continue };
                let (i, at) = (i as i64, at as i64);
                let left = if fw { at - i } else { at - (n - i - K as i64) };
                match &mut first {
                    None => first = Some(Placed { mol, fw, left, right: left + n, amin: at, amax: at }),
                    Some(p) if p.mol == mol && p.fw == fw && (p.left - left).abs() <= PLACE_TOLERANCE => {
                        p.amin = p.amin.min(at);
                        p.amax = p.amax.max(at);
                    }
                    Some(_) => {}
                }
            }
            if let Some(p) = first {
                by_name.entry(name_of(&r.id)).or_default().push(p);
            }
        })?;
    }
    report.placed_reads = by_name.values().map(|v| v.len() as u64).sum();
    let lens: Vec<i64> = molecules.iter().map(|m| m.sequence.len() as i64).collect();
    // fragments: (mol, amin, amax); a pair placed forward-then-reverse on one molecule is one
    // fragment, across the origin of a circle too; anything else counts read by read. A pair
    // is only a fragment when its span is one the library makes (checked below): mates placed
    // far apart are what a wrong join looks like from the other side, not support for one.
    let mut singles: Vec<(u32, i64, i64)> = Vec::new();
    let mut pairs: Vec<((u32, i64, i64), [(u32, i64, i64); 2])> = Vec::new();
    for v in by_name.values() {
        let mut paired = false;
        if let [a, b] = v[..] {
            let (f, r) = if a.fw { (a, b) } else { (b, a) };
            if f.mol == r.mol && f.fw && !r.fw {
                let l = lens[f.mol as usize];
                let ends = [(f.mol, f.amin, f.amax), (r.mol, r.amin, r.amax)];
                if r.right > f.left {
                    pairs.push(((f.mol, f.amin.min(r.amin), f.amax.max(r.amax)), ends));
                    paired = true;
                } else if molecules[f.mol as usize].circular && r.right + l > f.left {
                    pairs.push(((f.mol, f.amin, r.amax + l), ends));
                    paired = true;
                }
            }
        }
        if !paired {
            for p in v.iter().filter(|_| v.len() <= 2) {
                singles.push((p.mol, p.amin, p.amax));
            }
        }
    }
    // the library's spans: pairs and reads wholly within one unique step of at least
    // `MIN_COPY_STEP` bp (no wrong join lies inside a unique step)
    let unique_steps: Vec<(u32, usize, usize)> = all_steps
        .iter()
        .enumerate()
        .flat_map(|(mi, st)| st.iter().filter(|s| !repeat(&s.segment)).map(move |s| (mi as u32, s.start, s.end)))
        .collect();
    let inside = |&(fm, a, b): &(u32, i64, i64), min_len: usize| -> bool {
        unique_steps
            .iter()
            .any(|&(m, s, e)| m == fm && e - s >= min_len && a >= s as i64 && b + K as i64 <= e as i64)
    };
    let mut local: Vec<i64> = pairs.iter().map(|p| p.0).chain(singles.iter().copied()).filter(|f| inside(f, MIN_COPY_STEP)).map(|f| f.2 - f.1).collect();
    local.sort_unstable();
    let span_limit = local.get(local.len().saturating_sub(1) * 99 / 100).map(|&s| s + s / 10).unwrap_or(i64::MAX);
    let mut fragments = singles;
    for (whole, ends) in pairs {
        if whole.2 - whole.1 <= span_limit {
            fragments.push(whole);
        } else {
            fragments.extend(ends);
        }
    }
    report.fragments = fragments.len() as u64;

    // library: fragments wholly within one long unique step (they can be no longer than the
    // step, so only steps of at least twice the longest such span are used, iterated once)
    let within = |min_len: usize| -> (Vec<i64>, f64) {
        let mut spans = Vec::new();
        let mut positions = 0.0;
        for &(m, s, e) in unique_steps.iter().filter(|u| u.2 - u.1 >= min_len) {
            positions += (e - s) as f64;
            for &(fm, a, b) in &fragments {
                if fm == m && a >= s as i64 && b + K as i64 <= e as i64 {
                    spans.push(b - a);
                }
            }
        }
        (spans, positions)
    };
    let (first_spans, _) = within(MIN_COPY_STEP);
    let longest = first_spans.iter().copied().max().unwrap_or(0).max(K as i64) as usize;
    let (mut spans, positions) = within((2 * longest).max(MIN_COPY_STEP));
    spans.sort_unstable();
    report.library_fragments = spans.len();
    report.library_positions = positions;
    report.span_median = spans.get(spans.len() / 2).copied();
    report.span_max = spans.last().copied();
    let expected_at = |g: usize| -> f64 {
        if positions <= 0.0 {
            return 0.0;
        }
        spans.iter().map(|&s| (s - g as i64 - K as i64 + 1).max(0) as f64).sum::<f64>() / positions
    };
    // Coverage is lumpier than Poisson: the same count at points of unique sequence (a join
    // of CONTROL_GAP bases that is certainly right) gives how far below expectation a true join
    // can fall. A join is unsupported only when its count is below the lowest `ALPHA` share of
    // those control ratios as well as improbable under Poisson.
    let control_expected = expected_at(CONTROL_GAP);
    let mut ratios: Vec<f64> = Vec::new();
    if control_expected > 0.0 {
        let margin = report.span_max.unwrap_or(0).max(0) as usize + K;
        let mut by_mol: FxHashMap<u32, Vec<(i64, i64)>> = FxHashMap::default();
        for &(m, a, b) in &fragments {
            by_mol.entry(m).or_default().push((a, b));
        }
        for &(m, s, e) in &unique_steps {
            let Some(fr) = by_mol.get(&m) else { continue };
            if e < s + 2 * margin + CONTROL_GAP {
                continue;
            }
            let mut p = s + margin;
            while p + CONTROL_GAP + margin <= e {
                let (x, y) = (p as i64, (p + CONTROL_GAP) as i64);
                let n = fr.iter().filter(|&&(a, b)| a + (K as i64) <= x && b >= y).count();
                ratios.push(n as f64 / control_expected);
                p += CONTROL_STEP;
            }
        }
        ratios.sort_by(f64::total_cmp);
    }
    let ratio_floor = if ratios.len() >= 20 {
        ratios[((ratios.len() as f64 * ALPHA).floor() as usize).min(ratios.len() - 1)]
    } else {
        0.0
    };
    report.control_points = ratios.len();
    report.control_ratio_floor = ratio_floor;

    // joins
    for (mi, (m, steps)) in molecules.iter().zip(&all_steps).enumerate() {
        let n = steps.len();
        let uniq: Vec<usize> = (0..n).filter(|&i| !repeat(&steps[i].segment)).collect();
        if uniq.is_empty() {
            continue;
        }
        let pairs: Vec<(usize, usize, bool)> = uniq
            .windows(2)
            .map(|w| (w[0], w[1], false))
            .chain(m.circular.then(|| (*uniq.last().unwrap(), uniq[0], true)))
            .collect();
        let l = m.sequence.len();
        for (ia, ib, wrap) in pairs {
            let inner: Vec<usize> = if wrap { (ia + 1..n).chain(0..ib).collect() } else { (ia + 1..ib).collect() };
            if inner.is_empty() {
                continue; // two unique steps side by side: not a join through a repeat
            }
            let x = steps[ia].end;
            let y = steps[ib].start + if wrap { l } else { 0 };
            let g = y.saturating_sub(x);
            let observed = fragments
                .iter()
                .filter(|&&(fm, a, b)| fm as usize == mi && a + (K as i64) <= x as i64 && b >= y as i64)
                .count() as u64;
            let flank = segs[&steps[ia].segment].1.min(segs[&steps[ib].segment].1);
            let scale = base.map(|b| (flank / b).min(1.0)).unwrap_or(1.0);
            let expected = expected_at(g) * scale;
            let p = poisson_cdf(observed, expected);
            let below_controls = ratios.len() < 20 || (observed as f64) < ratio_floor * expected;
            let class = if p < ALPHA && below_controls {
                "unsupported"
            } else if observed == 0 {
                "undetermined_long_repeat"
            } else {
                "supported"
            };
            report.joins.push(Join {
                molecule: m.name.clone(),
                left: format!("{}{}", steps[ia].segment, if steps[ia].forward { '+' } else { '-' }),
                right: format!("{}{}", steps[ib].segment, if steps[ib].forward { '+' } else { '-' }),
                repeat: inner.iter().map(|&i| format!("{}{}", steps[i].segment, if steps[i].forward { '+' } else { '-' })).collect(),
                x,
                y,
                repeat_bp: g,
                crossing_fragments: observed,
                expected,
                p_at_most_observed: p,
                class,
            });
        }
    }
    report.unsupported = report.joins.iter().filter(|j| j.class == "unsupported").count();
    report.undetermined_long_repeat = report.joins.iter().filter(|j| j.class == "undetermined_long_repeat").count();
    report.supported = report.joins.iter().filter(|j| j.class == "supported").count();

    // copy numbers of every visited segment
    let seen: BTreeSet<&str> = visits.keys().copied().collect();
    let other: rustc_hash::FxHashSet<u64> = other_organelle
        .iter()
        .filter(|s| s.len() >= K)
        .flat_map(|s| s.windows(K).filter_map(canonical).collect::<Vec<_>>())
        .collect();
    let shared = |seg: &str| -> f64 {
        let s = &segs[seg].0;
        if other.is_empty() || s.len() < K {
            return 0.0;
        }
        let words: Vec<u64> = s.windows(K).filter_map(canonical).collect();
        words.iter().filter(|w| other.contains(w)).count() as f64 / words.len().max(1) as f64
    };
    if let (Some(_), Some(sd)) = (base, log_sd) {
        for seg in seen {
            let c = copy(seg).unwrap_or(0.0);
            let ci = (c * (-1.96 * sd).exp(), c * (1.96 * sd).exp());
            let v = visits[seg];
            let off = ci.0 >= v as f64 + 0.5 || ci.1 <= v as f64 - 0.5;
            let class = if off && shared(seg) >= SHARED_WITH_OTHER {
                "shared_with_other_organelle"
            } else if ci.0 >= v as f64 + 0.5 {
                "under_traversed"
            } else if ci.1 <= v as f64 - 0.5 {
                "over_traversed"
            } else if shared(seg) < SHARED_WITH_OTHER && c.round() > v as f64 {
                // the depth's best copy number is more than the visits, its interval reaches them
                "under_traversed_uncertain"
            } else if shared(seg) < SHARED_WITH_OTHER && c.round() < v as f64 {
                "over_traversed_uncertain"
            } else {
                "consistent"
            };
            report.copies.push(CopyCall {
                segment: seg.to_string(),
                length: segs[seg].0.len(),
                depth: segs[seg].1,
                copy_estimate: c,
                copy_ci95: ci,
                visits: v,
                class,
            });
        }
    }
    report.under_traversed = report.copies.iter().filter(|c| c.class == "under_traversed").count();
    report.over_traversed = report.copies.iter().filter(|c| c.class == "over_traversed").count();
    report.uncertain_copy_calls = report.copies.iter().filter(|c| c.class.ends_with("_uncertain")).count();
    Ok((report, all_steps))
}

/// Cut molecules at their unsupported joins (at the start of the unique step after the
/// repeat, so the repeat stays with the sequence before it). A circle with cuts becomes as
/// many linear pieces, starting after its first cut. Returns (name, sequence, circular, note).
pub fn cut(molecules: &[Molecule], joins: &[Join]) -> Vec<(String, Vec<u8>, bool, String)> {
    let mut out = Vec::new();
    for m in molecules {
        let l = m.sequence.len();
        let mut at: Vec<usize> = joins
            .iter()
            .filter(|j| j.molecule == m.name && j.class == "unsupported")
            .map(|j| j.y % l.max(1))
            .collect();
        at.sort_unstable();
        at.dedup();
        if at.is_empty() {
            out.push((m.name.clone(), m.sequence.clone(), m.circular, String::new()));
            continue;
        }
        let mut pieces: Vec<Vec<u8>> = Vec::new();
        if m.circular {
            for (i, &s) in at.iter().enumerate() {
                let e = at[(i + 1) % at.len()];
                let piece = if e > s { m.sequence[s..e].to_vec() } else { [&m.sequence[s..], &m.sequence[..e]].concat() };
                pieces.push(piece);
            }
        } else {
            let mut s = 0;
            for &e in at.iter().chain(std::iter::once(&l)) {
                if e > s {
                    pieces.push(m.sequence[s..e].to_vec());
                }
                s = e;
            }
        }
        let total = pieces.len();
        for (i, p) in pieces.into_iter().enumerate() {
            out.push((
                format!("{}.j{}", m.name, i + 1),
                p,
                false,
                format!("cut_from={} piece={}/{} cut_at=unsupported_join", m.name, i + 1, total),
            ));
        }
    }
    out
}

/// How the panel layer of a reconstruction was made.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum PanelMode {
    /// No panel report: the reads decided, or no panel was given.
    NoPanel,
    /// The panel has runs whose walks it can use.
    Panel,
    /// A panel was given but none of its runs is usable evidence (a graph merged before source
    /// boundaries were recorded leaves every unresolved walk out): nothing is attributed to it.
    /// Merge the panel again with this version.
    NoUsableRuns,
}

pub fn panel_mode(report_present: bool, usable_runs: u64) -> PanelMode {
    match (report_present, usable_runs) {
        (false, _) => PanelMode::NoPanel,
        (true, 0) => PanelMode::NoUsableRuns,
        (true, _) => PanelMode::Panel,
    }
}

/// Where a join's pairing comes from. `bridge`: the panel's (probability, informative, taken)
/// for it; `panel_adopted`: the molecule is the panel's representative.
pub fn join_source(class: &str, bridge: Option<(f64, bool, bool)>, panel_adopted: bool, mode: PanelMode) -> &'static str {
    match class {
        "supported" => "reads",
        "unsupported" => "contradicted_by_reads",
        _ if mode == PanelMode::Panel
            && panel_adopted
            && bridge.is_some_and(|(_, informative, taken)| informative && taken) =>
        {
            "panel"
        }
        _ => "undetermined",
    }
}

#[cfg(test)]
mod source_tests {
    use super::*;

    #[test]
    fn a_join_is_the_panels_only_when_its_adopted_representative_took_it() {
        let b = Some((0.8, true, true));
        let p = PanelMode::Panel;
        assert_eq!(join_source("supported", b, true, p), "reads");
        assert_eq!(join_source("unsupported", b, true, p), "contradicted_by_reads");
        assert_eq!(join_source("undetermined_long_repeat", b, true, p), "panel");
        assert_eq!(join_source("undetermined_long_repeat", b, false, p), "undetermined");
        assert_eq!(join_source("undetermined_long_repeat", Some((0.8, true, false)), true, p), "undetermined");
        assert_eq!(join_source("undetermined_long_repeat", Some((0.5, false, true)), true, p), "undetermined");
        assert_eq!(join_source("undetermined_long_repeat", None, true, p), "undetermined");
        assert_eq!(join_source("undetermined_long_repeat", b, true, PanelMode::NoUsableRuns), "undetermined");
    }

    #[test]
    fn a_panel_without_usable_runs_is_reported() {
        assert_eq!(panel_mode(true, 0), PanelMode::NoUsableRuns);
        assert_eq!(panel_mode(true, 47), PanelMode::Panel);
        assert_eq!(panel_mode(false, 0), PanelMode::NoPanel);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn random(n: usize, seed: &mut u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                *seed ^= *seed << 13;
                *seed ^= *seed >> 7;
                *seed ^= *seed << 17;
                b"ACGT"[(*seed % 4) as usize]
            })
            .collect()
    }

    /// Paired reads (2 x 100 bp, 300 bp fragments every `step` bp) of a circle `genome`.
    fn write_pairs(path: &Path, genome: &[u8], step: usize) {
        let mut f = std::fs::File::create(path).unwrap();
        let l = genome.len();
        let doubled = [genome, genome].concat();
        for (i, s) in (0..l).step_by(step).enumerate() {
            let frag = &doubled[s..s + 300];
            let r1 = &frag[..100];
            let r2 = revcomp(&frag[200..]);
            writeln!(f, "@p{i}/1\n{}\n+\n{}", std::str::from_utf8(r1).unwrap(), "I".repeat(100)).unwrap();
            writeln!(f, "@p{i}/2\n{}\n+\n{}", std::str::from_utf8(&r2).unwrap(), "I".repeat(100)).unwrap();
        }
    }

    fn graph(seed: u64, repeat_len: usize) -> BTreeMap<String, (Vec<u8>, f64)> {
        let mut s = seed;
        let mut g = BTreeMap::new();
        for n in ["A", "B", "C", "D"] {
            g.insert(n.to_string(), (random(3000, &mut s), 50.0));
        }
        g.insert("R".to_string(), (random(repeat_len, &mut s), 100.0));
        g
    }

    fn walk(g: &BTreeMap<String, (Vec<u8>, f64)>, names: &[&str]) -> (Vec<u8>, Vec<(String, bool)>) {
        let seq = names.iter().flat_map(|n| g[*n].0.clone()).collect();
        (seq, names.iter().map(|n| (n.to_string(), true)).collect())
    }

    fn run(repeat_len: usize, tag: &str) -> Report {
        let g = graph(7, repeat_len);
        // truth: A R C B R D (one circle); the molecule claims A R B ... D R C
        let (truth, _) = walk(&g, &["A", "R", "C", "B", "R", "D"]);
        let (seq, path) = walk(&g, &["A", "R", "B", "C", "R", "D"]);
        let dir = std::env::temp_dir().join(format!("ovasm-junction-{tag}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let reads = dir.join("r.fq");
        write_pairs(&reads, &truth, 5);
        let m = Molecule { name: "m".into(), sequence: seq, circular: true, path: Some(path) };
        let (r, _) = analyse(&[m], &g, &[reads], &[]).unwrap();
        std::fs::remove_dir_all(&dir).ok();
        r
    }

    #[test]
    fn a_short_repeat_the_pairs_could_cross_but_do_not_is_unsupported() {
        let r = run(80, "short");
        let ab = r.joins.iter().find(|j| j.left == "A+" && j.right == "B+").unwrap();
        assert_eq!(ab.class, "unsupported", "{ab:?}");
        assert_eq!(ab.crossing_fragments, 0);
        assert!(ab.expected > 10.0);
    }

    #[test]
    fn a_repeat_longer_than_the_fragments_is_undetermined_not_wrong() {
        let r = run(600, "long");
        assert!(r.joins.iter().all(|j| j.class == "undetermined_long_repeat"), "{:?}", r.joins);
        assert_eq!(r.unsupported, 0);
    }

    #[test]
    fn the_true_walk_is_supported_at_every_join() {
        let g = graph(11, 80);
        let (truth, path) = walk(&g, &["A", "R", "C", "B", "R", "D"]);
        let dir = std::env::temp_dir().join(format!("ovasm-junction-true-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let reads = dir.join("r.fq");
        write_pairs(&reads, &truth, 5);
        let m = Molecule { name: "m".into(), sequence: truth, circular: true, path: Some(path) };
        let (r, _) = analyse(&[m], &g, &[reads], &[]).unwrap();
        std::fs::remove_dir_all(&dir).ok();
        assert_eq!(r.joins.len(), 2);
        assert!(r.joins.iter().all(|j| j.class == "supported"), "{:?}", r.joins);
    }

    #[test]
    fn a_segment_at_twice_the_depth_visited_once_is_under_traversed() {
        let mut g = graph(13, 80);
        g.get_mut("C").unwrap().1 = 110.0;
        g.get_mut("D").unwrap().1 = 52.0;
        let (seq, path) = walk(&g, &["A", "R", "C", "B", "R", "D"]);
        let m = Molecule { name: "m".into(), sequence: seq, circular: true, path: Some(path) };
        let (r, _) = analyse(&[m], &g, &[], &[]).unwrap();
        let c = r.copies.iter().find(|c| c.segment == "C").unwrap();
        assert_eq!(c.class, "under_traversed", "{c:?}");
        assert_eq!(r.copies.iter().find(|c| c.segment == "R").unwrap().class, "consistent");
    }

    #[test]
    fn a_deep_segment_shared_with_the_other_organelle_is_not_called_under_traversed() {
        let mut g = graph(17, 80);
        g.get_mut("C").unwrap().1 = 300.0;
        let plastid = vec![g["C"].0.clone()];
        let (seq, path) = walk(&g, &["A", "R", "C", "B", "R", "D"]);
        let m = Molecule { name: "m".into(), sequence: seq, circular: true, path: Some(path) };
        let (r, _) = analyse(&[m], &g, &[], &plastid).unwrap();
        assert_eq!(r.copies.iter().find(|c| c.segment == "C").unwrap().class, "shared_with_other_organelle");
    }

    #[test]
    fn a_join_crossed_by_most_of_its_expected_fragments_is_supported() {
        // the true walk read by a denser library: control joins in unique sequence are counted
        // and the true joins stay within their spread
        let g = graph(19, 80);
        let (truth, path) = walk(&g, &["A", "R", "C", "B", "R", "D"]);
        let dir = std::env::temp_dir().join(format!("ovasm-junction-thin-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let reads = dir.join("r.fq");
        write_pairs(&reads, &truth, 3);
        let m = Molecule { name: "m".into(), sequence: truth, circular: true, path: Some(path) };
        let (r, _) = analyse(&[m], &g, &[reads], &[]).unwrap();
        std::fs::remove_dir_all(&dir).ok();
        assert!(r.control_points >= 20);
        assert!(r.joins.iter().all(|j| j.class == "supported"), "{:?}", r.joins);
    }

    #[test]
    fn cutting_a_circle_at_two_joins_gives_two_linear_pieces_of_all_its_bases() {
        let m = Molecule { name: "m".into(), sequence: b"AAAACCCCGGGGTTTT".to_vec(), circular: true, path: None };
        let j = |y| Join {
            molecule: "m".into(), left: String::new(), right: String::new(), repeat: vec![], x: 0, y,
            repeat_bp: 0, crossing_fragments: 0, expected: 20.0, p_at_most_observed: 0.0, class: "unsupported",
        };
        let out = cut(&[m], &[j(4), j(12)]);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].1, b"CCCCGGGG".to_vec());
        assert_eq!(out[1].1, b"TTTTAAAA".to_vec());
        assert!(out.iter().all(|p| !p.2));
    }

    #[test]
    fn steps_are_found_where_neighbouring_segments_overlap() {
        let mut g = BTreeMap::new();
        g.insert("a".to_string(), (b"ACGTACGTTTGACCA".to_vec(), 1.0));
        g.insert("b".to_string(), (b"GACCATTGGCAACGT".to_vec(), 1.0));
        let seq = b"ACGTACGTTTGACCATTGGCAACGT";
        let steps = locate_steps(seq, &[("a".into(), true), ("b".into(), true)], &g).unwrap();
        assert_eq!((steps[1].start, steps[1].end), (10, 25));
    }

    #[test]
    fn poisson_tail() {
        assert!((poisson_cdf(0, 10.0) - (-10.0f64).exp()).abs() < 1e-12);
        assert!((poisson_cdf(2, 1.0) - 0.9196986).abs() < 1e-6);
    }
}
