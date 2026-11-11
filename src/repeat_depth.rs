//! Local depth evidence for repeat-role experiments on an unchanged assembly graph.
use std::collections::{BTreeMap, BTreeSet};
use std::hash::Hasher;
use std::path::Path;

use anyhow::{bail, Context, Result};
use rustc_hash::FxHasher;
use serde::Serialize;

use crate::evidence::{flip, Graph, Oriented};
use crate::thread_evidence::ThreadCounter;

// The same minimum number of contributing records used by the default linearizer.
const MIN_FLANK_RECORDS: u64 = 2;

#[derive(Clone, Copy, Serialize)]
pub(crate) struct Bounds {
    pub(crate) minimum: f64,
    pub(crate) maximum: f64,
}

#[derive(Clone, Copy, Serialize)]
pub(crate) struct Measure {
    pub(crate) segment_depth: Option<f64>,
    pub(crate) anchor_range: Option<Bounds>,
    pub(crate) reciprocal_comparator: bool,
    pub(crate) repeated_visit_records: u64,
}

#[derive(Serialize)]
pub(crate) struct SegmentCalibration {
    left_anchors: Vec<String>,
    right_anchors: Vec<String>,
    reciprocal_anchors: Vec<String>,
    pub(crate) mean: Measure,
    pub(crate) median: Measure,
}

#[derive(Serialize)]
pub(crate) struct FlankCalibration {
    mode: &'static str,
    graph_fingerprint: String,
    minimum_distinct_records: u64,
    segments: BTreeMap<String, SegmentCalibration>,
}

fn graph_fingerprint(gfa: &str) -> String {
    let mut hash = FxHasher::default();
    hash.write(gfa.as_bytes());
    format!("gfa-bytes-fxhash-v1:{:016x}", hash.finish())
}

pub(crate) fn validate_report_graph(gfa: &Path, report: &serde_json::Value) -> Result<()> {
    let Some(calibration) = report.get("repeat_depth_calibration") else {
        return Ok(());
    };
    let fingerprint = calibration["graph_fingerprint"]
        .as_str()
        .context("repeat depth calibration lacks a graph fingerprint")?;
    let graph = std::fs::read_to_string(gfa)?;
    if fingerprint != graph_fingerprint(&graph) {
        bail!("repeat depth calibration does not match this GFA, including its depth tags");
    }
    Ok(())
}

type Link = (Oriented, Oriented);

fn canonical((a, b): Link) -> Link {
    (a, b).min((flip(b), flip(a)))
}

fn positive(value: Option<f64>) -> Option<f64> {
    value.filter(|x| x.is_finite() && *x > 0.0)
}

fn bounds(
    values: &[Option<f64>],
    left: &BTreeSet<usize>,
    right: &BTreeSet<usize>,
) -> Option<Bounds> {
    // Both physical ends must connect to different anchor segments. Two orientations of
    // the same segment give only one depth measurement.
    if left.is_empty() || right.is_empty() || left.union(right).count() < 2 {
        return None;
    }
    range(values, left.union(right))
}

fn range<'a>(
    values: &[Option<f64>],
    anchors: impl Iterator<Item = &'a usize>,
) -> Option<Bounds> {
    let mut low = f64::INFINITY;
    let mut high: f64 = 0.0;
    for index in anchors {
        let value = values[*index]?;
        low = low.min(value);
        high = high.max(value);
    }
    (high > 0.0).then_some(Bounds { minimum: low, maximum: high })
}

impl FlankCalibration {
    pub(crate) fn new(
        g: &Graph,
        gfa: &str,
        records: &ThreadCounter,
        repeat_fold: f64,
    ) -> (Self, Vec<Measure>, Vec<Measure>) {
        let tags = crate::unify::segment_tags(gfa);
        let (mean, median): (Vec<_>, Vec<_>) = g
            .names
            .iter()
            .zip(&g.seqs)
            .map(|(name, seq)| {
                let Some(segment) = tags.get(name) else {
                    return (None, None);
                };
                let mean = positive(crate::unify::depth(segment, seq.len()).map(|d| d.0));
                let median = positive(segment.get("dm").and_then(|d| d.parse().ok())).or(mean);
                (mean, median)
            })
            .unzip();

        // A link count is bound to oriented endpoints, so parallel overlaps cannot
        // supply a calibrated flank. Count each physical end of every supported link.
        let mut overlap_count: BTreeMap<Link, usize> = BTreeMap::new();
        for &edge in &g.links {
            *overlap_count.entry(canonical(edge)).or_default() += 1;
        }
        let mut degree = vec![[0_usize; 2]; g.names.len()];
        let mut self_link = vec![false; g.names.len()];
        let mut supported = Vec::new();
        for (edge, overlaps) in overlap_count {
            let ((a, a_dir), (b, b_dir)) = edge;
            let support = records.link_records(edge);
            if support == 0 {
                continue;
            }
            let a_side = usize::from(a_dir);
            let b_side = usize::from(!b_dir);
            degree[a][a_side] += overlaps;
            degree[b][b_side] += overlaps;
            if a == b {
                self_link[a] = true;
            } else if overlaps == 1 && support >= MIN_FLANK_RECORDS {
                supported.push((a, a_side, b, b_side));
            }
        }
        let anchor = |s: usize| !self_link[s] && degree[s].iter().all(|&n| n <= 1);
        let mut flanks: Vec<[BTreeSet<usize>; 2]> = vec![Default::default(); g.names.len()];
        for (a, a_side, b, b_side) in supported {
            if anchor(b) {
                flanks[a][a_side].insert(b);
            }
            if anchor(a) {
                flanks[b][b_side].insert(a);
            }
        }

        // A continuous record that revisits the central repeat can support looking at
        // the opposing candidate flank. Compare each flank to the other flank's depth,
        // leaving the tested flank and the repeat itself out of its reference.
        let repeated = records.revisited_segment_records(g.names.len());
        let mut reciprocal = vec![BTreeSet::new(); g.names.len()];
        for (center, [left, right]) in flanks.iter().enumerate() {
            if repeated[center] < MIN_FLANK_RECORDS
                || !median[center].zip(bounds(&median, left, right))
                    .is_some_and(|(d, r)| d >= repeat_fold * r.maximum)
            {
                continue;
            }
            for &a in left.difference(right) {
                reciprocal[a].extend(right.iter().copied());
            }
            for &a in right.difference(left) {
                reciprocal[a].extend(left.iter().copied());
            }
        }

        let mut mean_measurements = Vec::with_capacity(g.names.len());
        let mut median_measurements = Vec::with_capacity(g.names.len());
        let mut segments = BTreeMap::new();
        for (s, [left, right]) in flanks.iter().enumerate() {
            let measure = |v: &[Option<f64>]| {
                let direct = bounds(v, left, right);
                let opposite = range(v, reciprocal[s].iter());
                let anchor_range = match (direct, opposite) {
                    (Some(a), Some(b)) => Some(Bounds {
                        minimum: a.minimum.min(b.minimum),
                        maximum: a.maximum.max(b.maximum),
                    }),
                    (a, b) => a.or(b),
                };
                Measure {
                    segment_depth: v[s],
                    anchor_range,
                    reciprocal_comparator: opposite.is_some(),
                    repeated_visit_records: repeated[s],
                }
            };
            let names = |set: &BTreeSet<usize>| {
                let mut names: Vec<_> = set.iter().map(|&i| g.names[i].clone()).collect();
                names.sort_unstable();
                names
            };
            let mean = measure(&mean);
            let median = measure(&median);
            mean_measurements.push(mean);
            median_measurements.push(median);
            segments.insert(g.names[s].clone(), SegmentCalibration {
                left_anchors: names(left),
                right_anchors: names(right),
                reciprocal_anchors: names(&reciprocal[s]),
                mean,
                median,
            });
        }
        let report = Self {
            mode: "reciprocal_flank_depth_range_v2",
            graph_fingerprint: graph_fingerprint(gfa),
            minimum_distinct_records: MIN_FLANK_RECORDS,
            segments,
        };
        (report, mean_measurements, median_measurements)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(links: Vec<Link>) -> Graph {
        Graph {
            names: ["A", "R", "B", "C"].into_iter().map(String::from).collect(),
            seqs: vec![b"ACGT".to_vec(); 4],
            overlaps: vec![0; links.len()],
            links,
        }
    }

    const LEFT: Link = ((0, true), (1, true));
    const RIGHT: Link = ((1, true), (2, true));
    const ALTERNATE: Link = ((1, true), (3, true));
    const TAGS: &str = "S\tA\tACGT\tdp:f:100\tdm:f:95\n\
                        S\tR\tACGT\tdp:f:200\tdm:f:180\n\
                        S\tB\tACGT\tdp:f:110\tdm:f:105\n\
                        S\tC\tACGT\tdp:f:30\tdm:f:25\n";

    fn count(counter: &mut ThreadCounter, link: Link, id: &str) {
        counter.observe(id.as_bytes(), b"ACGT", &[(link.0, link.1, 0, 4)], false, &[], &[], false)
            .unwrap();
    }

    fn two_records(counter: &mut ThreadCounter, link: Link, prefix: &str) {
        for i in 0..2 {
            count(counter, link, &format!("{prefix}{i}"));
        }
    }

    #[test]
    fn different_depth_regimes_bound_the_ratio_and_duplicate_records_do_not_join() {
        let g = graph(vec![LEFT, RIGHT, ALTERNATE]);
        let mut records = ThreadCounter::default();
        two_records(&mut records, LEFT, "left");
        two_records(&mut records, RIGHT, "right");
        count(&mut records, ALTERNATE, "alternate");
        count(&mut records, ALTERNATE, "alternate");
        let (report, mean, median) = FlankCalibration::new(&g, TAGS, &records, 1.35);
        assert_eq!(report.segments["R"].right_anchors, ["B"]);
        let mean_range = mean[1].anchor_range.unwrap();
        assert_eq!((mean_range.minimum, mean_range.maximum), (100.0, 110.0));
        let median_range = median[1].anchor_range.unwrap();
        assert_eq!((median_range.minimum, median_range.maximum), (95.0, 105.0));

        count(&mut records, ALTERNATE, "second-alternate");
        let (_, mean, median) = FlankCalibration::new(&g, TAGS, &records, 1.35);
        let mean_range = mean[1].anchor_range.unwrap();
        assert_eq!((mean_range.minimum, mean_range.maximum), (30.0, 110.0));
        let median_range = median[1].anchor_range.unwrap();
        assert_eq!((median_range.minimum, median_range.maximum), (25.0, 105.0));
    }

    #[test]
    fn missing_or_parallel_flanks_do_not_supply_a_range() {
        let g = graph(vec![LEFT, RIGHT]);
        let mut records = ThreadCounter::default();
        two_records(&mut records, LEFT, "left");
        two_records(&mut records, RIGHT, "right");
        let partial = TAGS.replace("dp:f:110\tdm:f:105", "dp:f:NaN");
        let (_, mean, median) = FlankCalibration::new(&g, &partial, &records, 1.35);
        assert_eq!(mean[1].segment_depth, Some(200.0));
        assert!(mean[1].anchor_range.is_none());
        assert!(median[1].anchor_range.is_none());

        let mut parallel = graph(vec![LEFT, RIGHT, RIGHT]);
        parallel.overlaps[2] = 1;
        let (_, mean, median) = FlankCalibration::new(&parallel, TAGS, &records, 1.35);
        assert!(mean[1].anchor_range.is_none());
        assert!(median[1].anchor_range.is_none());
    }

    #[test]
    fn two_orientations_of_the_same_flank_are_one_measurement() {
        let back_to_left = ((1, true), (0, true));
        let g = graph(vec![LEFT, back_to_left]);
        let mut records = ThreadCounter::default();
        two_records(&mut records, LEFT, "left");
        two_records(&mut records, back_to_left, "back");
        let (report, _, median) = FlankCalibration::new(&g, TAGS, &records, 1.35);
        assert_eq!(report.segments["R"].left_anchors, ["A"]);
        assert_eq!(report.segments["R"].right_anchors, ["A"]);
        assert!(median[1].anchor_range.is_none());
    }
}


#[cfg(test)]
mod reciprocal_test {
    use super::*;

    #[test]
    fn reciprocal_flanks_require_a_complete_revisited_repeat_and_exclude_self() {
        let left = ((0, true), (1, true));
        let right = ((1, true), (2, true));
        let turn = ((2, true), (1, false));
        let g = Graph {
            names: ["A", "R", "B"].into_iter().map(String::from).collect(),
            seqs: vec![b"ACGT".to_vec(); 3],
            links: vec![left, right, turn],
            overlaps: vec![0; 3],
        };
        let tags = "S\tA\tACGT\tdp:f:100\nS\tR\tACGT\tdp:f:200\nS\tB\tACGT\tdp:f:110\n";
        let mut records = ThreadCounter::default();
        for id in [b"left1", b"left2"] {
            records.observe(id, b"ACGTACGT", &[(left.0, left.1, 0, 4)], false, &[], &[], false)
                .unwrap();
        }
        let walk = vec![(1, true), (2, true), (1, false)];
        let steps = [(right.0, right.1, 0, 4), (turn.0, turn.1, 4, 8)];
        for id in [b"one", b"one"] {
            records.observe(id, b"ACGTACGT", &steps, false, std::slice::from_ref(&walk), &[], false)
                .unwrap();
        }
        let (_, _, low) = FlankCalibration::new(&g, tags, &records, 1.35);
        assert!(low[0].anchor_range.is_none());
        records.observe(b"two", b"ACGTACGT", &steps, false, &[walk], &[], false).unwrap();
        let (report, _, enough) = FlankCalibration::new(&g, tags, &records, 1.35);
        assert_eq!(enough[1].repeated_visit_records, 2);
        assert_eq!(enough[0].anchor_range.unwrap().minimum, 110.0);
        assert_eq!(enough[2].anchor_range.unwrap().minimum, 100.0);
        assert_eq!(report.segments["A"].reciprocal_anchors, ["B"]);
        assert_eq!(report.segments["B"].reciprocal_anchors, ["A"]);
    }
}
