//! Opt-in diagnostics on caller-supplied inputs; never part of default validation.
use super::*;
use serde_json::json;
use std::collections::BTreeMap;
use std::fs;

#[test]
#[ignore = "requires OVASM_PEAK_AUDIT_MANIFEST and OVASM_PEAK_AUDIT_OUT"]
fn audit_nuclear_count_histogram() {
    let manifest = std::env::var("OVASM_PEAK_AUDIT_MANIFEST").unwrap();
    let output = std::env::var("OVASM_PEAK_AUDIT_OUT").unwrap();
    let inputs: Vec<PathBuf> =
        serde_json::from_str(&fs::read_to_string(manifest).unwrap()).unwrap();
    let params = DiscoverParams {
        k: 21,
        scale: 64,
        sketch_bases: 0,
        depth_fold: 3.0,
        min_sampled: 20,
        min_high_frac: 0.9,
        max_spread: 1.5,
        min_distinct: 0.8,
        min_cluster_reads: 20,
        target_depth: None,
        salt: 0,
        dump_profiles: false,
    };
    let mut reports = Vec::new();
    for input in inputs {
        let sketch = build_sketch(
            std::slice::from_ref(&input),
            &params,
            u64::MAX / params.scale,
        )
        .unwrap();
        let mut histogram = BTreeMap::<u32, u64>::new();
        for &count in sketch.counts.values() {
            *histogram.entry(count).or_default() += 1;
        }
        let report = json!({
            "input": input,
            "bases": sketch.bases,
            "sampled_distinct_kmers": sketch.counts.len(),
            "current_nuclear_peak": nuclear_peak(&sketch.counts),
            "histogram": histogram,
        });
        println!("{}: {}", report["input"], report["current_nuclear_peak"]);
        reports.push(report);
        fs::write(&output, serde_json::to_string_pretty(&reports).unwrap()).unwrap();
    }
}

#[test]
#[ignore = "requires OVASM_PEAK_AUDIT_MANIFEST and OVASM_PEAK_AUDIT_OUT"]
fn audit_read_context_nuclear_histogram() {
    let inputs: Vec<PathBuf> = serde_json::from_str(
        &fs::read_to_string(std::env::var("OVASM_PEAK_AUDIT_MANIFEST").unwrap()).unwrap(),
    )
    .unwrap();
    let output = std::env::var("OVASM_PEAK_AUDIT_OUT").unwrap();
    let params = DiscoverParams {
        k: 21,
        scale: 64,
        sketch_bases: 0,
        depth_fold: 3.0,
        min_sampled: 20,
        min_high_frac: 0.9,
        max_spread: 1.5,
        min_distinct: 0.8,
        min_cluster_reads: 20,
        target_depth: None,
        salt: 0,
        dump_profiles: false,
    };
    let mut reports = Vec::new();
    for input in inputs {
        let sketch = build_sketch(
            std::slice::from_ref(&input),
            &params,
            u64::MAX / params.scale,
        )
        .unwrap();
        let histogram = read_context_histogram(
            std::slice::from_ref(&input),
            &sketch,
            &params,
            u64::MAX / params.scale,
        )
        .unwrap();
        let peak = histogram
            .iter()
            .max_by_key(|&(&c, &n)| (n, std::cmp::Reverse(c)))
            .map_or(1, |(&c, _)| c);
        let report = json!({"input": input, "raw_nuclear_peak": nuclear_peak(&sketch.counts),
            "read_context_nuclear_peak": peak, "read_context_histogram": histogram});
        println!(
            "{}: {} -> {}",
            report["input"], report["raw_nuclear_peak"], peak
        );
        reports.push(report);
        fs::write(&output, serde_json::to_string_pretty(&reports).unwrap()).unwrap();
    }
}
