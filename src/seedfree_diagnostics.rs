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
            read_batches(std::slice::from_ref(&input)),
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
            read_batches(std::slice::from_ref(&input)),
            &params,
            u64::MAX / params.scale,
        )
        .unwrap();
        let histogram = read_context_histogram(
            read_batches(std::slice::from_ref(&input)),
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

/// `context_supported_keys` with the quartiles over every sampled k-mer of the read,
/// singletons included (the rule up to 0.2.1), for comparison.
fn supported_with_singletons(keys: Vec<u64>, counts: &FxHashMap<u64, u32>) -> Vec<u64> {
    let mut depths: Vec<u32> = keys.iter().map(|key| counts[key]).collect();
    if depths.is_empty() {
        return Vec::new();
    }
    depths.sort_unstable();
    let q = |f: f64| depths[((depths.len() - 1) as f64 * f).round() as usize];
    let (lo, hi) = (q(0.25), q(0.75));
    keys.into_iter()
        .filter(|key| {
            let count = counts[key];
            count >= 2 && count >= lo && count <= hi
        })
        .collect()
}

/// `context_supported_keys` with the quartiles over the read's recurrent k-mers only.
fn supported_recurrent(keys: Vec<u64>, counts: &FxHashMap<u64, u32>) -> Vec<u64> {
    let mut depths: Vec<u32> = keys
        .iter()
        .map(|key| counts[key])
        .filter(|&c| c >= 2)
        .collect();
    if depths.is_empty() {
        return Vec::new();
    }
    depths.sort_unstable();
    let q = |f: f64| depths[((depths.len() - 1) as f64 * f).round() as usize];
    let (lo, hi) = (q(0.25), q(0.75));
    keys.into_iter()
        .filter(|key| {
            let count = counts[key];
            count >= 2 && count >= lo && count <= hi
        })
        .collect()
}

fn mode_of(hist: &BTreeMap<u32, u64>) -> u32 {
    hist.iter()
        .max_by_key(|&(&c, &n)| (n, std::cmp::Reverse(c)))
        .map_or(1, |(&c, _)| c)
}

/// One dataset at a time (all of its files): the sketch's raw count histogram, the read-context
/// histogram with the quartiles over all sampled k-mers (up to 0.2.1) and over recurrent ones,
/// the histogram of per-read median depths, and every read's depth profile as TSV
/// (`<out>/<name>.reads.tsv`: id, length, sampled, recurrent, distinct, p10, q1, median, q3).
/// Manifest: {"name": ["reads.fastq.gz", ...], ...}.
#[test]
#[ignore = "requires OVASM_GATE_AUDIT_MANIFEST and OVASM_GATE_AUDIT_OUT"]
fn audit_depth_gate() {
    let manifest: BTreeMap<String, Vec<PathBuf>> = serde_json::from_str(
        &fs::read_to_string(std::env::var("OVASM_GATE_AUDIT_MANIFEST").unwrap()).unwrap(),
    )
    .unwrap();
    let out = PathBuf::from(std::env::var("OVASM_GATE_AUDIT_OUT").unwrap());
    fs::create_dir_all(&out).unwrap();
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
    let threshold = u64::MAX / params.scale;
    for (name, inputs) in manifest {
        let t0 = Instant::now();
        let sketch = build_sketch(read_batches(&inputs), &params, threshold).unwrap();
        let mut raw = BTreeMap::<u32, u64>::new();
        for &count in sketch.counts.values() {
            *raw.entry(count).or_default() += 1;
        }
        let (mut old, mut new) = (
            rustc_hash::FxHashSet::default(),
            rustc_hash::FxHashSet::default(),
        );
        let mut medians = BTreeMap::<u32, u64>::new();
        let mut tsv = std::io::BufWriter::new(
            fs::File::create(out.join(format!("{name}.reads.tsv"))).unwrap(),
        );
        use std::io::Write;
        writeln!(tsv, "id\tlength\tsampled\trecurrent\tdistinct\tp10\tq1\tmedian\tq3").unwrap();
        let (rx, reader) = read_batches(&inputs);
        for batch in rx {
            let rows: Vec<(Vec<u64>, Vec<u64>, String)> = batch
                .par_iter()
                .map(|r| {
                    let keys: Vec<u64> = sampled(&r.seq, params.k, threshold).collect();
                    let distinct = keys.iter().collect::<rustc_hash::FxHashSet<_>>().len();
                    let mut depths: Vec<u32> = keys
                        .iter()
                        .map(|k| sketch.counts[k])
                        .filter(|&c| c >= 2)
                        .collect();
                    depths.sort_unstable();
                    let q = |f: f64| {
                        if depths.is_empty() {
                            0
                        } else {
                            depths[((depths.len() - 1) as f64 * f).round() as usize]
                        }
                    };
                    let row = format!(
                        "{}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}",
                        String::from_utf8_lossy(&r.id),
                        r.seq.len(),
                        keys.len(),
                        depths.len(),
                        if keys.is_empty() {
                            0.0
                        } else {
                            distinct as f64 / keys.len() as f64
                        },
                        q(0.1),
                        q(0.25),
                        q(0.5),
                        q(0.75)
                    );
                    (
                        supported_with_singletons(keys.clone(), &sketch.counts),
                        supported_recurrent(keys, &sketch.counts),
                        row,
                    )
                })
                .collect();
            for (a, b, row) in rows {
                old.extend(a);
                new.extend(b);
                let median: u32 = row.split('\t').nth(7).unwrap().parse().unwrap();
                if median >= 2 {
                    *medians.entry(median).or_default() += 1;
                }
                writeln!(tsv, "{row}").unwrap();
            }
        }
        reader.join().unwrap().unwrap();
        tsv.flush().unwrap();
        let histogram = |keys: &rustc_hash::FxHashSet<u64>| {
            let mut h = BTreeMap::<u32, u64>::new();
            for key in keys {
                *h.entry(sketch.counts[key]).or_default() += 1;
            }
            h
        };
        let (old_h, new_h) = (histogram(&old), histogram(&new));
        let cut = |h: &BTreeMap<u32, u64>| -> BTreeMap<u32, u64> {
            h.range(..=3000).map(|(&c, &n)| (c, n)).collect()
        };
        let report = json!({
            "name": name,
            "inputs": inputs,
            "bases": sketch.bases,
            "sampled_distinct_kmers": sketch.counts.len(),
            "raw_mode_from_2": nuclear_peak(&sketch.counts),
            "context_peak_with_singletons": mode_of(&old_h),
            "context_peak_recurrent": mode_of(&new_h),
            "read_median_mode": mode_of(&medians),
            "raw_histogram": cut(&raw),
            "context_histogram_with_singletons": cut(&old_h),
            "context_histogram_recurrent": cut(&new_h),
            "read_median_histogram": cut(&medians),
            "seconds": t0.elapsed().as_secs_f64(),
        });
        println!(
            "{name}: raw {} | context (singletons) {} | context (recurrent) {} | read medians {} | {:.0}s",
            report["raw_mode_from_2"],
            report["context_peak_with_singletons"],
            report["context_peak_recurrent"],
            report["read_median_mode"],
            report["seconds"].as_f64().unwrap()
        );
        fs::write(
            out.join(format!("{name}.gate.json")),
            serde_json::to_string_pretty(&report).unwrap(),
        )
        .unwrap();
    }
}
