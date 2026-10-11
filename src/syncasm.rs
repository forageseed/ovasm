//! Native FFI binding to the attributed Syncasm graph kernel distributed with Oatk 1.0.
//! This stage consumes raw reads, performs run-coordinate correction and read-threaded
//! unzipping, and returns a graph restored to base coordinates. No executable or reference
//! assembly is invoked. Organelle selection and sequence export remain explicit stages.

use crate::{evidence::parse_gfa, recruit::peak_rss_mb};
use anyhow::{bail, Context, Result};
use serde::Serialize;
use std::{
    fs,
    path::{Path, PathBuf},
    time::Instant,
};

#[derive(Serialize)]
pub struct SyncasmReport {
    pub graph_engine: &'static str,
    pub upstream: &'static str,
    pub upstream_revision: &'static str,
    pub upstream_patches: &'static [&'static str],
    pub homopolymer_compressed_k: bool,
    pub k: usize,
    pub s: usize,
    pub k_tried: Vec<(usize, f64)>,
    pub min_count: u32,
    pub min_count_auto: bool,
    pub reads: u64,
    pub bases: u64,
    pub unitigs: usize,
    pub circular_unitigs: usize,
    pub links: usize,
    pub total_length: usize,
    pub node_depth: f64,
    pub n50: usize,
    pub longest: usize,
    pub correction_edit_fraction: f64,
    pub edge_fraction: f64,
    pub unzip_rounds: usize,
    pub elapsed_seconds: f64,
    pub peak_rss_mb: Option<f64>,
}

#[cfg(syncasm_native)]
extern "C" {
    fn ov_syncasm(
        reads: *const *const std::ffi::c_char,
        n: i32,
        k: i32,
        s: i32,
        cutoff: i32,
        threads: i32,
        prefix: *const std::ffi::c_char,
        read_count: *mut u64,
        base_count: *mut u64,
    ) -> i32;
}

pub fn run(
    inputs: Vec<PathBuf>,
    k: usize,
    s: usize,
    cutoff: u32,
    threads: usize,
    out: &Path,
    out_json: &Path,
) -> Result<SyncasmReport> {
    if k < 3 || k % 2 == 0 || k > i32::MAX as usize || s < 5 || s > 31 || s >= k {
        bail!("Syncasm requires odd k > s and s in 5..31");
    }
    if cutoff == 0 || cutoff > i32::MAX as u32 || !(1..=256).contains(&threads) {
        bail!("Syncasm requires positive --min-count and threads in 1..256");
    }
    if inputs.is_empty() || inputs.len() > i32::MAX as usize {
        bail!("at least one read input is required");
    }
    for p in &inputs {
        if !p.is_file() {
            bail!("read file not found: {}", p.display());
        }
    }
    #[cfg(not(syncasm_native))]
    {
        let _ = (out, out_json);
        bail!("the published syncmer graph kernel requires a Linux scientific provider; select the server");
    }
    #[cfg(syncasm_native)]
    {
        use std::ffi::CString;
        use std::os::unix::ffi::OsStrExt;
        let started = Instant::now();
        fs::create_dir_all(out.parent().unwrap_or(Path::new(".")))?;
        let prefix = out.with_extension("syncasm");
        let paths: Vec<CString> = inputs
            .iter()
            .map(|p| CString::new(p.as_os_str().as_bytes()))
            .collect::<std::result::Result<_, _>>()?;
        let pointers: Vec<_> = paths.iter().map(|p| p.as_ptr()).collect();
        let c_prefix = CString::new(prefix.as_os_str().as_bytes())?;
        // C owns all temporary graph/read structures for the duration of this call and
        // frees them before returning. C strings and the pointer array remain alive here.
        let (mut read_count, mut bases) = (0u64, 0u64);
        let rc = unsafe {
            ov_syncasm(
                pointers.as_ptr(),
                pointers.len() as i32,
                k as i32,
                s as i32,
                cutoff as i32,
                threads as i32,
                c_prefix.as_ptr(),
                &mut read_count,
                &mut bases,
            )
        };
        if rc != 0 {
            bail!("native Syncasm graph construction failed with status {rc}");
        }
        let final_graph = PathBuf::from(format!("{}.utg.final.gfa", prefix.display()));
        let graph = parse_gfa(&final_graph)
            .context("native Syncasm did not produce a valid base-coordinate graph")?;
        for (i, &(a, b)) in graph.links.iter().enumerate() {
            let overlap = graph.overlaps[i];
            let left = graph.oriented_seq(a);
            let right = graph.oriented_seq(b);
            if left[left.len() - overlap..] != right[..overlap] {
                bail!("native Syncasm link {i} has inconsistent restored-base overlap");
            }
        }
        fs::copy(&final_graph, out)?;
        let mut lengths: Vec<usize> = graph.seqs.iter().map(Vec::len).collect();
        lengths.sort_unstable_by(|a, b| b.cmp(a));
        let total: usize = lengths.iter().sum();
        let mut sum = 0;
        let n50 = lengths
            .iter()
            .find(|&&l| {
                sum += l;
                sum * 2 >= total
            })
            .copied()
            .unwrap_or(0);
        let text = fs::read_to_string(out)?;
        let mut depth = Vec::new();
        for l in text.lines().filter(|l| l.starts_with("S\t")) {
            let f: Vec<_> = l.split('\t').collect();
            if let Some(d) = f
                .iter()
                .find_map(|t| t.strip_prefix("KC:i:").and_then(|v| v.parse::<f64>().ok()))
            {
                depth.push(d / f[2].len() as f64);
            }
        }
        depth.sort_by(f64::total_cmp);
        let report = SyncasmReport {
            graph_engine: "native-syncasm",
            upstream: "Oatk 1.0 / Syncasm (MIT, Chenxi Zhou)",
            upstream_revision: "98cd343ff7e8ff506df7fd5e6e60ceaab7657c69",
            upstream_patches: &[
                "original observation support uses the correction bit, independent of strand",
                "shared HPC graph positions pool distinct original read positions before run-length restoration",
                "original read lengths are saved before HPC instead of rereading inputs",
                "empty LCS appends do not copy null pointers",
                "sequence word comparisons use memcpy for unaligned positions",
            ],
            homopolymer_compressed_k: true,
            k,
            s,
            k_tried: vec![],
            min_count: cutoff,
            min_count_auto: false,
            reads: read_count,
            bases,
            unitigs: graph.seqs.len(),
            circular_unitigs: graph.links.iter().filter(|(a, b)| a == b).count(),
            links: graph.links.len(),
            total_length: total,
            node_depth: depth.get(depth.len() / 2).copied().unwrap_or(0.0),
            n50,
            longest: lengths.first().copied().unwrap_or(0),
            correction_edit_fraction: 0.02,
            edge_fraction: 0.35,
            unzip_rounds: 3,
            elapsed_seconds: started.elapsed().as_secs_f64(),
            peak_rss_mb: peak_rss_mb(),
        };
        fs::create_dir_all(out_json.parent().unwrap_or(Path::new(".")))?;
        fs::write(out_json, serde_json::to_string_pretty(&report)?)?;
        Ok(report)
    }
}

#[cfg(all(test, syncasm_native))]
mod tests {
    use super::*;
    #[test]
    fn shared_graph_positions_have_one_distinct_read_consensus() {
        let dir = std::env::temp_dir().join(format!("ovasm-shared-runs-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let binary = dir.join("shared-runs");
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let compiled = std::process::Command::new("cc")
            .args([
                "-std=gnu11",
                "-ffunction-sections",
                "-fdata-sections",
                "-Wl,--gc-sections",
            ])
            .arg("-I")
            .arg(root.join("vendor/syncasm"))
            .arg(root.join("tests/shared_run_consensus.c"))
            .arg(root.join("vendor/syncasm/syncmer.c"))
            .arg(root.join("vendor/syncasm/kalloc.c"))
            .arg("-lm")
            .arg("-o")
            .arg(&binary)
            .status()
            .unwrap();
        assert!(compiled.success());
        assert!(std::process::Command::new(binary)
            .status()
            .unwrap()
            .success());
    }

    #[test]
    fn corrected_syncmer_support_is_independent_of_strand() {
        let dir = std::env::temp_dir().join(format!("ovasm-syncmer-flags-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let binary = dir.join("correction-flags");
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let compiled = std::process::Command::new("cc")
            .args([
                "-std=gnu11",
                "-ffunction-sections",
                "-fdata-sections",
                "-Wl,--gc-sections",
            ])
            .arg("-I")
            .arg(root.join("vendor/syncasm"))
            .arg(root.join("tests/syncmer_correction_flags.c"))
            .arg("-o")
            .arg(&binary)
            .status()
            .unwrap();
        assert!(compiled.success());
        assert!(std::process::Command::new(binary)
            .status()
            .unwrap()
            .success());
    }

    #[test]
    fn real_read_graph_recovers_eight_independent_circles_in_base_coordinates() {
        let dir =
            std::env::temp_dir().join(format!("ovasm-syncasm-circles-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let mut random = 17u64;
        let mut circles = Vec::new();
        let mut reads = Vec::new();
        for chr in 0..8 {
            let mut seq = Vec::new();
            for _ in 0..6_000 {
                random ^= random << 13;
                random ^= random >> 7;
                random ^= random << 17;
                seq.push(b"ACGT"[(random % 4) as usize]);
            }
            for i in 0..32 {
                let start = i * 173 % seq.len();
                let doubled = seq.repeat(3);
                let mut read = doubled[start..start + seq.len() + 1000].to_vec();
                if chr == 0 && i == 0 {
                    // A long ambiguous-base read exercises correction's word
                    // comparisons; the other reads still determine this circle.
                    read[3500] = b'N';
                }
                if i % 2 == 1 {
                    read = crate::evidence::revcomp(&read);
                }
                reads.extend_from_slice(format!(">c{chr}_r{i}\n").as_bytes());
                reads.extend_from_slice(&read);
                reads.push(b'\n');
            }
            circles.push(seq);
        }
        // The raw counters include this read although its HPC length is below
        // k. It also exercises ambiguous bases and the long-run side array.
        reads.extend_from_slice(b">unused_short_with_Ns\n");
        reads.extend_from_slice(&vec![b'A'; 300]);
        reads.extend_from_slice(b"NNNNCC\n");
        reads.extend_from_slice(b">short_raw_read\nACNTA\n");
        let input = dir.join("reads.fa");
        fs::write(&input, &reads).unwrap();
        let gfa = dir.join("assembly.gfa");
        let report = run(
            vec![input.clone()],
            101,
            31,
            5,
            2,
            &gfa,
            &dir.join("assembly.json"),
        )
        .unwrap();
        assert_eq!(report.unitigs, 8);
        assert_eq!(report.reads, 258);
        assert_eq!(report.bases, 8 * 32 * 7000 + 311);
        let out = dir.join("primary.fa");
        let primary =
            crate::primary::run(&gfa, &out, &dir.join("primary.json"), 10, 1_000_000, 0.9).unwrap();
        assert_eq!(primary.records.len(), 8);
        assert_eq!(primary.total_length, 48_000);
        let mut parser = needletail::parse_fastx_file(out).unwrap();
        let mut matched = vec![false; 8];
        while let Some(r) = parser.next() {
            let r = r.unwrap();
            let seq = r.seq();
            let rc = crate::evidence::revcomp(&seq);
            let found = circles
                .iter()
                .position(|s| {
                    let double = s.repeat(2);
                    double
                        .windows(s.len())
                        .any(|w| w == seq.as_ref() || w == rc.as_slice())
                })
                .unwrap();
            assert!(!matched[found]);
            matched[found] = true;
        }
        assert!(matched.into_iter().all(|x| x));
        assert_eq!(fs::read(input).unwrap(), reads);
    }
}
