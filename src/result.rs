//! `<out>/result.json`: the outcome of `ovasm run` for scripts and agents (schema
//! `ovasm.result.v1`). One small file with a stable shape, written whether the run worked or
//! not: the status, per organelle its molecules (name, length, circular) and files (paths
//! relative to `<out>`), the run's time and peak memory, and the warnings. `summary.json` keeps
//! everything else.
//!
//! Status of an organelle: `complete` (molecules, all circular, or one linear molecule whose
//! chromosome call is decisive), `partial` (molecules, some open), `no_molecule` (a graph and no
//! molecule) or `failed` (the run stopped; `error` says why). The run's status is the worst of its
//! organelles'.

use std::fs;
use std::path::Path;

use anyhow::Result;
use serde_json::{json, Value};

const FILES: [&str; 6] =
    ["molecules.fasta", "organelle.gfa", "assembly.gfa", "summary.json", "assessment.json", "evidence.json"];

fn molecules(dir: &Path) -> Vec<Value> {
    let Ok(text) = fs::read_to_string(dir.join("molecules.fasta")) else {
        return Vec::new();
    };
    text.lines()
        .filter_map(|l| l.strip_prefix('>'))
        .map(|h| {
            let mut f = h.split_whitespace();
            let name = f.next().unwrap_or("").to_string();
            let tags: Vec<&str> = f.collect();
            let tag = |k: &str| tags.iter().find_map(|t| t.strip_prefix(k));
            json!({
                "name": name,
                "length": tag("length=").and_then(|x| x.parse::<u64>().ok()),
                "circular": tag("circular=") == Some("true"),
            })
        })
        .collect()
}

fn organelle(out: &Path, dir: &Path, name: &str, summary: Option<&Value>, error: Option<String>) -> Value {
    let mols = molecules(dir);
    let all_circular = !mols.is_empty() && mols.iter().all(|m| m["circular"] == true);
    let decisive = summary.and_then(|s| s["decisive"].as_bool()).unwrap_or(false);
    let status = if error.is_some() {
        "failed"
    } else if all_circular || (mols.len() == 1 && decisive) {
        "complete"
    } else if !mols.is_empty() {
        "partial"
    } else {
        "no_molecule"
    };
    let rel = |f: &str| {
        let p = dir.join(f);
        p.exists().then(|| p.strip_prefix(out).unwrap_or(&p).display().to_string())
    };
    let files: serde_json::Map<String, Value> =
        FILES.iter().filter_map(|f| rel(f).map(|p| (f.to_string(), Value::String(p)))).collect();
    let mut warnings: Vec<Value> = Vec::new();
    if let Some(s) = summary {
        for key in ["warnings", "recruitment_warnings"] {
            if let Some(a) = s[key].as_array() {
                warnings.extend(a.iter().cloned());
            }
        }
        if let Some(a) = s["recruitment"]["warnings"].as_array() {
            warnings.extend(a.iter().cloned());
        }
        if s["k_accepted"] == false {
            warnings.push(json!("no k explained half of the graph with molecules (k_accepted: false)"));
        }
    }
    json!({
        "organelle": name,
        "status": status,
        "error": error,
        "molecules": mols,
        "total_length": mols.iter().filter_map(|m| m["length"].as_u64()).sum::<u64>(),
        "k": summary.map(|s| s["k"].clone()),
        "gene_calls": summary.map(|s| s["assessment"]["gene_calls"].clone()),
        "seconds": summary.map(|s| s["elapsed_seconds"].clone()),
        "files": files,
        "warnings": warnings,
    })
}

/// Write `<out>/result.json` from the pipeline's outcome; returns what was written.
pub fn write(out: &Path, report: &Result<Value>, seconds: f64) -> Result<Value> {
    fs::create_dir_all(out)?;
    let mut organelles = Vec::new();
    let mut peak = crate::recruit::peak_rss_mb().unwrap_or(0.0);
    let mut error = None;
    match report {
        Ok(v) if v["organelles"].is_object() => {
            for (name, s) in v["organelles"].as_object().expect("checked") {
                let e = s["error"].as_str().map(str::to_string);
                organelles.push(organelle(out, &out.join(name), name, Some(s), e));
            }
            peak = v["peak_rss_mb"].as_f64().unwrap_or(peak);
        }
        Ok(v) => {
            let name = v["organelle"].as_str().unwrap_or("organelle");
            organelles.push(organelle(out, out, name, Some(v), None));
            peak = v["peak_rss_mb"].as_f64().unwrap_or(peak);
        }
        Err(e) => error = Some(format!("{e:#}")),
    }
    let rank = |s: &str| match s {
        "complete" => 0,
        "partial" => 1,
        "no_molecule" => 2,
        _ => 3,
    };
    let status = if error.is_some() {
        "failed"
    } else {
        organelles
            .iter()
            .map(|o| o["status"].as_str().unwrap_or("failed"))
            .max_by_key(|s| rank(s))
            .unwrap_or("failed")
    };
    let result = json!({
        "schema": "ovasm.result.v1",
        "version": env!("CARGO_PKG_VERSION"),
        "command": std::env::args().collect::<Vec<_>>(),
        "status": status,
        "error": error,
        "organelles": organelles,
        "seconds": (seconds * 10.0).round() / 10.0,
        "peak_memory_mb": peak.round(),
    });
    fs::write(out.join("result.json"), serde_json::to_string_pretty(&result)? + "\n")?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn result_json_reports_molecules_status_and_failures() {
        let dir = std::env::temp_dir().join(format!("ovasm-result-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join("molecules.fasta"), ">linear.1 circular=true length=16014 path=u0+\nACGT\n").unwrap();
        let ok: Result<Value> = Ok(json!({"organelle": "mitochondrion", "k": 1001, "decisive": false,
            "elapsed_seconds": 1.5, "peak_rss_mb": 420.0, "assessment": {"gene_calls": ["complete"]}}));
        let r = write(&dir, &ok, 2.04).unwrap();
        assert_eq!(r["schema"], "ovasm.result.v1");
        assert_eq!(r["status"], "complete");
        assert_eq!(r["organelles"][0]["molecules"][0]["length"], 16014);
        assert_eq!(r["organelles"][0]["files"]["molecules.fasta"], "molecules.fasta");
        assert_eq!(r["peak_memory_mb"], 420.0);
        let failed: Result<Value> = Err(anyhow::anyhow!("no k gave the reads an assembly graph"));
        let r = write(&dir, &failed, 1.0).unwrap();
        assert_eq!(r["status"], "failed");
        assert!(r["error"].as_str().unwrap().contains("no k"));
        assert!(fs::read_to_string(dir.join("result.json")).unwrap().contains("ovasm.result.v1"));
        fs::remove_dir_all(&dir).ok();
    }
}
