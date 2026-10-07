//! The built-in land-plant seed library (`seeddb/`), compiled into the binary.
//!
//! `ovasm run` uses it for every organelle whose `--seeds` were not given, so a plain
//! `ovasm run --reads ... --organelle both --out ...` works without any reference file.
//! The files are written under `<out>/seeds/` (with their manifest) so a run records exactly
//! which seeds it used.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

/// Version of the built-in library (the `version` field of `seeddb/manifest.json`).
pub const VERSION: &str = "v20261007";

const MITOCHONDRION: &[u8] = include_bytes!("../seeddb/landplants.mitochondrion.fasta");
const PLASTID: &[u8] = include_bytes!("../seeddb/landplants.plastid.fasta");
const MANIFEST: &[u8] = include_bytes!("../seeddb/manifest.json");

fn fasta(organelle: &str) -> Option<(&'static str, &'static [u8])> {
    match organelle {
        "mitochondrion" => Some(("landplants.mitochondrion.fasta", MITOCHONDRION)),
        "plastid" => Some(("landplants.plastid.fasta", PLASTID)),
        _ => None,
    }
}

/// Write the seed FASTA of each requested organelle, and the manifest, into `dir`; return the
/// written files by organelle, in the shape `--seeds` produces.
pub fn write(dir: &Path, organelles: &[&str]) -> Result<BTreeMap<String, Vec<PathBuf>>> {
    fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    fs::write(dir.join("manifest.json"), MANIFEST)
        .with_context(|| format!("writing {}", dir.join("manifest.json").display()))?;
    let mut written: BTreeMap<String, Vec<PathBuf>> = BTreeMap::new();
    for organelle in organelles {
        let (name, bytes) = fasta(organelle)
            .with_context(|| format!("no built-in seeds for {organelle:?}"))?;
        let path = dir.join(name);
        fs::write(&path, bytes).with_context(|| format!("writing {}", path.display()))?;
        written.entry((*organelle).to_string()).or_default().push(path);
    }
    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn records(bytes: &[u8]) -> Vec<String> {
        String::from_utf8_lossy(bytes)
            .lines()
            .filter_map(|l| l.strip_prefix('>').map(str::to_string))
            .collect()
    }

    #[test]
    fn the_library_covers_both_organelles_and_includes_arabidopsis() {
        let mito = records(MITOCHONDRION);
        let plastid = records(PLASTID);
        assert_eq!((mito.len(), plastid.len()), (16, 44));
        assert!(mito.iter().any(|a| a == "NC_037304.1"), "Arabidopsis mitochondrion");
        assert!(plastid.iter().any(|a| a == "NC_000932.1"), "Arabidopsis plastid");
    }

    #[test]
    fn the_manifest_matches_the_compiled_in_files_and_version() {
        let manifest: serde_json::Value = serde_json::from_slice(MANIFEST).unwrap();
        assert_eq!(manifest["version"], VERSION);
        for (organelle, bytes) in [("mitochondrion", MITOCHONDRION), ("plastid", PLASTID)] {
            let db = &manifest["databases"][organelle];
            assert_eq!(db["reference_count"].as_u64().unwrap() as usize, records(bytes).len());
            let bases: usize = String::from_utf8_lossy(bytes)
                .lines()
                .filter(|l| !l.starts_with('>'))
                .map(str::len)
                .sum();
            assert_eq!(db["total_bases"].as_u64().unwrap() as usize, bases);
        }
    }

    #[test]
    fn write_materialises_only_the_requested_organelles() {
        let dir = std::env::temp_dir().join(format!("ovasm-seeddb-test-{}", std::process::id()));
        let written = write(&dir, &["mitochondrion"]).unwrap();
        assert_eq!(written.keys().collect::<Vec<_>>(), ["mitochondrion"]);
        assert!(dir.join("landplants.mitochondrion.fasta").is_file());
        assert!(dir.join("manifest.json").is_file());
        assert!(!dir.join("landplants.plastid.fasta").exists());
        assert!(write(&dir, &["nucleus"]).is_err());
        fs::remove_dir_all(&dir).unwrap();
    }
}
