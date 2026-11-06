//! The built-in seed libraries (`seeddb/`), compiled into the binary: land plants (mitochondrion
//! and plastid), animals and fungi (mitochondrion; neither has a plastid). `--lineage` picks one.
//!
//! `ovasm run` uses the chosen library for every organelle whose `--seeds` were not given, so a
//! plain `ovasm run hifi --reads ... --organelle both --out ...` works without any reference file.
//! The files are written under `<out>/seeds/` (with their manifest) so a run records exactly
//! which seeds it used.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

/// Version of the built-in land-plant library (the `version` field of `seeddb/manifest.json`).
pub const VERSION: &str = "v20261007";
/// Version of the animal and fungal libraries (`seeddb/{metazoa,fungi}.manifest.json`).
pub const LINEAGE_VERSION: &str = "v20261104";

const MITOCHONDRION: &[u8] = include_bytes!("../seeddb/landplants.mitochondrion.fasta");
const PLASTID: &[u8] = include_bytes!("../seeddb/landplants.plastid.fasta");
const MANIFEST: &[u8] = include_bytes!("../seeddb/manifest.json");
const METAZOA_MITOCHONDRION: &[u8] = include_bytes!("../seeddb/metazoa.mitochondrion.fasta");
const METAZOA_MANIFEST: &[u8] = include_bytes!("../seeddb/metazoa.manifest.json");
const FUNGI_MITOCHONDRION: &[u8] = include_bytes!("../seeddb/fungi.mitochondrion.fasta");
const FUNGI_MANIFEST: &[u8] = include_bytes!("../seeddb/fungi.manifest.json");

/// Which built-in seed library seeds the organelles without `--seeds`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum Lineage {
    /// Land plants: 16 mitochondria and 44 plastids.
    #[default]
    Landplants,
    /// Animals: 32 RefSeq mitochondria, sponge to human.
    Metazoa,
    /// Fungi: 19 RefSeq mitochondria, Blastocladiomycota to Basidiomycota.
    Fungi,
}

impl Lineage {
    pub fn name(self) -> &'static str {
        match self {
            Lineage::Landplants => "landplants",
            Lineage::Metazoa => "metazoa",
            Lineage::Fungi => "fungi",
        }
    }

    pub fn version(self) -> &'static str {
        match self {
            Lineage::Landplants => VERSION,
            Lineage::Metazoa | Lineage::Fungi => LINEAGE_VERSION,
        }
    }

    /// Whether the library holds plastid seeds (animals and fungi have no plastid).
    pub fn has_plastid(self) -> bool {
        self.fasta("plastid").is_some()
    }

    fn manifest(self) -> &'static [u8] {
        match self {
            Lineage::Landplants => MANIFEST,
            Lineage::Metazoa => METAZOA_MANIFEST,
            Lineage::Fungi => FUNGI_MANIFEST,
        }
    }

    pub(crate) fn fasta(self, organelle: &str) -> Option<(&'static str, &'static [u8])> {
        match (self, organelle) {
            (Lineage::Landplants, "mitochondrion") => {
                Some(("landplants.mitochondrion.fasta", MITOCHONDRION))
            }
            (Lineage::Landplants, "plastid") => Some(("landplants.plastid.fasta", PLASTID)),
            (Lineage::Metazoa, "mitochondrion") => {
                Some(("metazoa.mitochondrion.fasta", METAZOA_MITOCHONDRION))
            }
            (Lineage::Fungi, "mitochondrion") => {
                Some(("fungi.mitochondrion.fasta", FUNGI_MITOCHONDRION))
            }
            _ => None,
        }
    }
}

/// Give every organelle that needs seeds and has none the built-in library of `lineage`. The
/// organelles that need seeds are the one asked for (`single`; both when `None`) and, unless
/// `single_seed`, the other one too, so that each organelle claims its own reads. Seeds already in
/// `seeds` are left alone. An organelle asked for that the library lacks (a plastid of an animal or
/// a fungus) is an error; the other organelle the library lacks is left unseeded. Returns the
/// organelles that took the library (empty: nothing was written).
pub fn fill_missing(
    seeds: &mut BTreeMap<String, Vec<PathBuf>>,
    dir: &Path,
    single: Option<&str>,
    single_seed: bool,
    lineage: Lineage,
) -> Result<BTreeSet<String>> {
    let asked: Vec<&str> = match single {
        Some(organelle) => vec![organelle],
        None => vec!["mitochondrion", "plastid"],
    };
    let mut needed = asked.clone();
    if let (Some(organelle), false) = (single, single_seed) {
        needed.push(if organelle == "mitochondrion" {
            "plastid"
        } else {
            "mitochondrion"
        });
    }
    let mut missing = Vec::new();
    for organelle in needed.into_iter().filter(|o| !seeds.contains_key(*o)) {
        if lineage.fasta(organelle).is_some() {
            missing.push(organelle);
        } else if asked.contains(&organelle) {
            bail!(
                "the built-in {} seed library has no {organelle} seeds: give --seeds \
                 {organelle}=PATH, or --organelle mitochondrion",
                lineage.name()
            );
        }
    }
    if missing.is_empty() {
        return Ok(BTreeSet::new());
    }
    seeds.extend(write(dir, &missing, lineage)?);
    Ok(missing.into_iter().map(str::to_string).collect())
}

/// Write the seed FASTA of each requested organelle from `lineage`'s library, and its manifest (as
/// `manifest.json`), into `dir`; return the written files by organelle, in the shape `--seeds`
/// produces.
pub fn write(
    dir: &Path,
    organelles: &[&str],
    lineage: Lineage,
) -> Result<BTreeMap<String, Vec<PathBuf>>> {
    fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    fs::write(dir.join("manifest.json"), lineage.manifest())
        .with_context(|| format!("writing {}", dir.join("manifest.json").display()))?;
    let mut written: BTreeMap<String, Vec<PathBuf>> = BTreeMap::new();
    for organelle in organelles {
        let (name, bytes) = lineage.fasta(organelle).with_context(|| {
            format!("no built-in {} seeds for {organelle:?}", lineage.name())
        })?;
        let path = dir.join(name);
        fs::write(&path, bytes).with_context(|| format!("writing {}", path.display()))?;
        written
            .entry((*organelle).to_string())
            .or_default()
            .push(path);
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
        assert!(
            mito.iter().any(|a| a == "NC_037304.1"),
            "Arabidopsis mitochondrion"
        );
        assert!(
            plastid.iter().any(|a| a == "NC_000932.1"),
            "Arabidopsis plastid"
        );
    }

    #[test]
    fn the_animal_and_fungal_libraries_hold_their_references() {
        let animals = records(METAZOA_MITOCHONDRION);
        let fungi = records(FUNGI_MITOCHONDRION);
        assert_eq!((animals.len(), fungi.len()), (32, 19));
        for (lib, acc, what) in [
            (&animals, "NC_012920.1", "human"),
            (&animals, "NC_024511.2", "Drosophila melanogaster"),
            (&animals, "NC_008944.1", "the sponge Amphimedon"),
            (&fungi, "NC_001224.1", "Saccharomyces cerevisiae"),
            (&fungi, "NC_001715.1", "Allomyces macrogynus"),
        ] {
            assert!(lib.iter().any(|a| a == acc), "{what}");
        }
    }

    #[test]
    fn every_manifest_matches_its_compiled_in_files_and_version() {
        for lineage in [Lineage::Landplants, Lineage::Metazoa, Lineage::Fungi] {
            let manifest: serde_json::Value = serde_json::from_slice(lineage.manifest()).unwrap();
            assert_eq!(manifest["version"], lineage.version());
            for organelle in ["mitochondrion", "plastid"] {
                let Some((name, bytes)) = lineage.fasta(organelle) else {
                    assert!(manifest["databases"][organelle].is_null());
                    continue;
                };
                let db = &manifest["databases"][organelle];
                assert_eq!(db["file"], name);
                assert_eq!(
                    db["reference_count"].as_u64().unwrap() as usize,
                    records(bytes).len()
                );
                let bases: usize = String::from_utf8_lossy(bytes)
                    .lines()
                    .filter(|l| !l.starts_with('>'))
                    .map(str::len)
                    .sum();
                assert_eq!(db["total_bases"].as_u64().unwrap() as usize, bases);
            }
        }
    }

    fn tmp(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!("ovasm-seeddb-{name}-{}", std::process::id()))
    }

    fn names(set: &BTreeSet<String>) -> Vec<&str> {
        set.iter().map(String::as_str).collect()
    }

    #[test]
    fn nothing_given_means_both_organelles_take_the_library() {
        let dir = tmp("fill-both");
        let mut seeds = BTreeMap::new();
        let filled = fill_missing(&mut seeds, &dir, None, false, Lineage::Landplants).unwrap();
        assert_eq!(names(&filled), ["mitochondrion", "plastid"]);
        assert_eq!(seeds.len(), 2);
        assert!(dir.join("landplants.mitochondrion.fasta").is_file());
        assert!(dir.join("landplants.plastid.fasta").is_file());
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn an_animal_mitochondrion_takes_the_animal_library_and_no_plastid() {
        let dir = tmp("fill-metazoa");
        let mut seeds = BTreeMap::new();
        let filled =
            fill_missing(&mut seeds, &dir, Some("mitochondrion"), false, Lineage::Metazoa)
                .unwrap();
        assert_eq!(names(&filled), ["mitochondrion"]);
        assert!(dir.join("metazoa.mitochondrion.fasta").is_file());
        assert!(!dir.join("landplants.plastid.fasta").exists());
        let manifest: serde_json::Value =
            serde_json::from_slice(&fs::read(dir.join("manifest.json")).unwrap()).unwrap();
        assert_eq!(manifest["lineage"], "metazoa");
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn a_plastid_asked_of_the_fungal_library_is_an_error() {
        let dir = tmp("fill-fungi-pt");
        let mut seeds = BTreeMap::new();
        let err = fill_missing(&mut seeds, &dir, None, false, Lineage::Fungi).unwrap_err();
        assert!(err.to_string().contains("no plastid seeds"), "{err}");
        let mut seeds = BTreeMap::new();
        assert!(fill_missing(&mut seeds, &dir, Some("plastid"), true, Lineage::Fungi).is_err());
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn seeds_the_user_gave_are_kept_and_not_reported_as_the_library() {
        let dir = tmp("fill-mixed");
        let mine = PathBuf::from("my_mt.fasta");
        let mut seeds = BTreeMap::from([("mitochondrion".to_string(), vec![mine.clone()])]);
        let filled = fill_missing(&mut seeds, &dir, None, false, Lineage::Landplants).unwrap();
        assert_eq!(names(&filled), ["plastid"]);
        assert_eq!(seeds["mitochondrion"], vec![mine]);
        assert!(!dir.join("landplants.mitochondrion.fasta").exists());
        fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn one_organelle_still_seeds_the_other_unless_single_seed() {
        let (dir, dir_single) = (tmp("fill-one"), tmp("fill-single"));
        let mut seeds = BTreeMap::new();
        let filled =
            fill_missing(&mut seeds, &dir, Some("mitochondrion"), false, Lineage::Landplants)
                .unwrap();
        assert_eq!(names(&filled), ["mitochondrion", "plastid"]);

        let mut seeds = BTreeMap::new();
        let filled =
            fill_missing(&mut seeds, &dir_single, Some("plastid"), true, Lineage::Landplants)
                .unwrap();
        assert_eq!(names(&filled), ["plastid"]);
        assert!(!seeds.contains_key("mitochondrion"));
        fs::remove_dir_all(&dir).unwrap();
        fs::remove_dir_all(&dir_single).unwrap();
    }

    #[test]
    fn nothing_is_written_when_every_needed_organelle_has_seeds() {
        let dir = tmp("fill-none");
        let mut seeds = BTreeMap::from([
            ("mitochondrion".to_string(), vec![PathBuf::from("m.fa")]),
            ("plastid".to_string(), vec![PathBuf::from("p.fa")]),
        ]);
        assert!(fill_missing(&mut seeds, &dir, None, false, Lineage::Landplants)
            .unwrap()
            .is_empty());
        assert!(!dir.exists());
    }

    #[test]
    fn write_materialises_only_the_requested_organelles() {
        let dir = std::env::temp_dir().join(format!("ovasm-seeddb-test-{}", std::process::id()));
        let written = write(&dir, &["mitochondrion"], Lineage::Landplants).unwrap();
        assert_eq!(written.keys().collect::<Vec<_>>(), ["mitochondrion"]);
        assert!(dir.join("landplants.mitochondrion.fasta").is_file());
        assert!(dir.join("manifest.json").is_file());
        assert!(!dir.join("landplants.plastid.fasta").exists());
        assert!(write(&dir, &["nucleus"], Lineage::Landplants).is_err());
        fs::remove_dir_all(&dir).unwrap();
    }
}
