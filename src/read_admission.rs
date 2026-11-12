//! Exact record admission after ordinary bait recruitment.
use std::collections::BTreeMap;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{ensure, Context, Result};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Deserializer, Serialize};
use sha2::{Digest, Sha256};

use crate::bait::BaitSet;
use crate::recruit::{write_rec, Call, Rec, TargetCounts};

const SCHEMA: &str = "ovasm.read-admission.v1";

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Manifest {
    schema: String,
    scope: String,
    records: Vec<Declaration>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Declaration {
    target: String,
    read_id_hex: String,
    sequence_sha256: String,
    #[serde(deserialize_with = "nullable_hash")]
    quality_sha256: Option<String>,
    provenance: Provenance,
}

fn nullable_hash<'de, D: Deserializer<'de>>(de: D) -> Result<Option<String>, D::Error> {
    Option::<String>::deserialize(de)
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Provenance {
    path: PathBuf,
    sha256: String,
}

#[derive(Clone, Copy, Debug, Serialize)]
#[serde(rename_all = "snake_case")]
enum Status {
    Admitted,
    AlreadyRecruited,
    CompetingTarget,
    Downsampled,
    BaseCap,
}

#[derive(Clone, Debug, Serialize)]
struct RecordReport {
    #[serde(flatten)]
    declaration: Declaration,
    baseline: &'static str,
    status: Status,
}

#[derive(Clone, Debug, Serialize)]
pub(crate) struct Report {
    schema: &'static str,
    manifest: Provenance,
    scope: String,
    policy: &'static str,
    records: Vec<RecordReport>,
}

struct Entry {
    declaration: Declaration,
    target: usize,
    seen: bool,
    baseline: Call,
    original: Option<Rec>,
    status: Option<Status>,
}

pub(crate) struct Admission {
    manifest: Provenance,
    scope: String,
    index: FxHashMap<Vec<u8>, usize>,
    records: Vec<Entry>,
    scan_error: Option<anyhow::Error>,
}

fn hex_digit(b: u8) -> Option<u8> {
    match b {
        b'0'..=b'9' => Some(b - b'0'),
        b'a'..=b'f' => Some(b - b'a' + 10),
        _ => None,
    }
}

fn decode_id(s: &str) -> Option<Vec<u8>> {
    if s.is_empty() || s.len() % 2 != 0 {
        return None;
    }
    s.as_bytes().chunks_exact(2).map(|b| {
        Some(hex_digit(b[0])? << 4 | hex_digit(b[1])?)
    }).collect()
}

fn valid_hash(s: &str) -> bool {
    s.len() == 64 && s.bytes().all(|b| hex_digit(b).is_some())
}

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

struct BoundManifest {
    provenance: Provenance,
    scope: String,
    records: Vec<(Vec<u8>, usize, Declaration)>,
}

fn bind_manifest(path: &Path, bait: &BaitSet, schema: &str, label: &str) -> Result<BoundManifest> {
    let path = path.canonicalize()
        .with_context(|| format!("{label} manifest {}", path.display()))?;
    let bytes = std::fs::read(&path)?;
    let manifest: Manifest = serde_json::from_slice(&bytes)
        .with_context(|| format!("invalid {label} manifest {}", path.display()))?;
    ensure!(manifest.schema == schema, "unsupported {label} schema {:?}", manifest.schema);
    ensure!(!manifest.scope.trim().is_empty(), "{label} requires a scope");
    ensure!(!manifest.records.is_empty(), "{label} manifest has no records");

    let mut records = Vec::new();
    let mut proofs = BTreeMap::new();
    for mut decl in manifest.records {
        let id = decode_id(&decl.read_id_hex)
            .with_context(|| format!("{label} read_id_hex must encode a nonempty full header in lowercase hex"))?;
        let target = bait.targets.iter().position(|t| t.name == decl.target)
            .with_context(|| format!("unknown {label} target {:?}", decl.target))?;
        ensure!(valid_hash(&decl.sequence_sha256)
            && decl.quality_sha256.as_deref().map_or(true, valid_hash)
            && valid_hash(&decl.provenance.sha256),
            "invalid SHA-256 in {label} record {}", decl.read_id_hex);
        let proof = path.parent().expect("canonical manifest has a parent").join(&decl.provenance.path);
        let proof = proof.canonicalize()
            .with_context(|| format!("{label} provenance {}", proof.display()))?;
        let actual = match proofs.entry(proof.clone()) {
            std::collections::btree_map::Entry::Occupied(e) => e.into_mut(),
            std::collections::btree_map::Entry::Vacant(e) => {
                let mut hash = Sha256::new();
                std::io::copy(&mut File::open(&proof)?, &mut hash)?;
                e.insert(format!("{:x}", hash.finalize()))
            }
        };
        ensure!(*actual == decl.provenance.sha256, "{label} provenance SHA-256 mismatch for {}", proof.display());
        decl.provenance.path = proof;
        records.push((id, target, decl));
    }
    Ok(BoundManifest { provenance: Provenance { path, sha256: digest(&bytes) }, scope: manifest.scope, records })
}

impl Admission {
    pub(crate) fn load(path: &Path, bait: &BaitSet) -> Result<Self> {
        let bound = bind_manifest(path, bait, SCHEMA, "read admission")?;
        let mut index = FxHashMap::default();
        let mut records = Vec::new();
        for (id, target, declaration) in bound.records {
            ensure!(index.insert(id, records.len()).is_none(),
                "duplicate read admission manifest read ID {}", declaration.read_id_hex);
            records.push(Entry { declaration, target, seen: false, baseline: Call::Other,
                original: None, status: None });
        }
        Ok(Self { manifest: bound.provenance, scope: bound.scope, index, records, scan_error: None })
    }

    pub(crate) fn check_exclusions(&self, exclusions: &Exclusions) -> Result<()> {
        for (id, &i) in &self.index {
            let entry = &self.records[i];
            ensure!(exclusions.mask(id) & (1 << entry.target) == 0,
                "conflicting admission and source exclusion for read {} and target {:?}",
                entry.declaration.read_id_hex, entry.declaration.target);
        }
        Ok(())
    }

    pub(crate) fn begin_scan(&mut self) {
        for entry in &mut self.records {
            entry.seen = false;
        }
    }

    pub(crate) fn observe(&mut self, rec: &Rec, call: Call, original_input: bool) {
        let Some(&i) = self.index.get(rec.id.as_slice()) else {
            return;
        };
        let entry = &mut self.records[i];
        if self.scan_error.is_some() {
            return;
        }
        let id = &entry.declaration.read_id_hex;
        if entry.seen {
            self.scan_error = Some(anyhow::anyhow!("duplicate read admission ID in original input {id}"));
            return;
        }
        entry.seen = true;
        entry.baseline = call;
        if original_input {
            let decl = &entry.declaration;
            if digest(&rec.seq) != decl.sequence_sha256
                || rec.qual.as_deref().map(digest) != decl.quality_sha256
            {
                self.scan_error = Some(anyhow::anyhow!("read admission sequence or quality SHA-256 mismatch for {id}"));
                return;
            }
            entry.original = Some(rec.clone());
        }
    }

    pub(crate) fn finish_scan(&mut self) -> Result<()> {
        if let Some(error) = self.scan_error.take() {
            return Err(error);
        }
        for entry in &self.records {
            ensure!(entry.seen && entry.original.is_some(),
                "read admission ID missing from input {}", entry.declaration.read_id_hex);
        }
        Ok(())
    }

    pub(crate) fn append(
        &mut self,
        target: usize,
        writer: &mut impl Write,
        eligible: impl Fn(&Rec) -> bool,
        available_bases: u64,
    ) -> Result<TargetCounts> {
        let mut added = TargetCounts::default();
        for entry in self.records.iter_mut().filter(|entry| entry.target == target) {
            let rec = entry.original.as_ref().expect("validated original record");
            let bases = rec.seq.len() as u64;
            let status = match entry.baseline {
                Call::Organelle(mask) if mask & (1 << target) != 0 => Status::AlreadyRecruited,
                Call::Organelle(_) => Status::CompetingTarget,
                _ if !eligible(rec) => Status::Downsampled,
                _ if bases > available_bases - added.bases => Status::BaseCap,
                _ => {
                    write_rec(writer, rec)?;
                    added.reads += 1;
                    added.bases += bases;
                    Status::Admitted
                }
            };
            entry.status = Some(status);
        }
        Ok(added)
    }

    pub(crate) fn report(self) -> Report {
        Report {
            schema: SCHEMA,
            manifest: self.manifest,
            scope: self.scope,
            policy: "externally supplied source evidence; only explicitly bound original records can enter final reads; no bait feedback",
            records: self.records.into_iter().map(|entry| RecordReport {
                declaration: entry.declaration,
                baseline: match entry.baseline {
                    Call::Other => "other",
                    Call::Organelle(_) => "organelle",
                    Call::NumtCandidate(_) => "partial",
                },
                status: entry.status.expect("all targets were written"),
            }).collect(),
        }
    }
}


const EXCLUSION_SCHEMA: &str = "ovasm.read-source-exclusions.v1";

#[derive(Debug, Clone, Serialize)]
struct ExcludedRecord {
    #[serde(flatten)]
    declaration: Declaration,
    #[serde(skip)]
    target: usize,
    organelle_match_passes: Vec<usize>,
    partial_match_passes: Vec<usize>,
}

#[derive(Default)]
struct ExcludedOriginal {
    mask: u8,
    rows: Vec<usize>,
    seen: bool,
}

#[derive(Debug, Clone, Serialize)]
pub(crate) struct ExclusionReport {
    schema: &'static str,
    manifest: Provenance,
    scope: String,
    policy: &'static str,
    library_walk_skipped_targets: Vec<String>,
    records: Vec<ExcludedRecord>,
}

pub(crate) struct Exclusions {
    manifest: Provenance,
    scope: String,
    originals: FxHashMap<Vec<u8>, ExcludedOriginal>,
    records: Vec<ExcludedRecord>,
    targets: u8,
}

impl Exclusions {
    pub(crate) fn load(path: &Path, bait: &BaitSet) -> Result<Self> {
        let bound = bind_manifest(path, bait, EXCLUSION_SCHEMA, "read source exclusion")?;
        let mut originals = FxHashMap::<Vec<u8>, ExcludedOriginal>::default();
        let mut records: Vec<ExcludedRecord> = Vec::new();
        let mut targets = 0;
        for (id, target, declaration) in bound.records {
            let entry = originals.entry(id).or_default();
            let bit = 1 << target;
            ensure!(entry.mask & bit == 0, "duplicate source exclusion for read {} and target {:?}",
                declaration.read_id_hex, declaration.target);
            if let Some(&i) = entry.rows.first() {
                let original = &records[i].declaration;
                ensure!(original.sequence_sha256 == declaration.sequence_sha256
                    && original.quality_sha256 == declaration.quality_sha256,
                    "conflicting original hashes for read source exclusion {}", declaration.read_id_hex);
            }
            entry.rows.push(records.len());
            entry.mask |= bit;
            targets |= bit;
            records.push(ExcludedRecord { declaration, target, organelle_match_passes: Vec::new(), partial_match_passes: Vec::new() });
        }
        Ok(Self { manifest: bound.provenance, scope: bound.scope, originals, records, targets })
    }

    pub(crate) fn mask(&self, id: &[u8]) -> u8 {
        self.originals.get(id).map_or(0, |entry| entry.mask)
    }

    pub(crate) fn protects(&self, target: usize) -> bool {
        self.targets & (1 << target) != 0
    }

    pub(crate) fn begin_scan(&mut self) {
        for entry in self.originals.values_mut() {
            entry.seen = false;
        }
    }

    /// Replay omits quality bytes. The original-input scan before replay and
    /// the final original-output scan both verify the exact quality digest.
    pub(crate) fn observe(&mut self, rec: &Rec, call: Call, original_input: bool, pass: Option<usize>) -> Result<u8> {
        let Some(entry) = self.originals.get_mut(rec.id.as_slice()) else { return Ok(0); };
        let decl = &self.records[entry.rows[0]].declaration;
        ensure!(!entry.seen, "duplicate read source exclusion ID in input {}", decl.read_id_hex);
        ensure!(digest(&rec.seq) == decl.sequence_sha256
            && (!original_input || rec.qual.as_deref().map(digest) == decl.quality_sha256),
            "read source exclusion sequence or quality SHA-256 mismatch for {}", decl.read_id_hex);
        entry.seen = true;
        if let Some(pass) = pass {
            for &i in &entry.rows {
                let row = &mut self.records[i];
                match call {
                    Call::Organelle(mask) if mask & (1 << row.target) != 0 => row.organelle_match_passes.push(pass),
                    Call::NumtCandidate(target) if target == row.target => row.partial_match_passes.push(pass),
                    _ => {}
                }
            }
        }
        Ok(entry.mask)
    }

    pub(crate) fn finish_scan(&self) -> Result<()> {
        for entry in self.originals.values() {
            ensure!(entry.seen, "read source exclusion ID missing from input {}",
                self.records[entry.rows[0]].declaration.read_id_hex);
        }
        Ok(())
    }

    pub(crate) fn report(self, bait: &BaitSet, kmer_gate: bool) -> ExclusionReport {
        let library_walk_skipped_targets = bait.targets.iter().enumerate()
            .filter(|(i, _)| kmer_gate && self.protects(*i))
            .map(|(_, t)| t.name.clone()).collect();
        ExclusionReport { schema: EXCLUSION_SCHEMA, manifest: self.manifest, scope: self.scope,
            policy: "prohibited target membership and read-derived bait after initial seed loading; original diagnostics remain available",
            library_walk_skipped_targets, records: self.records }
    }
}
