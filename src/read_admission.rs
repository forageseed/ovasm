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

impl Admission {
    pub(crate) fn load(path: &Path, bait: &BaitSet) -> Result<Self> {
        let path = path.canonicalize()
            .with_context(|| format!("read admission manifest {}", path.display()))?;
        let bytes = std::fs::read(&path)?;
        let manifest: Manifest = serde_json::from_slice(&bytes)
            .with_context(|| format!("invalid read admission manifest {}", path.display()))?;
        ensure!(manifest.schema == SCHEMA, "unsupported read admission schema {:?}", manifest.schema);
        ensure!(!manifest.scope.trim().is_empty(), "read admission requires a scope");
        ensure!(!manifest.records.is_empty(), "read admission manifest has no records");

        let mut index = FxHashMap::default();
        let mut records = Vec::new();
        let mut proofs = BTreeMap::new();
        for mut decl in manifest.records {
            let id = decode_id(&decl.read_id_hex)
                .context("read admission read_id_hex must encode a nonempty full header in lowercase hex")?;
            ensure!(index.insert(id, records.len()).is_none(),
                "duplicate read admission manifest read ID {}", decl.read_id_hex);
            let target = bait.targets.iter().position(|t| t.name == decl.target)
                .with_context(|| format!("unknown read admission target {:?}", decl.target))?;
            ensure!(valid_hash(&decl.sequence_sha256)
                && decl.quality_sha256.as_deref().map_or(true, valid_hash)
                && valid_hash(&decl.provenance.sha256),
                "invalid SHA-256 in read admission record {}", decl.read_id_hex);
            let proof = path.parent().expect("canonical manifest has a parent").join(&decl.provenance.path);
            let proof = proof.canonicalize()
                .with_context(|| format!("read admission provenance {}", proof.display()))?;
            let actual = match proofs.entry(proof.clone()) {
                std::collections::btree_map::Entry::Occupied(e) => e.into_mut(),
                std::collections::btree_map::Entry::Vacant(e) => {
                    let mut hash = Sha256::new();
                    std::io::copy(&mut File::open(&proof)?, &mut hash)?;
                    e.insert(format!("{:x}", hash.finalize()))
                }
            };
            ensure!(*actual == decl.provenance.sha256,
                "read admission provenance SHA-256 mismatch for {}", proof.display());
            decl.provenance.path = proof;
            records.push(Entry {
                declaration: decl,
                target,
                seen: false,
                baseline: Call::Other,
                original: None,
                status: None,
            });
        }
        Ok(Self {
            manifest: Provenance { path, sha256: digest(&bytes) },
            scope: manifest.scope,
            index,
            records,
            scan_error: None,
        })
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
