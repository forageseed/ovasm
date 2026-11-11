//! Verify a user-supplied seed source audit before using its retained FASTA in a run.

use anyhow::{ensure, Context, Result};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::fs::{self, OpenOptions};
use std::io::{ErrorKind, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

const SCHEMA: &str = "ovasm.seed-source-audit.v1";

#[derive(Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
struct Counts {
    target_supported: usize,
    competing_supported: usize,
    unresolved: usize,
}

#[derive(Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
enum State {
    TargetSupported,
    CompetingSupported,
    Unresolved,
}

#[derive(Deserialize)]
struct Row {
    state: State,
    read_id_hex: String,
    seed_sha256: String,
}

#[derive(Deserialize)]
struct Audit {
    schema: String,
    target: String,
    scope: String,
    manifest_sha256: String,
    retained_fasta_sha256: String,
    retained_records: usize,
    counts: Counts,
    records: Vec<Row>,
}

#[derive(Serialize)]
pub(crate) struct Provenance {
    report: PathBuf,
    report_sha256: String,
    manifest_sha256: String,
    retained_fasta: PathBuf,
    retained_fasta_sha256: String,
    scope: String,
    counts: Counts,
}

pub(crate) struct VerifiedView {
    audit: Audit,
    report: Vec<u8>,
    fasta: Vec<u8>,
}

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn is_sha256(s: &str) -> bool {
    s.len() == 64 && s.bytes().all(|b| matches!(b, b'0'..=b'9' | b'a'..=b'f'))
}

fn verify_bytes(report: Vec<u8>, fasta: Vec<u8>, target: &str) -> Result<VerifiedView> {
    let audit: Audit = serde_json::from_slice(&report).context("reading seed source audit")?;
    ensure!(audit.schema == SCHEMA, "unrecognized seed source audit schema");
    ensure!(audit.target == target, "seed source audit target is {}, expected {target}", audit.target);
    ensure!(!audit.scope.trim().is_empty() && is_sha256(&audit.manifest_sha256),
            "seed source audit has no bound manifest and scope");
    ensure!(audit.retained_fasta_sha256 == sha256(&fasta), "seed source view FASTA SHA-256 mismatch");
    ensure!(fasta.starts_with(b">"), "seed source view must contain a nonempty plain FASTA");

    let mut counts = Counts::default();
    let mut retained = 0;
    let mut reader = needletail::parse_fastx_reader(fasta.as_slice())?;
    for row in &audit.records {
        match row.state {
            State::CompetingSupported => {
                counts.competing_supported += 1;
                continue;
            }
            State::TargetSupported => counts.target_supported += 1,
            State::Unresolved => counts.unresolved += 1,
        }
        let record = reader.next().context("seed source view omits an audited retained crop")??;
        let name = record.id().split(u8::is_ascii_whitespace).next().unwrap_or_default();
        let name_hex: String = name.iter().map(|byte| format!("{byte:02x}")).collect();
        ensure!(!name.is_empty() && name_hex == row.read_id_hex
                && sha256(&record.seq().to_ascii_uppercase()) == row.seed_sha256,
                "seed source view differs from audited retained crop {}", retained + 1);
        retained += 1;
    }
    ensure!(reader.next().is_none(), "seed source view includes an unaudited crop");
    ensure!(counts == audit.counts && retained == audit.retained_records && retained > 0,
            "seed source audit counts do not match its retained crops");
    drop(reader);
    Ok(VerifiedView { audit, report, fasta })
}

pub(crate) fn verify(report: &Path, fasta: &Path, target: &str) -> Result<VerifiedView> {
    let report_bytes = fs::read(report).with_context(|| format!("reading {}", report.display()))?;
    let fasta_bytes = fs::read(fasta).with_context(|| format!("reading {}", fasta.display()))?;
    verify_bytes(report_bytes, fasta_bytes, target).with_context(|| {
        format!("checking {} against {}", report.display(), fasta.display())
    })
}

fn write_artifact(path: &Path, bytes: &[u8]) -> Result<()> {
    static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);
    let temp = path.with_file_name(format!(
        ".ovasm-source-seed-{}-{}.tmp",
        std::process::id(), NEXT_TEMP.fetch_add(1, Ordering::Relaxed),
    ));
    let mut file = OpenOptions::new().write(true).create_new(true).open(&temp)?;
    let result = (|| -> Result<()> {
        file.write_all(bytes)?;
        drop(file);
        match fs::hard_link(&temp, path) {
            Ok(()) => Ok(()),
            Err(error) if error.kind() == ErrorKind::AlreadyExists => {
                ensure!(fs::symlink_metadata(path)?.is_file() && fs::read(path)? == bytes,
                        "existing seed source artifact differs: {}", path.display());
                Ok(())
            }
            Err(error) => Err(error.into()),
        }
    })();
    let cleanup = fs::remove_file(&temp);
    result.with_context(|| format!("writing {}", path.display()))?;
    cleanup?;
    Ok(())
}

impl VerifiedView {
    pub(crate) fn materialize(self, output: &Path) -> Result<(PathBuf, Provenance)> {
        let dir = output.join("source-seeds");
        fs::create_dir_all(&dir)?;
        let report_sha256 = sha256(&self.report);
        let report = dir.join(format!("{}.{}.audit.json", self.audit.target, report_sha256));
        let fasta = dir.join(format!("{}.{}.fasta", self.audit.target, self.audit.retained_fasta_sha256));
        write_artifact(&fasta, &self.fasta)?;
        write_artifact(&report, &self.report)?;
        Ok((fasta.clone(), Provenance {
            report,
            report_sha256,
            manifest_sha256: self.audit.manifest_sha256,
            retained_fasta: fasta,
            retained_fasta_sha256: self.audit.retained_fasta_sha256,
            scope: self.audit.scope,
            counts: self.audit.counts,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    fn fixture(fasta: &[u8]) -> Vec<u8> {
        serde_json::to_vec(&json!({
            "schema": SCHEMA,
            "target": "mitochondrion",
            "scope": "frozen source test",
            "manifest_sha256": sha256(b"manifest"),
            "retained_fasta_sha256": sha256(fasta),
            "retained_records": 2,
            "counts": {"target_supported": 1, "competing_supported": 1, "unresolved": 1},
            "records": [
                {"state": "competing_supported", "read_id_hex": "62", "seed_sha256": sha256(b"TGCA")},
                {"state": "target_supported", "read_id_hex": "61", "seed_sha256": sha256(b"ACGT")},
                {"state": "unresolved", "read_id_hex": "61", "seed_sha256": sha256(b"ACGA")},
            ],
        })).unwrap()
    }

    #[test]
    fn retained_order_accepts_unresolved_and_repeated_identifiers() {
        let fasta = b">a span=0-4\r\nac\r\ngt\r\n>a span=4-8\r\nACGA\r\n";
        let view = verify_bytes(fixture(fasta), fasta.to_vec(), "mitochondrion").unwrap();
        assert_eq!(view.fasta, fasta);
    }

    #[test]
    fn changing_the_file_digest_does_not_authorize_dropping_an_unresolved_crop() {
        let fasta = b">a span=0-4\nACGT\n";
        let result = verify_bytes(fixture(fasta), fasta.to_vec(), "mitochondrion");
        assert!(result.err().unwrap().to_string().contains("omits"));
    }

    #[test]
    fn crop_reordering_or_wrong_target_are_rejected() {
        let reordered = b">a span=4-8\nACGA\n>a span=0-4\nACGT\n";
        let result = verify_bytes(fixture(reordered), reordered.to_vec(), "mitochondrion");
        assert!(result.err().unwrap().to_string().contains("differs"));
        let result = verify_bytes(fixture(reordered), reordered.to_vec(), "plastid");
        assert!(result.err().unwrap().to_string().contains("target"));
    }
}
