//! Verify a user-supplied seed source audit before using its retained FASTA in a run.

use anyhow::{ensure, Context, Result};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::BTreeMap;
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

#[derive(Clone, Copy, Default, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
enum RetentionPolicy {
    #[default]
    ExcludeCompeting,
    TargetSupportedOnly,
}

impl RetentionPolicy {
    fn is_default(&self) -> bool {
        *self == Self::ExcludeCompeting
    }

    fn retains(self, state: State) -> bool {
        match self {
            Self::ExcludeCompeting => state != State::CompetingSupported,
            Self::TargetSupportedOnly => state == State::TargetSupported,
        }
    }
}

type Position = (usize, usize, String);

#[derive(Deserialize)]
struct Observation {
    window: String,
    role: String,
    original_positions: Vec<Position>,
    inside_crop: Vec<Position>,
    qualifies: bool,
}

#[derive(Deserialize)]
struct Window {
    interval: [usize; 2],
    boundary: usize,
    role: String,
    sequence_sha256: String,
    provenance_sha256: String,
    unique_in_named_panel: bool,
}

#[derive(Deserialize)]
struct RawInput {
    path: PathBuf,
    sha256: String,
}

#[derive(Deserialize)]
struct Inputs {
    raw_reads: Vec<RawInput>,
}

#[derive(Deserialize)]
struct Row {
    state: State,
    read_id_hex: String,
    seed_sha256: String,
    crop: Option<[usize; 2]>,
    original_input_index: Option<usize>,
    original_record_index: Option<usize>,
    original_sha256: Option<String>,
    #[serde(default)]
    observations: Vec<Observation>,
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
    #[serde(default)]
    retention_policy: RetentionPolicy,
    inputs: Option<Inputs>,
    #[serde(default)]
    windows: BTreeMap<String, Window>,
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
    #[serde(skip_serializing_if = "RetentionPolicy::is_default")]
    retention_policy: RetentionPolicy,
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

fn crop_span(header: &[u8]) -> Option<[usize; 2]> {
    let mut spans = header.split(u8::is_ascii_whitespace).skip(1)
        .filter_map(|field| field.strip_prefix(b"span="));
    let (start, end) = std::str::from_utf8(spans.next()?).ok()?.split_once('-')?;
    if spans.next().is_some() || [start, end].into_iter().any(|n| {
        n.is_empty() || !n.bytes().all(|b| b.is_ascii_digit())
    }) {
        return None;
    }
    Some([start.parse().ok()?, end.parse().ok()?])
}

fn verify_strict_crop(audit: &Audit, row: &Row, header: &[u8], sequence: &[u8]) -> Result<()> {
    let crop = row.crop.context("strict seed source crop has no audited span")?;
    let [start, end] = crop;
    ensure!(start < end && end - start == sequence.len() && crop_span(header) == Some(crop),
            "strict seed source crop span differs from the audited original coordinates");
    let source = audit.inputs.as_ref().and_then(|inputs| inputs.raw_reads.get(row.original_input_index?))
        .context("strict seed source crop has no bound raw input")?;
    ensure!(!source.path.as_os_str().is_empty() && is_sha256(&source.sha256)
            && row.original_record_index.is_some() && row.original_sha256.as_deref().is_some_and(is_sha256),
            "strict seed source crop has no bound original record");
    ensure!(row.observations.iter().all(|observation| observation.role == "target"
        && audit.windows.get(&observation.window).is_some_and(|window| window.role == observation.role)),
        "strict seed source crop has a conflicting or invalid source role");
    let mut supported = false;
    for observation in row.observations.iter().filter(|observation| observation.qualifies) {
        let window = audit.windows.get(&observation.window).context("strict seed source crop names an unknown window")?;
        let [position] = observation.original_positions.as_slice() else {
            anyhow::bail!("strict seed source window has no unique original placement");
        };
        let (lo, hi, strand) = position;
        ensure!(observation.inside_crop.as_slice() == [position.clone()]
                && start <= *lo && lo < hi && *hi <= end,
                "strict seed source window is not fully inside the crop");
        let [ref_lo, ref_hi] = window.interval;
        ensure!(ref_lo < window.boundary && window.boundary < ref_hi && ref_hi - ref_lo == hi - lo
                && window.role == "target" && window.unique_in_named_panel
                && is_sha256(&window.provenance_sha256) && is_sha256(&window.sequence_sha256),
                "strict seed source window lacks unique target evidence");
        let word = &sequence[lo - start..hi - start];
        let reverse: Vec<_> = word.iter().rev().map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        }).collect();
        let oriented = match strand.as_str() {
            "+" => word,
            "-" => reverse.as_slice(),
            _ => anyhow::bail!("invalid strict seed source strand"),
        };
        ensure!(!reverse.contains(&b'N') && reverse != word && sha256(oriented) == window.sequence_sha256
                && sequence.windows(word.len()).filter(|&slice| slice == word || slice == reverse).take(2).count() == 1,
                "strict seed source crop does not uniquely spell its audited window");
        supported = true;
    }
    ensure!(supported, "strict seed source crop has no qualifying target window");
    Ok(())
}

fn verify_bytes(report: Vec<u8>, fasta: Vec<u8>, target: &str) -> Result<VerifiedView> {
    let audit: Audit = serde_json::from_slice(&report).context("reading seed source audit")?;
    ensure!(audit.schema == SCHEMA, "unrecognized seed source audit schema");
    ensure!(audit.target == target, "seed source audit target is {}, expected {target}", audit.target);
    ensure!(!audit.scope.trim().is_empty() && is_sha256(&audit.manifest_sha256),
            "seed source audit has no bound manifest and scope");
    ensure!(audit.retained_fasta_sha256 == sha256(&fasta), "seed source view FASTA SHA-256 mismatch");
    if audit.retention_policy == RetentionPolicy::TargetSupportedOnly {
        ensure!(audit.retained_records > 0, "target-supported-only audit retained zero crops for {target}");
    }
    ensure!(fasta.starts_with(b">"), "seed source view must contain a nonempty plain FASTA");

    let mut counts = Counts::default();
    let mut retained = 0;
    let mut reader = needletail::parse_fastx_reader(fasta.as_slice())?;
    for row in &audit.records {
        match row.state {
            State::CompetingSupported => counts.competing_supported += 1,
            State::TargetSupported => counts.target_supported += 1,
            State::Unresolved => counts.unresolved += 1,
        }
        if !audit.retention_policy.retains(row.state) {
            continue;
        }
        let record = reader.next().context("seed source view omits an audited retained crop")??;
        let name = record.id().split(u8::is_ascii_whitespace).next().unwrap_or_default();
        let name_hex: String = name.iter().map(|byte| format!("{byte:02x}")).collect();
        let sequence = record.seq().to_ascii_uppercase();
        ensure!(!name.is_empty() && name_hex == row.read_id_hex && sha256(&sequence) == row.seed_sha256,
                "seed source view differs from audited retained crop {}", retained + 1);
        if audit.retention_policy == RetentionPolicy::TargetSupportedOnly {
            verify_strict_crop(&audit, row, record.id(), &sequence)?;
        }
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
            retention_policy: self.audit.retention_policy,
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

    fn strict_fixture(fasta: &[u8]) -> serde_json::Value {
        json!({
            "schema": SCHEMA, "target": "mitochondrion",
            "scope": "one strict crop and two exclusions",
            "manifest_sha256": sha256(b"manifest"),
            "retained_fasta_sha256": sha256(fasta),
            "retention_policy": "target_supported_only", "retained_records": 1,
            "counts": {"target_supported": 1, "competing_supported": 1, "unresolved": 1},
            "inputs": {"raw_reads": [{"path":"originals.fastq","sha256":sha256(b"raw")}]},
            "windows": {"mt": {"role":"target","interval":[1,5],"boundary":3,
                "sequence_sha256":sha256(b"GACT"),"provenance_sha256":sha256(b"proof"),
                "unique_in_named_panel":true}},
            "records": [
                {"state":"competing_supported","read_id_hex":"62","seed_sha256":sha256(b"TGCA")},
                {"state":"unresolved","read_id_hex":"61","seed_sha256":sha256(b"ACGA")},
                {"state":"target_supported","read_id_hex":"61","seed_sha256":sha256(b"AGTCCT"),
                 "crop":[10,16],"original_input_index":0,"original_record_index":1,"original_sha256":sha256(b"original"),
                 "observations":[{"role":"target","window":"mt","original_positions":[[10,14,"-"]],
                    "inside_crop":[[10,14,"-"]],"qualifies":true}]}
            ],
        })
    }

    fn strict_check(value: serde_json::Value, fasta: &[u8]) -> Result<VerifiedView> {
        verify_bytes(serde_json::to_vec(&value).unwrap(), fasta.to_vec(), "mitochondrion")
    }

    #[test]
    fn strict_source_retains_only_the_bound_reverse_window() {
        let fasta = b">a span=10-16\r\nagtcct\r\n";
        let view = strict_check(strict_fixture(fasta), fasta).unwrap();
        assert_eq!(view.fasta, fasta);
        assert!(view.audit.retention_policy == RetentionPolicy::TargetSupportedOnly);
    }

    #[test]
    fn strict_source_rejects_shifted_or_repeated_spans_even_with_matching_digest() {
        for fasta in [
            b">a span=11-17\nAGTCCT\n".as_slice(),
            b">a span=10-16 span=10-16\nAGTCCT\n".as_slice(),
            b">a\nAGTCCT\n".as_slice(),
        ] {
            let error = strict_check(strict_fixture(fasta), fasta).err().unwrap();
            assert!(error.to_string().contains("span"), "{error}");
        }
    }

    #[test]
    fn strict_source_rejects_an_unresolved_extra_record_and_an_empty_view() {
        let fasta = b">a span=10-16\nAGTCCT\n>a span=4-8\nACGA\n";
        let error = strict_check(strict_fixture(fasta), fasta).err().unwrap();
        assert!(error.to_string().contains("unaudited crop"), "{error}");
        let mut empty = strict_fixture(b"");
        empty["retained_records"] = json!(0);
        let error = strict_check(empty, b"").err().unwrap();
        assert!(error.to_string().contains("zero crops"), "{error}");
    }

    #[test]
    fn strict_source_rejects_invalid_original_or_literal_window_evidence() {
        let fasta = b">a span=10-16\nAGTCCT\n";
        let mut bad_input = strict_fixture(fasta);
        bad_input["records"][2]["original_input_index"] = json!(2);
        let error = strict_check(bad_input, fasta).err().unwrap();
        assert!(error.to_string().contains("bound raw input"), "{error}");
        let mut wrong_word = strict_fixture(fasta);
        wrong_word["windows"]["mt"]["sequence_sha256"] = json!(sha256(b"AGTC"));
        let error = strict_check(wrong_word, fasta).err().unwrap();
        assert!(error.to_string().contains("spell"), "{error}");
        let mut mixed = strict_fixture(fasta);
        mixed["records"][2]["observations"].as_array_mut().unwrap().push(
            json!({"role":"competing","window":"nuc","original_positions":[[0,4,"+"]],
                   "inside_crop":[],"qualifies":false}));
        let mut competing = mixed["windows"]["mt"].clone();
        competing["role"] = json!("competing");
        mixed["windows"]["nuc"] = competing;
        let error = strict_check(mixed.clone(), fasta).err().unwrap();
        assert!(error.to_string().contains("conflicting"), "{error}");
        mixed["records"][2]["observations"][1]["role"] = json!("target");
        let error = strict_check(mixed, fasta).err().unwrap();
        assert!(error.to_string().contains("conflicting"), "{error}");
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
