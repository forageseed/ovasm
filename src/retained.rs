//! Retained source candidates. Bait compatibility and exact homology are not source certification.
use crate::bait::BaitSet;
use crate::evidence::{parse_gfa, revcomp, Graph};
use crate::kmer::CanonicalKmers;
use crate::recruit::{classify, for_each_rec, write_rec, Call, ClassifyParams, Rec};
use anyhow::Result;
use rustc_hash::{FxHashMap, FxHashSet, FxHasher};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::hash::Hasher;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

fn save(path: &Path, value: &Value) -> Result<()> {
    fs::write(path, serde_json::to_string_pretty(value)? + "\n")?;
    Ok(())
}

/// Per-target compatible blocks use precisely the recruitment gap rule.
pub(crate) fn blocks(seq: &[u8], bait: &BaitSet, gap: usize) -> Vec<Vec<(usize, usize)>> {
    let mut result: Vec<Vec<(usize, usize)>> = vec![Vec::new(); bait.targets.len()];
    for (pos, key) in CanonicalKmers::new(seq, bait.k) {
        let mask = bait.get(key);
        for (t, list) in result.iter_mut().enumerate() {
            if mask & (1 << t) == 0 {
                continue;
            }
            if let Some(last) = list.last_mut() {
                if pos <= last.1 + gap {
                    last.1 = last.1.max(pos + bait.k);
                    continue;
                }
            }
            list.push((pos, pos + bait.k));
        }
    }
    result
}

const FRONTIER_AUDIT_EXAMPLES: usize = 8;
const FRONTIER_AUDIT_ID_BYTES: usize = 1024;
const FRONTIER_AUDIT_SEQUENCE_BYTES: usize = 65536;
const FRONTIER_AUDIT_INTERVALS: usize = 16;

#[derive(Clone, Copy)]
pub(crate) enum FrontierReason {
    ExtensionSubsampled,
    KmerDepthGate,
    ReadDepthProfile,
    BelowTargetDepth,
    AboveTargetDepth,
    UnbaitedAfterCommit,
}

impl FrontierReason {
    const ALL: [Self; 6] = [
        Self::ExtensionSubsampled,
        Self::KmerDepthGate,
        Self::ReadDepthProfile,
        Self::BelowTargetDepth,
        Self::AboveTargetDepth,
        Self::UnbaitedAfterCommit,
    ];

    fn name(self) -> &'static str {
        match self {
            Self::ExtensionSubsampled => "extension_subsampled",
            Self::KmerDepthGate => "kmer_depth_gate",
            Self::ReadDepthProfile => "read_depth_profile",
            Self::BelowTargetDepth => "below_target_depth",
            Self::AboveTargetDepth => "above_target_depth",
            Self::UnbaitedAfterCommit => "unbaited_after_commit",
        }
    }
}

struct FrontierWitness {
    key: u64,
    id: Vec<u8>,
    sequence: Vec<u8>,
    pass: usize,
    in_final_partial_pool: bool,
}

#[derive(Default)]
struct FrontierEvents {
    extension_records: u64,
    within_pass_records: u64,
    oversized_records: u64,
    witnesses: Vec<FrontierWitness>,
}

pub(crate) struct FrontierAudit {
    events: Vec<[FrontierEvents; FrontierReason::ALL.len()]>,
    witness_index: FxHashMap<Vec<u8>, Vec<(usize, usize, usize)>>,
}

impl FrontierAudit {
    pub(crate) fn new(targets: usize) -> Self {
        Self {
            events: (0..targets).map(|_| Default::default()).collect(),
            witness_index: FxHashMap::default(),
        }
    }

    pub(crate) fn record(
        &mut self,
        target: usize,
        pass: usize,
        reason: FrontierReason,
        name_key: u64,
        read: &Rec,
    ) {
        let events = &mut self.events[target][reason as usize];
        events.extension_records += 1;
        if read.id.len() > FRONTIER_AUDIT_ID_BYTES
            || read.seq.len() > FRONTIER_AUDIT_SEQUENCE_BYTES
        {
            events.oversized_records += 1;
            return;
        }
        let rank = (name_key, read.id.as_slice(), read.seq.as_slice());
        let position = events.witnesses.binary_search_by(|witness| {
            (witness.key, witness.id.as_slice(), witness.sequence.as_slice()).cmp(&rank)
        });
        let Err(position) = position else { return };
        if position >= FRONTIER_AUDIT_EXAMPLES {
            return;
        }
        if events.witnesses.len() == FRONTIER_AUDIT_EXAMPLES {
            events.witnesses.pop();
        }
        events.witnesses.insert(position, FrontierWitness {
            key: name_key,
            id: read.id.clone(),
            sequence: read.seq.clone(),
            pass,
            in_final_partial_pool: false,
        });
    }

    pub(crate) fn count(
        &mut self, target: usize, reason: FrontierReason, within_pass: bool, records: u64,
    ) {
        let events = &mut self.events[target][reason as usize];
        if within_pass {
            events.within_pass_records += records;
        } else {
            events.extension_records += records;
        }
    }

    fn prepare_pool_membership(&mut self) {
        for (t, reasons) in self.events.iter().enumerate() {
            for (r, events) in reasons.iter().enumerate() {
                for (w, witness) in events.witnesses.iter().enumerate() {
                    self.witness_index.entry(witness.id.clone()).or_default().push((t, r, w));
                }
            }
        }
    }

    fn match_pool_record(&mut self, read: &Rec) {
        if let Some(indices) = self.witness_index.get(&read.id) {
            for &(t, r, w) in indices {
                let witness = &mut self.events[t][r].witnesses[w];
                if witness.sequence == read.seq {
                    witness.in_final_partial_pool = true;
                }
            }
        }
    }

    fn write(&self, dir: &Path, bait: &BaitSet, gap: usize) -> Result<Value> {
        let path = dir.join("frontier-evidence.jsonl");
        let mut rows = BufWriter::new(File::create(&path)?);
        let mut targets = Vec::new();
        for (t, reasons) in self.events.iter().enumerate() {
            let name = &bait.targets[t].name;
            let mut summary = Vec::new();
            for (reason, events) in FrontierReason::ALL.into_iter().zip(reasons) {
                for w in &events.witnesses {
                    let compatible = blocks(&w.sequence, bait, gap);
                    let intervals = &compatible[t];
                    let shown = intervals.len().min(FRONTIER_AUDIT_INTERVALS);
                    let id = String::from_utf8_lossy(&w.id);
                    let exact_id = std::str::from_utf8(&w.id).is_err().then(|| {
                        w.id.iter().map(|byte| format!("{byte:02x}")).collect::<String>()
                    });
                    let mut fingerprint = FxHasher::default();
                    for &base in &w.sequence {
                        fingerprint.write_u8(base);
                    }
                    writeln!(rows, "{}", json!({
                        "target": name, "reason": reason.name(), "pass": w.pass,
                        "read_id": id, "read_id_hex_if_non_utf8": exact_id,
                        "read_length": w.sequence.len(),
                        "read_sequence_fxhash64": format!("{:016x}", fingerprint.finish()),
                        "in_final_partial_pool": w.in_final_partial_pool,
                        "final_bait_compatible_blocks": &intervals[..shown],
                        "final_bait_intervals_omitted": intervals.len() - shown,
                        "flank_origin": "unknown"
                    }))?;
                }
                summary.push(json!({
                    "reason": reason.name(), "extension_records": events.extension_records,
                    "within_pass_records": events.within_pass_records,
                    "oversized_extension_records": events.oversized_records,
                    "sampled": events.witnesses.len()
                }));
            }
            targets.push(json!({"target": name, "reasons": summary}));
        }
        rows.flush()?;
        Ok(json!({
            "schema": "ovasm.recruitment-frontier-evidence.v1",
            "read_evidence": path,
            "scope": "held or unevaluated records in configured extension; source origin is unknown",
            "counting": "record events by stage and reason, possibly repeated across passes",
            "sampling": "bounded main-extension examples by name hash, raw ID and exact sequence; within-pass and post-commit records are aggregate only",
            "extension_subsampled": "not evaluated for extension",
            "pool_matching": "full read ID and sequence, not an independent molecule identity",
            "coordinates": "0-based half-open intervals of final-bait compatibility, not alignments",
            "max_examples_per_target_reason": FRONTIER_AUDIT_EXAMPLES,
            "max_sample_id_bytes": FRONTIER_AUDIT_ID_BYTES,
            "max_sample_sequence_bytes": FRONTIER_AUDIT_SEQUENCE_BYTES,
            "max_intervals_per_example": FRONTIER_AUDIT_INTERVALS,
            "targets": targets
        }))
    }
}

pub(crate) fn recruitment(
    out: &Path,
    pool: &Path,
    bait: &BaitSet,
    params: &ClassifyParams,
    mut frontier: Option<&mut FrontierAudit>,
    exclusions: Option<&crate::read_admission::Exclusions>,
) -> Result<Value> {
    let start = Instant::now();
    let dir = out.join("source_candidates");
    fs::create_dir_all(&dir)?;
    let mut writers = bait
        .targets
        .iter()
        .map(|target| {
            Ok(BufWriter::new(File::create(
                dir.join(format!("{}.fastq", target.name)),
            )?))
        })
        .collect::<Result<Vec<_>>>()?;
    let mut counts = vec![0u64; bait.targets.len()];
    let mut rows = BufWriter::new(File::create(
        dir.join("nuclear-insertion-candidates.jsonl"),
    )?);
    if let Some(audit) = frontier.as_deref_mut() {
        audit.prepare_pool_membership();
    }
    let mut failure = Ok(());
    for_each_rec(pool, |rec| {
        if failure.is_err() {
            return;
        }
        failure = (|| -> Result<()> {
            if let Some(audit) = frontier.as_deref_mut() {
                audit.match_pool_record(&rec);
            }
            let (call, score) = classify(&rec.seq, bait, params);
            let compatible = blocks(&rec.seq, bait, params.block_gap);
            let mut targets = Vec::new();
            for (t, intervals) in compatible.into_iter().enumerate() {
                let longest = intervals.iter().map(|&(a, b)| b - a).max().unwrap_or(0);
                if longest < params.numt_min_block {
                    continue;
                }
                write_rec(&mut writers[t], &rec)?;
                counts[t] += 1;
                targets.push(json!({"target": bait.targets[t].name, "compatible_blocks": intervals,
                    "longest_block": longest,
                    "candidate_label": match bait.targets[t].name.as_str() { "mitochondrion" => "NUMT", "plastid" => "NUPT", _ => "partial_target_match" }}));
            }
            let best = match call {
                Call::NumtCandidate(t) => Some(&bait.targets[t].name),
                _ => None,
            };
            writeln!(
                rows,
                "{}",
                json!({"read_id": String::from_utf8_lossy(&rec.id), "read_length": rec.seq.len(),
                "best_bait_target": best, "cover_fraction": score.cover_frac, "hit_fraction": score.hit_frac,
                "targets": targets, "nuclear_origin_certified": false, "flank_origin": "unknown"})
            )?;
            Ok(())
        })();
    })?;
    failure?;
    let mut failure = Ok(());
    for writer in &mut writers {
        writer.flush()?;
    }
    rows.flush()?;
    // These are all accepted competing-target ties, before optional target-depth downsampling.
    let mut shared = BufWriter::new(File::create(dir.join("mt-pt-shared-reads.fastq"))?);
    let mt = bait.targets.iter().position(|t| t.name == "mitochondrion");
    let pt = bait.targets.iter().position(|t| t.name == "plastid");
    let mut shared_count = 0u64;
    if let (Some(mt), Some(pt)) = (mt, pt) {
        for_each_rec(&out.join(".pass.mitochondrion.fastq"), |rec| {
            if failure.is_err() {
                return;
            }
            if let Call::Organelle(mask) = classify(&rec.seq, bait, params).0 {
                let allowed = mask & !exclusions.map_or(0, |plan| plan.mask(&rec.id));
                if allowed & (1 << mt) != 0 && allowed & (1 << pt) != 0 {
                    failure = write_rec(&mut shared, &rec);
                    shared_count += 1;
                }
            }
        })?;
        failure?;
    }
    shared.flush()?;
    let mut report = json!({"scope": "partial bait matches and competing-target ties; no nuclear-reference or insertion-junction certification",
        "coordinates": "0-based half-open read intervals; gap-bridged bait-compatible blocks, not alignments",
        "original_merged_pool": pool, "original_pool_preserved": true, "classification_parameters": params,
        "short_read_limit": "reads shorter than numt_min_block cannot enter this candidate pool; an empty pool does not exclude NUMTs/NUPTs",
        "target_files": bait.targets.iter().enumerate().map(|(t,target)| json!({"target": target.name, "reads": counts[t],
            "fastq": dir.join(format!("{}.fastq",target.name))})).collect::<Vec<_>>(),
        "read_evidence": dir.join("nuclear-insertion-candidates.jsonl"),
        "mt_pt_shared_reads": shared_count, "mt_pt_shared_fastq": dir.join("mt-pt-shared-reads.fastq"),
        "transfer_direction": Value::Null, "elapsed_seconds": start.elapsed().as_secs_f64()});
    if let Some(audit) = frontier {
        report["frontier_evidence"] = audit.write(&dir, bait, params.block_gap)?;
    }
    save(&dir.join("summary.json"), &report)?;
    Ok(report)
}

#[derive(Clone, serde::Serialize)]
struct Position {
    segment: String,
    start: usize,
    end: usize,
    canonical_strand: &'static str,
}
fn index(graph: &Graph, k: usize) -> FxHashMap<u64, Vec<Position>> {
    let mut map: FxHashMap<u64, Vec<Position>> = FxHashMap::default();
    for (s, seq) in graph.seqs.iter().enumerate() {
        for (start, key) in CanonicalKmers::new(seq, k) {
            let word = &seq[start..start + k];
            let strand = if word <= revcomp(word).as_slice() {
                "+"
            } else {
                "-"
            };
            map.entry(key).or_default().push(Position {
                segment: graph.names[s].clone(),
                start,
                end: start + k,
                canonical_strand: strand,
            });
        }
    }
    map
}

fn support(path: &Path, keys: &FxHashMap<u64, usize>, k: usize) -> Result<Vec<u64>> {
    let mut counts = vec![0; keys.len()];
    for_each_rec(path, |rec| {
        let seen: FxHashSet<_> = CanonicalKmers::new(&rec.seq, k)
            .filter_map(|(_, key)| keys.get(&key).copied())
            .collect();
        for i in seen {
            counts[i] += 1;
        }
    })?;
    Ok(counts)
}

/// Exact shared words retain every graph occurrence, including repeats. This intentionally
/// makes no claim to exhaust divergent homology or identify transfer direction.
pub fn cross_organelle(assembly: &Path, out: &Path, k: usize) -> Result<Value> {
    let start = Instant::now();
    let dir = out.join("source_candidates");
    fs::create_dir_all(&dir)?;
    let mt_path = assembly.join("mitochondrion/organelle.gfa");
    let pt_path = assembly.join("plastid/organelle.gfa");
    if !mt_path.is_file() || !pt_path.is_file() {
        let report = json!({"status": "not_evaluated", "reason": "both organelle graphs required", "transfer_direction": Value::Null});
        save(&dir.join("mt-pt-sharing.json"), &report)?;
        return Ok(report);
    }
    let mt_graph = parse_gfa(&mt_path)?;
    let pt_graph = parse_gfa(&pt_path)?;
    let mt = index(&mt_graph, k);
    let pt = index(&pt_graph, k);
    let mut words: Vec<_> = mt
        .keys()
        .filter(|key| pt.contains_key(key))
        .copied()
        .collect();
    words.sort_unstable();
    let lookup: FxHashMap<_, _> = words.iter().enumerate().map(|(i, &key)| (key, i)).collect();
    let mut counts = BTreeMap::new();
    for target in ["mitochondrion", "plastid"] {
        let path = assembly.join("recruit").join(format!("{target}.fastq"));
        counts.insert(
            target,
            if path.is_file() {
                Some(support(&path, &lookup, k)?)
            } else {
                None
            },
        );
    }
    let mut diagonals: BTreeMap<(String, String, bool, i64), Vec<usize>> = BTreeMap::new();
    for key in &words {
        if let ([a], [b]) = (mt[key].as_slice(), pt[key].as_slice()) {
            let forward = a.canonical_strand == b.canonical_strand;
            let pt_len =
                pt_graph.seqs[pt_graph.names.iter().position(|n| n == &b.segment).unwrap()].len();
            let pos = if forward { b.start } else { pt_len - b.end };
            diagonals
                .entry((
                    a.segment.clone(),
                    b.segment.clone(),
                    forward,
                    pos as i64 - a.start as i64,
                ))
                .or_default()
                .push(a.start);
        }
    }
    let block_file = dir.join("mt-pt-unique-seed-blocks.jsonl");
    let mut block_rows = BufWriter::new(File::create(&block_file)?);
    let mut block_count = 0u64;
    for ((mt_segment, pt_segment, forward, diagonal), mut starts) in diagonals {
        starts.sort_unstable();
        let mut ranges: Vec<(usize, usize)> = Vec::new();
        for pos in starts {
            if let Some(last) = ranges.last_mut() {
                if pos <= last.1 {
                    last.1 = last.1.max(pos + k);
                    continue;
                }
            }
            ranges.push((pos, pos + k));
        }
        let pt_len = pt_graph.seqs[pt_graph
            .names
            .iter()
            .position(|n| n == &pt_segment)
            .unwrap()]
        .len() as i64;
        for (a, b) in ranges {
            let (c, d) = if forward {
                (a as i64 + diagonal, b as i64 + diagonal)
            } else {
                (pt_len - b as i64 - diagonal, pt_len - a as i64 - diagonal)
            };
            writeln!(
                block_rows,
                "{}",
                json!({"mitochondrion_segment": mt_segment,"mitochondrion_start": a,"mitochondrion_end": b,
                "plastid_segment": pt_segment,"plastid_start":c,"plastid_end":d,"relative_strand":if forward { "+" } else { "-" },
                "exact_length":b-a,"transfer_direction":Value::Null,"physical_origin_certified":false})
            )?;
            block_count += 1;
        }
    }
    block_rows.flush()?;
    let output = dir.join("mt-pt-shared-words.jsonl");
    let mut rows = BufWriter::new(File::create(&output)?);
    for (i, key) in words.iter().enumerate() {
        writeln!(
            rows,
            "{}",
            json!({"canonical_kmer_hex": format!("{key:016x}"), "mitochondrion": mt[key], "plastid": pt[key],
            "mitochondrion_recruited_reads": counts["mitochondrion"].as_ref().map(|c|c[i]),
            "plastid_recruited_reads": counts["plastid"].as_ref().map(|c|c[i]),
            "transfer_direction": Value::Null, "physical_origin_certified": false})
        )?;
    }
    rows.flush()?;
    let report = json!({"status": "evaluated", "scope": "exact shared k-mers on both emitted graphs, including every repeat occurrence; not exhaustive divergent homology",
        "k": k, "shared_words": words.len(), "coordinates": "0-based half-open GFA segment coordinates; relative strand follows canonical_strand agreement",
        "word_evidence": output, "unique_seed_blocks": block_file, "unique_seed_block_count":block_count,
        "block_scope":"collinear overlapping exact words unique in each graph; repeated and unresolved words remain in word_evidence, not forced into a unique block",
        "mitochondrion_graph": mt_path, "plastid_graph": pt_path,
        "read_support_scope": "distinct records containing each exact word, per recruited target file; overlapping target files must not be summed; not insertion-flank or repeat-phase evidence",
        "transfer_direction": Value::Null, "physical_origin_certified": false, "elapsed_seconds": start.elapsed().as_secs_f64()});
    save(&dir.join("mt-pt-sharing.json"), &report)?;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frontier_samples_are_bounded_and_match_exact_records() {
        let mut audit = FrontierAudit::new(1);
        let record = |id: &[u8], seq: &[u8]| Rec { id: id.to_vec(), seq: seq.to_vec(), qual: None };
        let reason = FrontierReason::ReadDepthProfile;
        for key in (0..20).rev() {
            let name = format!("read-{key}");
            audit.record(0, 1, reason, key, &record(name.as_bytes(), b"AAAA"));
        }
        audit.record(0, 2, reason, 0, &record(b"read-0", b"AAAA"));
        audit.record(0, 2, reason, 0, &record(b"read-0", b"CCCC"));
        audit.record(0, 2, reason, 0, &record(&vec![b'x'; FRONTIER_AUDIT_ID_BYTES + 1], b"A"));
        audit.count(0, reason, true, 3);
        audit.prepare_pool_membership();
        audit.match_pool_record(&record(b"read-0", b"CCCC"));
        audit.match_pool_record(&record(b"read-1", b"GGGG"));

        let events = &audit.events[0][reason as usize];
        assert_eq!(events.extension_records, 23);
        assert_eq!(events.within_pass_records, 3);
        assert_eq!(events.oversized_records, 1);
        assert_eq!(events.witnesses.len(), FRONTIER_AUDIT_EXAMPLES);
        let first: Vec<_> = events.witnesses.iter().take(3).map(|w| {
            (w.key, w.sequence.as_slice(), w.in_final_partial_pool)
        }).collect();
        assert_eq!(first, vec![(0, b"AAAA".as_slice(), false), (0, b"CCCC".as_slice(), true),
                              (1, b"AAAA".as_slice(), false)]);
    }

    #[test]
    fn shared_words_keep_repeat_coordinates_and_reverse_orientation() {
        let seq = b"ACGTTGCAAGGCTTACGATCTAGCGGA".to_vec();
        let graph = Graph {
            names: vec!["p1".into(), "p2".into()],
            seqs: vec![seq.clone(), revcomp(&seq)],
            links: vec![],
            overlaps: vec![],
        };
        let idx = index(&graph, 17);
        let key = CanonicalKmers::new(&seq, 17).next().unwrap().1;
        assert_eq!(idx[&key].len(), 2);
        assert_eq!(idx[&key][0].start, 0);
        assert_eq!(idx[&key][1].start, seq.len() - 17);
        assert_ne!(idx[&key][0].canonical_strand, idx[&key][1].canonical_strand);
    }
    #[test]
    fn nuclear_candidate_targets_are_retained_without_certifying_flanks() {
        let dir = std::env::temp_dir().join(format!("ovasm-retained-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let mut bait = BaitSet::new(5);
        let t = bait.target("plastid").unwrap();
        bait.add_sequence_counted(t, b"ACGTTGCAAGGCTTACGATCTAGCGGA");
        let seq = b"NNNNNNNNNNACGTTGCAAGGCTTACGATCTAGCGGANNNNNNNNNN";
        let pool = dir.join("numt_candidates.fastq");
        let data = format!(
            "@x\n{}\n+\n{}\n",
            String::from_utf8_lossy(seq),
            "#".repeat(seq.len())
        );
        fs::write(&pool, &data).unwrap();
        let params = ClassifyParams {
            min_distinct: 0.0,
            min_cover: 0.8,
            min_hit_frac: 0.8,
            numt_min_block: 15,
            block_gap: 0,
        };
        let report = recruitment(&dir, &pool, &bait, &params, None, None).unwrap();
        assert_eq!(report["target_files"][0]["reads"], 1);
        assert_eq!(
            fs::read_to_string(dir.join("source_candidates/plastid.fastq")).unwrap(),
            data
        );
        let row: Value = serde_json::from_str(
            fs::read_to_string(dir.join("source_candidates/nuclear-insertion-candidates.jsonl"))
                .unwrap()
                .trim(),
        )
        .unwrap();
        assert_eq!(row["targets"][0]["candidate_label"], "NUPT");
        assert_eq!(row["nuclear_origin_certified"], false);
        assert_eq!(fs::read_to_string(&pool).unwrap(), data);
        fs::remove_dir_all(dir).unwrap();
    }
    #[test]
    fn sharing_reports_coordinates_repeats_and_distinct_record_support() {
        let dir = std::env::temp_dir().join(format!("ovasm-sharing-{}", std::process::id()));
        for d in ["mitochondrion", "plastid", "recruit"] {
            fs::create_dir_all(dir.join(d)).unwrap();
        }
        let seq = b"ACGTTGCAAGGCTTACGATCTAGCGGA";
        fs::write(
            dir.join("mitochondrion/organelle.gfa"),
            format!("S\tm\t{}\n", String::from_utf8_lossy(seq)),
        )
        .unwrap();
        fs::write(
            dir.join("plastid/organelle.gfa"),
            format!("S\tp\t{}\n", String::from_utf8_lossy(&revcomp(seq))),
        )
        .unwrap();
        fs::write(
            dir.join("recruit/mitochondrion.fastq"),
            format!(
                "@mread\n{}\n+\n{}\n",
                String::from_utf8_lossy(seq),
                "!".repeat(seq.len())
            ),
        )
        .unwrap();
        fs::write(dir.join("recruit/plastid.fastq"), "").unwrap();
        let report = cross_organelle(&dir, &dir, 17).unwrap();
        assert_eq!(report["status"], "evaluated");
        let blocks: Vec<Value> =
            fs::read_to_string(dir.join("source_candidates/mt-pt-unique-seed-blocks.jsonl"))
                .unwrap()
                .lines()
                .map(|s| serde_json::from_str(s).unwrap())
                .collect();
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0]["relative_strand"], "-");
        assert_eq!(blocks[0]["mitochondrion_start"], 0);
        assert_eq!(blocks[0]["plastid_start"], 0);
        assert_eq!(blocks[0]["exact_length"], seq.len());
        assert!(blocks[0]["transfer_direction"].is_null());
        for row in fs::read_to_string(dir.join("source_candidates/mt-pt-shared-words.jsonl"))
            .unwrap()
            .lines()
        {
            let row: Value = serde_json::from_str(row).unwrap();
            assert_eq!(row["mitochondrion_recruited_reads"], 1);
            assert_eq!(row["plastid_recruited_reads"], 0);
        }
        // Empty evidence is distinct from missing input: missing reports null, never zero.
        fs::remove_file(dir.join("recruit/plastid.fastq")).unwrap();
        cross_organelle(&dir, &dir, 17).unwrap();
        let first =
            fs::read_to_string(dir.join("source_candidates/mt-pt-shared-words.jsonl")).unwrap();
        let row: Value = serde_json::from_str(first.lines().next().unwrap()).unwrap();
        assert!(row["plastid_recruited_reads"].is_null());
        fs::remove_dir_all(dir).unwrap();
    }
}
