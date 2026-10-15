//! Retained source candidates. Bait compatibility and exact homology are not source certification.
use crate::bait::BaitSet;
use crate::evidence::{parse_gfa, revcomp, Graph};
use crate::kmer::CanonicalKmers;
use crate::recruit::{classify, for_each_rec, write_rec, Call, ClassifyParams};
use anyhow::Result;
use rustc_hash::{FxHashMap, FxHashSet};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

fn save(path: &Path, value: &Value) -> Result<()> {
    fs::write(path, serde_json::to_string_pretty(value)? + "\n")?;
    Ok(())
}

/// Per-target compatible blocks use precisely the recruitment gap rule.
fn blocks(seq: &[u8], bait: &BaitSet, gap: usize) -> Vec<Vec<(usize, usize)>> {
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

pub fn recruitment(
    out: &Path,
    pool: &Path,
    bait: &BaitSet,
    params: &ClassifyParams,
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
    let mut failure = Ok(());
    for_each_rec(pool, |rec| {
        if failure.is_err() {
            return;
        }
        failure = (|| -> Result<()> {
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
                if mask & (1 << mt) != 0 && mask & (1 << pt) != 0 {
                    failure = write_rec(&mut shared, &rec);
                    shared_count += 1;
                }
            }
        })?;
        failure?;
    }
    shared.flush()?;
    let report = json!({"scope": "partial bait matches and competing-target ties; no nuclear-reference or insertion-junction certification",
        "coordinates": "0-based half-open read intervals; gap-bridged bait-compatible blocks, not alignments",
        "original_merged_pool": pool, "original_pool_preserved": true, "classification_parameters": params,
        "short_read_limit": "reads shorter than numt_min_block cannot enter this candidate pool; an empty pool does not exclude NUMTs/NUPTs",
        "target_files": bait.targets.iter().enumerate().map(|(t,target)| json!({"target": target.name, "reads": counts[t],
            "fastq": dir.join(format!("{}.fastq",target.name))})).collect::<Vec<_>>(),
        "read_evidence": dir.join("nuclear-insertion-candidates.jsonl"),
        "mt_pt_shared_reads": shared_count, "mt_pt_shared_fastq": dir.join("mt-pt-shared-reads.fastq"),
        "transfer_direction": Value::Null, "elapsed_seconds": start.elapsed().as_secs_f64()});
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
            min_cover: 0.8,
            min_hit_frac: 0.8,
            numt_min_block: 15,
            block_gap: 0,
        };
        let report = recruitment(&dir, &pool, &bait, &params).unwrap();
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
