//! Read-threaded repeat phases (step U1 of docs/unified-assembly-design.zh-CN.md section 6).
//!
//! Diagnostic only: this never changes a graph. For every repeat node of an assembly graph
//! (a node with at least two distinct entries or at least two distinct exits, decided from the
//! graph alone) it counts the complete original reads that walk entry -> node -> exit, one count
//! per read, and the reads that only touch one side. A phase (entry, exit) is supported when at
//! least `min_count` distinct reads show it, the same "distinct original reads" rule the
//! junction stage uses. The reads' alignment, competitive link assignment and read-order walks
//! are those of `evidence`; nothing here adds a cutoff of its own.
//!
//! Classes of a repeat node, from the supported phases only:
//! - `unresolved`: no phase reaches `min_count` (reads shorter than the repeat, or too few);
//! - `resolved`: the supported phases pair every entry with exactly one exit and vice versa;
//! - `partial`: supported phases pair entries with exits one to one but leave some end unpaired;
//! - `conflicted`: a one-to-one supported pairing has below-threshold counterexamples; it must
//!   not be treated as resolved simply because the counterexamples miss `min_count`;
//! - `mixed`: some entry (or exit) has two or more supported partners, as with the master and the
//!   recombinant configuration of a recombining repeat. Both are reported, none is chosen.
//!
//! Reads below `min_count` for a phase are reported as `reads_below_min_count`, not discarded.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::hash::{Hash, Hasher};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Result};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::Serialize;

use crate::evidence::{
    assign, build_index, crossings, flip, parse_gfa, read_steps, traversals, EvidenceParams, Graph,
    Index, Oriented,
};
use crate::recruit::{peak_rss_mb, read_batches};

/// Read ids kept per phase as witnesses.
const WITNESS_CAP: usize = 3;

type Triple = (Oriented, Oriented, Oriented);
type Adjacency = (Oriented, Oriented);

/// What one read says: the repeat traversals X -> S -> Y it completes (S forward) and the
/// directed adjacencies it crosses (both strand views).
pub(crate) struct ReadThread {
    pub triples: BTreeSet<Triple>,
    pub adjacencies: BTreeSet<Adjacency>,
    /// Link assignment was ambiguous for this read.
    pub ambiguous: bool,
    /// The read carries a non-ACGT base inside the span its crossings cover: its traversals
    /// are not counted as complete (its adjacencies still are).
    pub ambiguous_base: bool,
}

pub(crate) fn thread_read(g: &Graph, idx: &Index, p: &EvidenceParams, seq: &[u8]) -> ReadThread {
    let (kept, ambiguous) = assign(g, crossings(seq, idx, p), p.spacing_tol);
    let steps = read_steps(g, &kept);
    let mut adjacencies = BTreeSet::new();
    let (mut lo, mut hi) = (i64::MAX, i64::MIN);
    for &(a, b, start, end) in &steps {
        adjacencies.insert((a, b));
        adjacencies.insert((flip(b), flip(a)));
        lo = lo.min(start.min(end));
        hi = hi.max(start.max(end));
    }
    let ambiguous_base = !steps.is_empty() && {
        let (lo, hi) = (lo.max(0) as usize, (hi.max(0) as usize).min(seq.len()));
        seq[lo.min(hi)..hi]
            .iter()
            .any(|b| !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
    };
    let mut triples = BTreeSet::new();
    // A competing, unresolved anchor in this complete read prevents a unique phase call.
    // Its unambiguous adjacencies remain available as partial evidence.
    if !ambiguous_base && ambiguous.is_empty() {
        for (x, s, y) in traversals(g, &steps, p.spacing_tol) {
            triples.insert(if s.1 {
                (x, s, y)
            } else {
                (flip(y), flip(s), flip(x))
            });
        }
    }
    ReadThread {
        triples,
        adjacencies,
        ambiguous: !ambiguous.is_empty(),
        ambiguous_base,
    }
}

struct Ends {
    entries: Vec<Oriented>,
    exits: Vec<Oriented>,
}

/// Repeat nodes by graph structure alone: forward node `s` with at least two distinct entries
/// or at least two distinct exits (a link and its reverse complement are one adjacency).
fn repeat_nodes(g: &Graph) -> BTreeMap<usize, Ends> {
    let mut pred: BTreeMap<Oriented, BTreeSet<Oriented>> = BTreeMap::new();
    let mut succ: BTreeMap<Oriented, BTreeSet<Oriented>> = BTreeMap::new();
    for &(a, b) in &g.links {
        succ.entry(a).or_default().insert(b);
        pred.entry(b).or_default().insert(a);
        succ.entry(flip(b)).or_default().insert(flip(a));
        pred.entry(flip(a)).or_default().insert(flip(b));
    }
    let mut out = BTreeMap::new();
    for s in 0..g.seqs.len() {
        let f = (s, true);
        let entries: Vec<Oriented> = pred.get(&f).into_iter().flatten().copied().collect();
        let exits: Vec<Oriented> = succ.get(&f).into_iter().flatten().copied().collect();
        if entries.len() >= 2 || exits.len() >= 2 {
            out.insert(s, Ends { entries, exits });
        }
    }
    out
}

#[derive(Default)]
pub(crate) struct Tally {
    phases: BTreeMap<Triple, (u64, Vec<String>)>,
    adjacency: BTreeMap<Adjacency, u64>,
    pub reads_scanned: u64,
    pub unique_records: u64,
    pub reads_with_thread: u64,
    pub ambiguous_reads: u64,
    pub ambiguous_base_reads: u64,
}

impl Tally {
    pub(crate) fn add(&mut self, id: &str, t: &ReadThread) {
        self.unique_records += 1;
        self.ambiguous_reads += u64::from(t.ambiguous);
        self.ambiguous_base_reads += u64::from(t.ambiguous_base);
        self.reads_with_thread += u64::from(!t.triples.is_empty());
        for &a in &t.adjacencies {
            *self.adjacency.entry(a).or_default() += 1;
        }
        for &tr in &t.triples {
            let e = self.phases.entry(tr).or_default();
            e.0 += 1;
            if e.1.len() < WITNESS_CAP {
                e.1.push(id.to_string());
            }
        }
    }
}

#[derive(Serialize)]
pub struct EndReport {
    pub segment: String,
    /// Distinct reads crossing the link between this end and the node.
    pub adjacency_reads: u64,
}

#[derive(Serialize)]
pub struct PhaseReport {
    pub entry: String,
    pub exit: String,
    pub reads: u64,
    pub supported: bool,
    pub witnesses: Vec<String>,
}

#[derive(Serialize)]
pub struct RepeatReport {
    pub segment: String,
    pub length: usize,
    pub entries: Vec<EndReport>,
    pub exits: Vec<EndReport>,
    pub phases: Vec<PhaseReport>,
    pub classification: &'static str,
    pub reads_below_min_count: u64,
}

#[derive(Serialize)]
pub struct UnzipReport {
    pub method: &'static str,
    /// Structural identity, independent of GFA tags, line order and mirrored links.
    pub graph_fingerprint: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cleaning_support: Option<serde_json::Value>,
    pub min_count: u32,
    pub reads_scanned: u64,
    pub unique_records: u64,
    pub reads_with_complete_traversal: u64,
    pub ambiguous_reads: u64,
    pub reads_with_ambiguous_base_in_span: u64,
    pub repeat_nodes: Vec<RepeatReport>,
    pub graph_changed: bool,
    pub physical_molecule_count: Option<u32>,
    pub physical_topology_identified: bool,
    pub elapsed_seconds: f64,
    pub peak_rss_mb: Option<f64>,
}

/// Pinned rustc-hash fingerprint for binding evidence to sequence and overlap topology.
/// This is an identity check for accidental stale reports, not a cryptographic checksum.
pub(crate) fn graph_fingerprint(g: &Graph) -> String {
    let mut segments: Vec<_> = g.names.iter().zip(&g.seqs).collect();
    segments.sort();
    let mut links: Vec<_> = g
        .links
        .iter()
        .zip(&g.overlaps)
        .map(|(&(a, b), &ov)| {
            let forward = (g.label(a), g.label(b), ov);
            let reverse = (g.label(flip(b)), g.label(flip(a)), ov);
            forward.min(reverse)
        })
        .collect();
    links.sort();
    links.dedup();
    let mut h = rustc_hash::FxHasher::default();
    segments.hash(&mut h);
    links.hash(&mut h);
    format!("fxhash-v1:{:016x}", h.finish())
}

fn classify(
    entries: &[Oriented],
    exits: &[Oriented],
    supported: &BTreeSet<(Oriented, Oriented)>,
    has_counterexamples: bool,
) -> &'static str {
    if supported.is_empty() {
        return "unresolved";
    }
    let mut per_entry: BTreeMap<Oriented, usize> = BTreeMap::new();
    let mut per_exit: BTreeMap<Oriented, usize> = BTreeMap::new();
    for &(x, y) in supported {
        *per_entry.entry(x).or_default() += 1;
        *per_exit.entry(y).or_default() += 1;
    }
    if per_entry.values().chain(per_exit.values()).any(|&n| n >= 2) {
        "mixed"
    } else if has_counterexamples {
        "conflicted"
    } else if entries.iter().all(|e| per_entry.contains_key(e))
        && exits.iter().all(|x| per_exit.contains_key(x))
    {
        "resolved"
    } else {
        "partial"
    }
}

fn report_nodes(g: &Graph, tally: &Tally, min_count: u32) -> Vec<RepeatReport> {
    let min = u64::from(min_count);
    repeat_nodes(g)
        .into_iter()
        .map(|(s, ends)| {
            let node = (s, true);
            let adj = |a: Adjacency| tally.adjacency.get(&a).copied().unwrap_or(0);
            let mut phases: Vec<PhaseReport> = tally
                .phases
                .iter()
                .filter(|((_, mid, _), _)| *mid == node)
                .map(|(&(x, _, y), (n, w))| PhaseReport {
                    entry: g.label(x),
                    exit: g.label(y),
                    reads: *n,
                    supported: *n >= min,
                    witnesses: w.clone(),
                })
                .collect();
            phases.sort_by(|a, b| (&a.entry, &a.exit).cmp(&(&b.entry, &b.exit)));
            let supported: BTreeSet<(Oriented, Oriented)> = tally
                .phases
                .iter()
                .filter(|((_, mid, _), (n, _))| *mid == node && *n >= min)
                .map(|(&(x, _, y), _)| (x, y))
                .collect();
            RepeatReport {
                segment: g.names[s].clone(),
                length: g.seqs[s].len(),
                entries: ends
                    .entries
                    .iter()
                    .map(|&x| EndReport {
                        segment: g.label(x),
                        adjacency_reads: adj((x, node)),
                    })
                    .collect(),
                exits: ends
                    .exits
                    .iter()
                    .map(|&y| EndReport {
                        segment: g.label(y),
                        adjacency_reads: adj((node, y)),
                    })
                    .collect(),
                classification: classify(
                    &ends.entries,
                    &ends.exits,
                    &supported,
                    phases.iter().any(|p| !p.supported && p.reads > 0),
                ),
                reads_below_min_count: phases
                    .iter()
                    .filter(|p| !p.supported)
                    .map(|p| p.reads)
                    .sum(),
                phases,
            }
        })
        .collect()
}

/// Thread every original read through `g` and report the phases of its repeat nodes.
pub fn run(
    gfa: &Path,
    inputs: Vec<PathBuf>,
    p: EvidenceParams,
    min_count: u32,
    out_json: &Path,
) -> Result<UnzipReport> {
    run_audited(gfa, inputs, p, min_count, out_json, None)
}

pub fn run_audited(
    gfa: &Path,
    inputs: Vec<PathBuf>,
    p: EvidenceParams,
    min_count: u32,
    out_json: &Path,
    cleaning_audit: Option<&Path>,
) -> Result<UnzipReport> {
    if min_count == 0 || inputs.is_empty() {
        bail!("unzip needs reads and a positive --min-count");
    }
    if out_json.exists() {
        let out = fs::canonicalize(out_json)?;
        if out == fs::canonicalize(gfa)?
            || inputs
                .iter()
                .any(|p| fs::canonicalize(p).ok().as_ref() == Some(&out))
        {
            bail!("unzip report cannot overwrite the graph or input reads");
        }
    }
    let started = Instant::now();
    let g = parse_gfa(gfa)?;
    let mut audit: Option<serde_json::Value> = cleaning_audit
        .map(|path| -> Result<_> {
            let v: serde_json::Value = serde_json::from_str(&fs::read_to_string(path)?)?;
            if v["restored_graph_fingerprint"].as_str() != Some(graph_fingerprint(&g).as_str()) {
                bail!("cleaning audit is not bound to this restored graph");
            }
            Ok(v)
        })
        .transpose()?;
    let idx = build_index(&g, &p);
    let repeats = repeat_nodes(&g);
    let mut tally = Tally::default();
    // A read named twice counts once; the same name for different sequences is an error.
    // This identity ledger is O(unique read IDs), unlike the bounded per-phase witnesses.
    // Sequences and complete read threads are retained only for the current bounded batch.
    let mut seen: FxHashMap<Vec<u8>, (usize, u64)> = FxHashMap::default();
    let (rx, reader) = read_batches(&inputs);
    for batch in rx {
        tally.reads_scanned += batch.len() as u64;
        let threads: Vec<_> = batch
            .par_iter()
            .map(|r| {
                let id =
                    r.id.split(|c| c.is_ascii_whitespace())
                        .next()
                        .unwrap_or(&[])
                        .to_vec();
                let mut h = rustc_hash::FxHasher::default();
                h.write(&r.seq.to_ascii_uppercase());
                let mut thread = thread_read(&g, &idx, &p, &r.seq);
                thread
                    .triples
                    .retain(|(_, s, _)| repeats.contains_key(&s.0));
                (id, (r.seq.len(), h.finish()), thread)
            })
            .collect();
        for (id, fingerprint, thread) in threads {
            if id.is_empty() {
                bail!("unzip requires nonempty original read IDs");
            }
            if let Some(previous) = seen.insert(id.clone(), fingerprint) {
                if previous != fingerprint {
                    bail!("the same read ID names different complete records");
                }
                continue;
            }
            tally.add(&String::from_utf8_lossy(&id), &thread);
        }
    }
    reader.join().expect("reader thread panicked")?;
    if let Some(v) = audit.as_mut() {
        let reports = report_nodes(&g, &tally, min_count);
        let proposals = v["proposals"]
            .as_array_mut()
            .ok_or_else(|| anyhow::anyhow!("cleaning audit has no proposals"))?;
        for proposal in proposals {
            let mut rows = Vec::new();
            for name in proposal["retained_backbone_segments"]
                .as_array()
                .ok_or_else(|| anyhow::anyhow!("missing retained segments"))?
            {
                let name = name
                    .as_str()
                    .ok_or_else(|| anyhow::anyhow!("invalid retained segment"))?;
                if let Some(s) = g.names.iter().position(|n| n == name) {
                    let n = (s, true);
                    let entries:Vec<_>=tally.adjacency.iter().filter(|((_,b),_)|*b==n).map(|(&(a,_),&count)|serde_json::json!({"segment":g.label(a),"reads":count})).collect();
                    let exits:Vec<_>=tally.adjacency.iter().filter(|((a,_),_)|*a==n).map(|(&(_,b),&count)|serde_json::json!({"segment":g.label(b),"reads":count})).collect();
                    rows.push(serde_json::json!({"segment":name,"status":"measured","entries":entries,"exits":exits,
                        "repeat_phase":reports.iter().find(|r|r.segment==name)}));
                } else {
                    rows.push(serde_json::json!({"segment":name,"status":"not_present_after_junction_restoration"}));
                }
            }
            proposal["original_read_evidence"] = serde_json::json!(rows);
        }
    }
    let report = UnzipReport {
        method: "complete original-read traversals of structurally branching nodes; one count per read; no graph change",
        graph_fingerprint: graph_fingerprint(&g),
        cleaning_support: audit,
        min_count,
        reads_scanned: tally.reads_scanned,
        unique_records: tally.unique_records,
        reads_with_complete_traversal: tally.reads_with_thread,
        ambiguous_reads: tally.ambiguous_reads,
        reads_with_ambiguous_base_in_span: tally.ambiguous_base_reads,
        repeat_nodes: report_nodes(&g, &tally, min_count),
        graph_changed: false,
        physical_molecule_count: None,
        physical_topology_identified: false,
        elapsed_seconds: started.elapsed().as_secs_f64(),
        peak_rss_mb: peak_rss_mb(),
    };
    if let Some(parent) = out_json.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(out_json, serde_json::to_string_pretty(&report)?)?;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::revcomp;

    fn rng_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut s = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        (0..n)
            .map(|_| {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                b"ACGT"[(s % 4) as usize]
            })
            .collect()
    }

    fn params() -> EvidenceParams {
        EvidenceParams {
            short_read_adjacency: false,
            k: 21,
            flank: 500,
            min_anchor: 200,
            diag_tol: 20,
            min_side_hits: 5,
            spacing_tol: 0.05,
            bootstrap: 0,
            seed: 1,
            max_repeat_len: 20_000,
        }
    }

    /// Entries A, C -> repeat R -> exits B, D.
    fn fixture() -> (Graph, [Vec<u8>; 5]) {
        let [a, b, c, d, r] = [11, 12, 13, 14, 15].map(|i| rng_seq(2000, i));
        let g = Graph {
            names: ["A", "B", "C", "D", "R"].map(String::from).to_vec(),
            seqs: vec![a.clone(), b.clone(), c.clone(), d.clone(), r.clone()],
            links: vec![
                ((0, true), (4, true)),
                ((2, true), (4, true)),
                ((4, true), (1, true)),
                ((4, true), (3, true)),
            ],
            overlaps: vec![0; 4],
        };
        (g, [a, b, c, d, r])
    }

    /// Entry tail + repeat + exit head.
    fn through(entry: &[u8], r: &[u8], exit: &[u8]) -> Vec<u8> {
        [&entry[1400..], r, &exit[..600]].concat()
    }

    fn tally(g: &Graph, reads: &[Vec<u8>]) -> Tally {
        let p = params();
        let idx = build_index(g, &p);
        let mut t = Tally::default();
        for (i, r) in reads.iter().enumerate() {
            t.add(&format!("r{i}"), &thread_read(g, &idx, &p, r));
        }
        t
    }

    fn node<'a>(report: &'a [RepeatReport], name: &str) -> &'a RepeatReport {
        report.iter().find(|r| r.segment == name).unwrap()
    }

    fn phase(r: &RepeatReport, entry: &str, exit: &str) -> u64 {
        r.phases
            .iter()
            .find(|p| p.entry == entry && p.exit == exit)
            .map_or(0, |p| p.reads)
    }

    #[test]
    fn only_structurally_branching_nodes_are_repeat_nodes() {
        let (g, _) = fixture();
        let nodes = repeat_nodes(&g);
        assert_eq!(nodes.keys().copied().collect::<Vec<_>>(), vec![4]);
        assert_eq!(nodes[&4].entries, vec![(0, true), (2, true)]);
        assert_eq!(nodes[&4].exits, vec![(1, true), (3, true)]);
    }

    #[test]
    fn both_pairings_with_enough_reads_are_mixed_and_both_are_reported() {
        let (g, [a, b, c, d, r]) = fixture();
        let mut reads = Vec::new();
        for _ in 0..4 {
            reads.push(through(&a, &r, &b));
            reads.push(through(&c, &r, &d));
        }
        for _ in 0..3 {
            reads.push(through(&a, &r, &d));
        }
        let rep = report_nodes(&g, &tally(&g, &reads), 3);
        let r = node(&rep, "R");
        assert_eq!(r.classification, "mixed");
        assert_eq!(phase(r, "A+", "B+"), 4);
        assert_eq!(phase(r, "C+", "D+"), 4);
        assert_eq!(phase(r, "A+", "D+"), 3);
        assert_eq!(r.reads_below_min_count, 0);
    }

    #[test]
    fn a_below_threshold_counterexample_prevents_a_resolved_call() {
        let (g, [a, b, c, d, r]) = fixture();
        let mut reads = Vec::new();
        for _ in 0..4 {
            reads.push(through(&a, &r, &b));
            reads.push(through(&c, &r, &d));
        }
        reads.push(through(&a, &r, &d)); // one read of the other pairing: below min_count 3
        let rep = report_nodes(&g, &tally(&g, &reads), 3);
        let r = node(&rep, "R");
        assert_eq!(r.classification, "conflicted");
        assert_eq!(r.reads_below_min_count, 1);
        assert!(
            !r.phases
                .iter()
                .find(|p| p.entry == "A+" && p.exit == "D+")
                .unwrap()
                .supported
        );
    }

    #[test]
    fn a_pairing_that_leaves_an_end_unpaired_is_partial() {
        let (g, [a, b, _, _, r]) = fixture();
        let reads: Vec<_> = (0..4).map(|_| through(&a, &r, &b)).collect();
        let rep = report_nodes(&g, &tally(&g, &reads), 3);
        assert_eq!(node(&rep, "R").classification, "partial");
    }

    #[test]
    fn reads_that_do_not_span_the_repeat_leave_it_unresolved_but_count_adjacency() {
        let (g, [a, _, c, _, r]) = fixture();
        // each read enters the repeat and stops inside it: no exit is seen
        let reads: Vec<_> = (0..5)
            .map(|i| [&(if i % 2 == 0 { &a } else { &c })[1400..], &r[..900]].concat())
            .collect();
        let rep = report_nodes(&g, &tally(&g, &reads), 3);
        let r = node(&rep, "R");
        assert_eq!(r.classification, "unresolved");
        assert!(r.phases.is_empty());
        let adjacency = |s: &str| {
            r.entries
                .iter()
                .find(|e| e.segment == s)
                .unwrap()
                .adjacency_reads
        };
        assert_eq!((adjacency("A+"), adjacency("C+")), (3, 2));
    }

    #[test]
    fn strands_give_the_same_phases() {
        let (g, [a, b, ..]) = fixture();
        let p = params();
        let idx = build_index(&g, &p);
        let read = through(&a, &g.seqs[4], &b);
        let fwd = thread_read(&g, &idx, &p, &read);
        let rev = thread_read(&g, &idx, &p, &revcomp(&read));
        assert_eq!(fwd.triples, rev.triples);
        assert_eq!(fwd.triples.len(), 1);
        assert_eq!(fwd.adjacencies, rev.adjacencies);
    }

    #[test]
    fn shared_44_base_exits_and_48_base_entries_keep_distinct_phases() {
        let (mut g, [mut a, mut b, mut c, mut d, r]) = fixture();
        a[1952..].copy_from_slice(&r[..48]);
        c[1952..].copy_from_slice(&r[..48]);
        b[..44].copy_from_slice(&r[1956..]);
        d[..44].copy_from_slice(&r[1956..]);
        g.seqs = vec![a.clone(), b.clone(), c.clone(), d.clone(), r.clone()];
        g.overlaps = vec![48, 48, 44, 44];
        let spell = |x: &[u8], y: &[u8]| [&x[1400..], &r[48..], &y[44..600]].concat();
        let reads = vec![spell(&a, &b), spell(&c, &d), revcomp(&spell(&a, &d))];
        let rep = report_nodes(&g, &tally(&g, &reads), 1);
        let repeat = node(&rep, "R");
        assert_eq!(repeat.classification, "mixed");
        assert_eq!(phase(repeat, "A+", "B+"), 1);
        assert_eq!(phase(repeat, "C+", "D+"), 1);
        assert_eq!(phase(repeat, "A+", "D+"), 1);
        assert_eq!(phase(repeat, "C+", "B+"), 0);
    }

    #[test]
    fn an_unresolved_anchor_in_a_read_prevents_a_complete_phase_call() {
        let (mut g, [a, b, _, _, r]) = fixture();
        let x = rng_seq(2000, 71);
        let y = rng_seq(2000, 72);
        g.names.extend(["X", "Y", "Z"].map(String::from));
        g.seqs.extend([x.clone(), y.clone(), y.clone()]);
        g.links
            .extend([((5, true), (6, true)), ((5, true), (7, true))]);
        g.overlaps.extend([0, 0]);
        let read = [through(&a, &r, &b), x[1400..].to_vec(), y[..600].to_vec()].concat();
        let p = params();
        let t = thread_read(&g, &build_index(&g, &p), &p, &read);
        assert!(t.ambiguous);
        assert!(t.triples.is_empty());
        assert!(!t.adjacencies.is_empty());
    }

    #[test]
    fn original_reads_materialize_resolved_copies_with_shared_overlaps_on_both_strands() {
        let (mut g, [mut a, mut b, mut c, mut d, r]) = fixture();
        a[1952..].copy_from_slice(&r[..48]);
        c[1952..].copy_from_slice(&r[..48]);
        b[..44].copy_from_slice(&r[1956..]);
        d[..44].copy_from_slice(&r[1956..]);
        g.seqs = vec![a.clone(), b.clone(), c.clone(), d.clone(), r.clone()];
        g.overlaps = vec![48, 48, 44, 44];
        // Store one physical link in its reverse-complement representation.
        g.links[0] = (flip(g.links[0].1), flip(g.links[0].0));
        let dir =
            std::env::temp_dir().join(format!("ovasm-unzip-materialize-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let mut text = String::new();
        for (name, seq) in g.names.iter().zip(&g.seqs) {
            text.push_str(&format!(
                "S\t{name}\t{}\tdp:f:100\n",
                String::from_utf8_lossy(seq)
            ));
        }
        for (&(x, y), ov) in g.links.iter().zip(&g.overlaps) {
            text.push_str(&format!(
                "L\t{}\t{}\t{}\t{}\t{ov}M\n",
                g.names[x.0],
                if x.1 { '+' } else { '-' },
                g.names[y.0],
                if y.1 { '+' } else { '-' }
            ));
        }
        let graph = dir.join("in.gfa");
        fs::write(&graph, &text).unwrap();
        let mut reads = String::new();
        for i in 0..80 {
            for (j, (x, y)) in [(&a, &b), (&c, &d)].iter().enumerate() {
                let forward = [&x[1400..], &r[48..], &y[44..600]].concat();
                let sequence = if i % 2 == 0 {
                    forward
                } else {
                    revcomp(&forward)
                };
                reads.push_str(&format!(
                    ">r{i}_{j}\n{}\n",
                    String::from_utf8_lossy(&sequence)
                ));
            }
        }
        let input = dir.join("reads.fa");
        fs::write(&input, reads).unwrap();
        let evidence = dir.join("phases.json");
        let report = run(&graph, vec![input], params(), 80, &evidence).unwrap();
        assert_eq!(report.unique_records, 160);
        assert_eq!(node(&report.repeat_nodes, "R").classification, "resolved");
        let out = dir.join("split.gfa");
        assert!(
            crate::phase::export(&graph, &evidence, None, Some(&out))
                .unwrap()
                .graph_changed
        );
        let split = parse_gfa(&out).unwrap();
        let links: BTreeSet<_> = split
            .links
            .iter()
            .flat_map(|&(a, b)| [(a, b), (flip(b), flip(a))])
            .collect();
        for (entry, exit) in [(0, 1), (2, 3)] {
            let first = (
                split
                    .names
                    .iter()
                    .position(|n| n == &g.names[entry])
                    .unwrap(),
                true,
            );
            let middle = links.iter().find(|(a, _)| *a == first).unwrap().1;
            let last = links.iter().find(|(a, _)| *a == middle).unwrap().1;
            assert_eq!(split.names[last.0], g.names[exit]);
            assert_eq!(
                crate::linearize::spell(&split, &[first, middle, last], false).unwrap(),
                crate::linearize::spell(&g, &[(entry, true), (4, true), (exit, true)], false)
                    .unwrap()
            );
        }
        assert_eq!(fs::read_to_string(&graph).unwrap(), text);
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn separated_junction_fragments_do_not_make_a_complete_traversal() {
        let (g, [a, b, _, _, r]) = fixture();
        // Both links are real but their placement in this artificial read cannot spell R.
        let read = [&a[1400..], &r[..600], &r[1400..], &b[..600]].concat();
        let p = params();
        let t = thread_read(&g, &build_index(&g, &p), &p, &read);
        assert_eq!(t.adjacencies.len(), 4);
        assert!(t.triples.is_empty());
    }

    #[test]
    fn an_ambiguous_base_inside_the_spanned_region_is_not_a_complete_read() {
        let (g, [a, b, _, _, r]) = fixture();
        let p = params();
        let idx = build_index(&g, &p);
        let mut read = through(&a, &r, &b);
        let middle = read.len() / 2;
        read[middle] = b'N';
        let t = thread_read(&g, &idx, &p, &read);
        assert!(t.ambiguous_base && t.triples.is_empty());
        assert!(
            !t.adjacencies.is_empty(),
            "its adjacencies are still evidence"
        );
        // outside the span the N does not matter
        let mut tail = through(&a, &r, &b);
        tail.extend_from_slice(b"NNNNNNNN");
        assert_eq!(thread_read(&g, &idx, &p, &tail).triples.len(), 1);
    }

    #[test]
    fn one_read_counts_once_however_many_traversals_it_repeats() {
        let (g, [a, b, ..]) = fixture();
        let r = g.seqs[4].clone();
        let one = through(&a, &r, &b);
        let t = tally(&g, &[one.clone(), one.clone(), revcomp(&one)]);
        // three records, three reads: count 3 (ids differ); within one record, once
        assert_eq!(t.phases.values().map(|v| v.0).collect::<Vec<_>>(), vec![3]);
        let p = params();
        let idx = build_index(&g, &p);
        let doubled = [one.clone(), one].concat();
        assert!(thread_read(&g, &idx, &p, &doubled).triples.len() <= 1);
    }

    #[test]
    fn witnesses_are_capped() {
        let (g, [a, b, ..]) = fixture();
        let r = g.seqs[4].clone();
        let reads: Vec<_> = (0..10).map(|_| through(&a, &r, &b)).collect();
        let t = tally(&g, &reads);
        let (n, w) = t.phases.values().next().unwrap();
        assert_eq!((*n, w.len()), (10, WITNESS_CAP));
    }

    #[test]
    fn a_graph_without_branching_has_no_repeat_nodes() {
        let (mut g, _) = fixture();
        g.links = vec![((0, true), (4, true)), ((4, true), (1, true))];
        g.overlaps = vec![0; 2];
        assert!(repeat_nodes(&g).is_empty());
    }

    #[test]
    fn the_run_writes_a_report_and_leaves_the_graph_untouched() {
        let (g, [a, b, c, d, r]) = fixture();
        let dir = std::env::temp_dir().join(format!("ovasm-unzip-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let mut gfa = String::new();
        for (n, s) in g.names.iter().zip(&g.seqs) {
            gfa.push_str(&format!("S\t{n}\t{}\n", String::from_utf8_lossy(s)));
        }
        for &(x, y) in &g.links {
            gfa.push_str(&format!(
                "L\t{}\t+\t{}\t+\t0M\n",
                g.names[x.0], g.names[y.0]
            ));
        }
        let gfa_path = dir.join("g.gfa");
        fs::write(&gfa_path, &gfa).unwrap();
        let mut fa = String::new();
        for i in 0..4 {
            for (j, read) in [through(&a, &r, &b), through(&c, &r, &d)]
                .iter()
                .enumerate()
            {
                fa.push_str(&format!(">r{i}_{j}\n{}\n", String::from_utf8_lossy(read)));
            }
        }
        // the same record twice, and a different sequence under a new name
        fa.push_str(&format!(
            ">r0_0\n{}\n",
            String::from_utf8_lossy(&through(&a, &r, &b))
        ));
        let reads = dir.join("reads.fa");
        fs::write(&reads, &fa).unwrap();
        let report = run(&gfa_path, vec![reads], params(), 3, &dir.join("unzip.json")).unwrap();
        assert_eq!((report.reads_scanned, report.unique_records), (9, 8));
        assert_eq!(node(&report.repeat_nodes, "R").classification, "resolved");
        assert!(!report.graph_changed && report.physical_molecule_count.is_none());
        assert_eq!(fs::read_to_string(&gfa_path).unwrap(), gfa);
        assert!(dir.join("unzip.json").exists());
        fs::remove_dir_all(dir).unwrap();
    }
}
