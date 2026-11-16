//! `ovasm pan cohort`: one population summary from many `pan reconstruct` runs on one panel.
//!
//! Each sample's run holds its organelle graph against the panel graph: the panel segments its
//! reads support, its own new segments from the de novo assembly, its molecules as walks of
//! de novo unitigs (the representative configuration, `configuration=1`, where the panel layer
//! lists several), and each repeat junction with what decided it (`junctions.tsv`: `source` =
//! reads, contradicted_by_reads, panel or undetermined). The cohort merges them per organelle
//! into:
//!
//! - `cohort.gfa`: the panel segments (their own names) and the new segments (`n<i>`, one per
//!   distinct sequence in either orientation), the panel's links and the links the samples'
//!   read-decided walks take, and one path per read-decided piece of each sample's molecules.
//!   Joins only the panel decided are cut, not merged, so the panel's own inference is not fed
//!   back as data.
//! - `presence.tsv` / `presence.calls.tsv`: segments x samples, depth (or present/absent).
//! - `genotypes.tsv` / `genotypes.vcf`: the panel's variant sites x samples (major allele).
//! - `junctions.tsv`: each repeat junction's configurations across the cohort, counted apart by
//!   what decided them: reads (supported or not) or the panel.
//! - `pairings.tsv`: every pairing the graph allowed at each anchor end of each sample's run
//!   (its `pairings.tsv`), one row per pairing and sample: fragments crossing it and the
//!   library's expectation, so a repeat's rearrangement can be measured across the population.
//! - `nearest.tsv`: each sample's nearest panel runs and nearest cohort sample by major alleles,
//!   and its run's own closest panel members by posterior weight (reconstruction.json).
//!
//! Presence calls are present or absent only where presence.tsv says so and the segment has at
//! least MIN_OWN_KMERS own k-mers; otherwise undetermined. Samples whose median read length is
//! below `min_read_length` (or that are data_limited) keep their own columns and paths but are
//! left out of every count and frequency (GFA NS/RS, VCF AF/NS, junctions.tsv).
//!
//! The panel graph is the one each run's `reference.json` names, or `--reference`; either way
//! every run's `presence.tsv` must list the panel's segments at their lengths.
use std::collections::BTreeMap;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};
use rustc_hash::{FxHashMap, FxHashSet};

/// A panel segment with fewer own k-mers than this cannot be called present or absent.
pub const MIN_OWN_KMERS: usize = 20;
/// Default shortest median read length a sample needs to count in the cohort's frequencies.
pub const MIN_READ_LENGTH: f64 = 90.0;

pub const ORGANELLES: [&str; 2] = ["mitochondrion", "plastid"];

/// One oriented step of a walk in the cohort graph.
#[derive(Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Step {
    pub seg: String,
    pub fwd: bool,
}

impl Step {
    fn flip(&self) -> Step {
        Step { seg: self.seg.clone(), fwd: !self.fwd }
    }
    fn text(&self) -> String {
        format!("{}{}", self.seg, if self.fwd { '+' } else { '-' })
    }
}

fn parse_oriented(s: &str) -> Result<(String, bool)> {
    let s = s.trim();
    let (name, o) = s.split_at(s.len().checked_sub(1).ok_or_else(|| anyhow!("empty step"))?);
    match o {
        "+" => Ok((name.to_string(), true)),
        "-" => Ok((name.to_string(), false)),
        _ => bail!("step without orientation: {s:?}"),
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            x => x,
        })
        .collect()
}

/// A tab-separated table read by column name.
pub struct Table {
    pub cols: FxHashMap<String, usize>,
    pub rows: Vec<Vec<String>>,
}

impl Table {
    pub fn read(path: &Path) -> Result<Table> {
        let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        Table::parse(&text).with_context(|| format!("parsing {}", path.display()))
    }

    pub fn parse(text: &str) -> Result<Table> {
        let mut lines = text.lines().filter(|l| !l.is_empty());
        let header = lines.next().ok_or_else(|| anyhow!("empty table"))?;
        let cols = header.split('\t').enumerate().map(|(i, c)| (c.to_string(), i)).collect();
        let rows = lines.map(|l| l.split('\t').map(str::to_string).collect()).collect();
        Ok(Table { cols, rows })
    }

    pub fn col(&self, name: &str) -> Result<usize> {
        self.cols.get(name).copied().ok_or_else(|| anyhow!("no column {name:?}"))
    }

    pub fn get<'a>(&self, row: &'a [String], name: &str) -> Result<&'a str> {
        let i = self.col(name)?;
        row.get(i).map(String::as_str).ok_or_else(|| anyhow!("short row for column {name:?}"))
    }
}

/// The panel graph: segment names in file order, sequences, links and paths.
pub struct Panel {
    pub names: Vec<String>,
    pub seqs: FxHashMap<String, Vec<u8>>,
    pub links: FxHashSet<(Step, Step)>,
    pub paths: BTreeMap<String, Vec<Step>>,
}

impl Panel {
    pub fn read(path: &Path) -> Result<Panel> {
        let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
        Panel::parse(&text)
    }

    pub fn parse(text: &str) -> Result<Panel> {
        let mut p = Panel { names: Vec::new(), seqs: FxHashMap::default(), links: FxHashSet::default(), paths: BTreeMap::new() };
        for line in text.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            match f.first() {
                Some(&"S") if f.len() >= 3 => {
                    p.names.push(f[1].to_string());
                    p.seqs.insert(f[1].to_string(), f[2].as_bytes().to_vec());
                }
                Some(&"L") if f.len() >= 5 => {
                    let a = Step { seg: f[1].to_string(), fwd: f[2] == "+" };
                    let b = Step { seg: f[3].to_string(), fwd: f[4] == "+" };
                    p.links.insert(canonical_link(&a, &b));
                }
                Some(&"P") if f.len() >= 3 => {
                    let steps = f[2]
                        .split(',')
                        .map(|s| parse_oriented(s).map(|(seg, fwd)| Step { seg, fwd }))
                        .collect::<Result<Vec<_>>>()?;
                    p.paths.insert(f[1].to_string(), steps);
                }
                _ => {}
            }
        }
        if p.names.is_empty() {
            bail!("panel graph has no segments");
        }
        Ok(p)
    }
}

/// A link and its reverse are one link: the lesser of the two spellings.
pub fn canonical_link(a: &Step, b: &Step) -> (Step, Step) {
    let fwd = (a.clone(), b.clone());
    let rev = (b.flip(), a.flip());
    if (fwd.0.text(), fwd.1.text()) <= (rev.0.text(), rev.1.text()) { fwd } else { rev }
}

/// The cohort's new segments: one per distinct sequence in either orientation.
#[derive(Default)]
pub struct NewSegments {
    pub ids: FxHashMap<Vec<u8>, String>,
    pub seqs: Vec<(String, Vec<u8>)>,
    pub samples: FxHashMap<String, FxHashSet<String>>,
}

impl NewSegments {
    /// The step for `seq` read forwards, adding its segment the first time.
    pub fn step(&mut self, seq: &[u8], sample: &str) -> Step {
        let fwd = seq.to_ascii_uppercase();
        let rev = revcomp(&fwd);
        let (canon, is_fwd) = if fwd <= rev { (fwd, true) } else { (rev, false) };
        let id = match self.ids.get(&canon) {
            Some(id) => id.clone(),
            None => {
                let id = format!("n{}", self.seqs.len() + 1);
                self.ids.insert(canon.clone(), id.clone());
                self.seqs.push((id.clone(), canon));
                id
            }
        };
        self.samples.entry(id.clone()).or_default().insert(sample.to_string());
        Step { seg: id, fwd: is_fwd }
    }
}

/// What one individual-graph piece stands for in the cohort graph.
#[derive(Clone, Debug)]
enum Piece {
    /// A stretch of a panel segment, in the piece's forward orientation.
    Panel(Step),
    /// Sequence of the sample's own (forward orientation).
    New(Vec<u8>),
}

/// One sample's run for one organelle, reduced to cohort terms.
pub struct SampleOrganelle {
    /// Read-decided pieces of its molecules as cohort walks.
    pub walks: Vec<(String, Vec<Step>)>,
    /// Panel segment -> (depth, call), from presence.tsv; the call is present, absent or
    /// undetermined (presence.tsv's low or no_own_kmers, or fewer than MIN_OWN_KMERS own k-mers).
    pub presence: FxHashMap<String, (Option<f64>, String)>,
    /// The run's closest panel members (run, posterior weight), from reconstruction.json.
    pub relatives: Vec<(String, f64)>,
    /// New segments it holds.
    pub new_segments: FxHashSet<String>,
    /// site -> (major allele, minor alleles), from genotypes.tsv.
    pub calls: BTreeMap<u64, (Option<String>, Vec<String>)>,
    /// Junctions: (key, source, class).
    pub junctions: Vec<(String, String, String)>,
    /// Pairings the graph allowed (key -> the run's pairings.tsv columns from repeat_bp on).
    pub pairings: BTreeMap<String, String>,
    /// Pieces its segment table names that its individual graph holds no sequence for.
    pub missing_pieces: usize,
    /// De novo unitigs the individual graph does not map (each taken as one new segment).
    pub unmapped_unitigs: usize,
    /// Whether the run wrote junctions.tsv.
    pub has_junctions: bool,
}

/// The first `key` anywhere under `v`, depth first.
fn find_key<'a>(v: &'a serde_json::Value, key: &str) -> Option<&'a serde_json::Value> {
    match v {
        serde_json::Value::Object(m) => m.get(key).or_else(|| m.values().find_map(|x| find_key(x, key))),
        serde_json::Value::Array(a) => a.iter().find_map(|x| find_key(x, key)),
        _ => None,
    }
}

/// The run's closest panel members for the organelle, by posterior weight (reconstruction.json,
/// `closest_relatives`: [{run, weight}]); empty when the run lists none.
fn closest_relatives(dir: &Path, o: &str) -> Result<Vec<(String, f64)>> {
    let p = dir.join("reconstruction.json");
    if !p.exists() {
        return Ok(Vec::new());
    }
    let v: serde_json::Value = serde_json::from_str(&fs::read_to_string(&p)?).with_context(|| format!("parsing {}", p.display()))?;
    let Some(list) = find_key(&v["organelles"][o], "closest_relatives").and_then(|x| x.as_array()) else { return Ok(Vec::new()) };
    Ok(list.iter().filter_map(|e| Some((e["run"].as_str()?.to_string(), e["weight"].as_f64()?))).collect())
}

/// Reads per-sample read lengths (a table with `sample`, `read_length_median` and, if present,
/// `data_limited`, such as a run directory's summary.tsv): sample -> (median, data_limited).
pub fn read_sample_info(path: &Path) -> Result<FxHashMap<String, (Option<f64>, bool)>> {
    let t = Table::read(path)?;
    let limited = t.col("data_limited").ok();
    let mut out = FxHashMap::default();
    for row in &t.rows {
        let len = t.get(row, "read_length_median")?.parse::<f64>().ok();
        let lim = limited.is_some_and(|i| row.get(i).is_some_and(|x| x == "1" || x == "true"));
        out.insert(t.get(row, "sample")?.to_string(), (len, lim));
    }
    Ok(out)
}

/// Whether a junction's join was decided by reads.
fn reads_decided(source: &str) -> bool {
    source == "reads"
}

/// Cuts a molecule's unitig walk where a junction the reads did not decide joins its repeat to
/// its right anchor: the walk keeps left anchor and repeat, and a new walk starts at the right
/// anchor. `junctions`: (left, repeat steps, right) in the walk's unitig names.
pub fn cut_walk(walk: &[String], undecided: &[(String, Vec<String>, String)]) -> Vec<Vec<String>> {
    let mut cuts: Vec<usize> = Vec::new();
    for (left, repeat, right) in undecided {
        let mut pattern = vec![left.clone()];
        pattern.extend(repeat.iter().cloned());
        pattern.push(right.clone());
        let k = pattern.len();
        if walk.len() < k {
            continue;
        }
        for i in 0..=walk.len() - k {
            if walk[i..i + k] == pattern[..] {
                cuts.push(i + k - 1);
            }
        }
        // the same join read backwards
        let rev: Vec<String> = pattern.iter().rev().map(|s| flip_name(s)).collect();
        for i in 0..=walk.len() - k {
            if walk[i..i + k] == rev[..] {
                cuts.push(i + 1);
            }
        }
    }
    cuts.sort_unstable();
    cuts.dedup();
    let mut out = Vec::new();
    let mut start = 0;
    for c in cuts {
        if c > start && c < walk.len() {
            out.push(walk[start..c].to_vec());
            start = c;
        }
    }
    out.push(walk[start..].to_vec());
    out.into_iter().filter(|w| !w.is_empty()).collect()
}

fn flip_name(s: &str) -> String {
    match s.as_bytes().last() {
        Some(b'+') => format!("{}-", &s[..s.len() - 1]),
        Some(b'-') => format!("{}+", &s[..s.len() - 1]),
        _ => s.to_string(),
    }
}

/// Merges consecutive steps on the same panel segment and orientation (the pieces a sample's
/// graph split it into) into one step.
pub fn merge_steps(steps: Vec<Step>) -> Vec<Step> {
    let mut out: Vec<Step> = Vec::new();
    for s in steps {
        if out.last() == Some(&s) && !s.seg.starts_with('n') {
            continue;
        }
        out.push(s);
    }
    out
}

/// Reads one sample's run for organelle `o` and puts it in cohort terms.
pub fn read_sample(
    sample: &str,
    dir: &Path,
    o: &str,
    panel: &Panel,
    new: &mut NewSegments,
) -> Result<Option<SampleOrganelle>> {
    let od = dir.join(o);
    if !od.join("individual.segments.tsv").exists() {
        return Ok(None);
    }
    // individual graph piece sequences
    let gfa = fs::read_to_string(od.join("individual.gfa"))
        .with_context(|| format!("{sample}: reading {o}/individual.gfa"))?;
    let mut node_seq: FxHashMap<String, Vec<u8>> = FxHashMap::default();
    for line in gfa.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.first() == Some(&"S") && f.len() >= 3 {
            node_seq.insert(f[1].to_string(), f[2].as_bytes().to_vec());
        }
    }
    // individual.segments.tsv names the pieces of the merge before parallel nodes were
    // coalesced; individual.parallel-nodes.tsv takes each to the node individual.gfa keeps
    // (pan_node -> individual_node, reversed)
    let mut piece_seq: FxHashMap<String, Vec<u8>> = FxHashMap::default();
    let parallel = od.join("individual.parallel-nodes.tsv");
    if parallel.exists() {
        let t = Table::read(&parallel)?;
        for row in &t.rows {
            let (pan, ind) = (t.get(row, "pan_node")?, t.get(row, "individual_node")?);
            if let Some(seq) = node_seq.get(ind) {
                let seq = if t.get(row, "reversed")? == "true" { revcomp(seq) } else { seq.clone() };
                piece_seq.insert(pan.to_string(), seq);
            }
        }
    }
    for (id, seq) in node_seq {
        piece_seq.entry(id).or_insert(seq);
    }
    // pieces of panel segments, and each de novo unitig's pieces
    let segs = Table::read(&od.join("individual.segments.tsv"))?;
    let mut piece: FxHashMap<String, Piece> = FxHashMap::default();
    let mut unitig_pieces: FxHashMap<String, Vec<(String, bool)>> = FxHashMap::default();
    for row in &segs.rows {
        let seg = segs.get(row, "segment")?;
        let pieces = segs
            .get(row, "pieces")?
            .split(',')
            .filter(|s| !s.is_empty())
            .map(parse_oriented)
            .collect::<Result<Vec<_>>>()?;
        if let Some(rest) = seg.strip_prefix("ref.") {
            let n: usize = rest.split('.').next().unwrap_or("").parse().with_context(|| format!("{sample}: segment {seg}"))?;
            let name = panel.names.get(n).ok_or_else(|| anyhow!("{sample}: {seg} beyond the panel's {} segments", panel.names.len()))?;
            for (p, fwd) in pieces {
                piece.insert(p, Piece::Panel(Step { seg: name.clone(), fwd }));
            }
        } else if let Some((_, unitig)) = seg.rsplit_once('.') {
            unitig_pieces.insert(unitig.to_string(), pieces);
        }
    }
    // De novo unitigs the individual graph does not map: the molecules name assembly.gfa
    // unitigs, and the individual graph is built from organelle.gfa, which can leave some out.
    // Each becomes one piece of its own sequence (a new segment).
    let mut unmapped_unitigs = 0usize;
    if let Ok(text) = fs::read_to_string(od.join("de-novo").join("assembly.gfa")) {
        for line in text.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.first() == Some(&"S") && f.len() >= 3 && !unitig_pieces.contains_key(f[1]) {
                let id = format!("unitig:{}", f[1]);
                piece_seq.insert(id.clone(), f[2].as_bytes().to_vec());
                unitig_pieces.insert(f[1].to_string(), vec![(id, true)]);
                unmapped_unitigs += 1;
            }
        }
    }
    let mut new_segments = FxHashSet::default();
    let mut missing: FxHashSet<String> = FxHashSet::default();
    // a piece without sequence in the individual graph is left out of the walk (counted)
    let mut step_of = |p: &str, fwd: bool| -> Option<Step> {
        let s = match piece.get(p) {
            Some(Piece::Panel(s)) => s.clone(),
            Some(Piece::New(_)) | None => {
                let Some(seq) = piece_seq.get(p) else {
                    missing.insert(p.to_string());
                    return None;
                };
                let s = new.step(seq, sample);
                new_segments.insert(s.seg.clone());
                piece.entry(p.to_string()).or_insert_with(|| Piece::New(seq.clone()));
                s
            }
        };
        Some(if fwd { s } else { s.flip() })
    };
    // junctions the reads did not decide, per molecule
    // a run whose molecules never reached the junction check writes no junctions.tsv: no
    // junctions, and its walks are not cut
    let junctions_file = od.join("junctions.tsv");
    let has_junctions = junctions_file.exists();
    let jt = if has_junctions { Table::read(&junctions_file)? } else { Table::parse("molecule\n")? };
    let mut undecided: FxHashMap<String, Vec<(String, Vec<String>, String)>> = FxHashMap::default();
    let mut junction_rows = Vec::new();
    for row in &jt.rows {
        let molecule = jt.get(row, "molecule")?.to_string();
        let molecule_of_row = molecule.clone();
        let left = jt.get(row, "left")?.to_string();
        let repeat: Vec<String> = jt.get(row, "repeat")?.split(',').map(str::to_string).collect();
        let right = jt.get(row, "right")?.to_string();
        let source = jt.get(row, "source").unwrap_or("reads").to_string();
        let class = jt.get(row, "class")?.to_string();
        if !reads_decided(&source) {
            undecided.entry(molecule).or_default().push((left.clone(), repeat.clone(), right.clone()));
        }
        junction_rows.push((molecule_of_row, left, repeat, right, source, class));
    }
    // molecules as de novo unitig walks: the reads layer (<o>.reads.fasta, cut where the reads
    // did not decide a join) where every one of its molecules names its walk, else the final
    // layer, cut here at the same joins
    let reads_layer = dir.join(format!("{o}.reads.fasta"));
    let fasta = match fs::read_to_string(&reads_layer) {
        Ok(text) if text.lines().any(|l| l.starts_with('>'))
            && text.lines().filter(|l| l.starts_with('>')).all(|l| l.split_whitespace().any(|w| w.starts_with("path="))) => text,
        // a run that kept no molecule of this organelle writes no <o>.fasta: no paths, but its
        // presence and genotypes still count
        _ => {
            let final_layer = dir.join(format!("{o}.fasta"));
            if final_layer.exists() {
                fs::read_to_string(&final_layer).with_context(|| format!("{sample}: reading {o}.fasta"))?
            } else {
                String::new()
            }
        }
    };
    let mut walks = Vec::new();
    let mut other_configuration: FxHashSet<String> = FxHashSet::default();
    for header in fasta.lines().filter(|l| l.starts_with('>')) {
        let name = header[1..].split_whitespace().next().unwrap_or("").to_string();
        // the panel layer can list several configurations (configuration=i of=n): the cohort
        // takes the representative one, so a sample counts once
        if header.split_whitespace().any(|w| w.strip_prefix("configuration=").is_some_and(|c| c != "1")) {
            // a piece names the molecule it was cut from, which junctions.tsv names
            if let Some(from) = header.split_whitespace().find_map(|w| w.strip_prefix("cut_from=")) {
                other_configuration.insert(from.to_string());
            }
            other_configuration.insert(name);
            continue;
        }
        let Some(path) = header.split_whitespace().find_map(|w| w.strip_prefix("path=")) else {
            continue;
        };
        let walk: Vec<String> = path.split(',').map(str::to_string).collect();
        let none = Vec::new();
        let from = header.split_whitespace().find_map(|w| w.strip_prefix("cut_from=")).unwrap_or(&name);
        let pieces_of_walk = cut_walk(&walk, undecided.get(from).unwrap_or(&none));
        for (i, w) in pieces_of_walk.iter().enumerate() {
            let mut steps = Vec::new();
            for u in w {
                let (unitig, fwd) = parse_oriented(u)?;
                let pieces = unitig_pieces
                    .get(&unitig)
                    .ok_or_else(|| anyhow!("{sample}: unitig {unitig} of {name} not in individual.segments.tsv"))?;
                if fwd {
                    for (p, pf) in pieces {
                        steps.extend(step_of(p, *pf));
                    }
                } else {
                    for (p, pf) in pieces.iter().rev() {
                        steps.extend(step_of(p, !*pf));
                    }
                }
            }
            let label = if pieces_of_walk.len() > 1 { format!("{name}.r{}", i + 1) } else { name.clone() };
            walks.push((label, merge_steps(steps)));
        }
    }
    // junction keys in cohort terms
    let mut junctions = Vec::new();
    for (molecule, left, repeat, right, source, class) in junction_rows {
        if other_configuration.contains(&molecule) {
            continue;
        }
        let mut parts = Vec::new();
        for u in std::iter::once(&left).chain(repeat.iter()).chain(std::iter::once(&right)) {
            let (unitig, fwd) = parse_oriented(u)?;
            let Some(pieces) = unitig_pieces.get(&unitig) else { continue };
            let mut steps = Vec::new();
            if fwd {
                for (p, pf) in pieces {
                    steps.extend(step_of(p, *pf));
                }
            } else {
                for (p, pf) in pieces.iter().rev() {
                    steps.extend(step_of(p, !*pf));
                }
            }
            parts.push(merge_steps(steps));
        }
        if parts.len() < 3 {
            continue;
        }
        let l = parts[0].last().cloned();
        let r = parts[parts.len() - 1].first().cloned();
        let (Some(l), Some(r)) = (l, r) else { continue };
        let rep: Vec<Step> = merge_steps(parts[1..parts.len() - 1].concat());
        junctions.push((junction_key(&l, &rep, &r), source, class));
    }
    // every pairing the graph allowed, in panel terms (one row per key: a bridge and its
    // reverse are one pairing)
    let mut pairings = BTreeMap::new();
    let pairings_file = od.join("pairings.tsv");
    if pairings_file.exists() {
        let pt = Table::read(&pairings_file)?;
        for row in &pt.rows {
            let repeat: Vec<String> = pt.get(row, "repeat")?.split(',').filter(|s| !s.is_empty()).map(str::to_string).collect();
            let named: Vec<String> = std::iter::once(pt.get(row, "left")?.to_string())
                .chain(repeat)
                .chain(std::iter::once(pt.get(row, "right")?.to_string()))
                .collect();
            let mut parts = Vec::new();
            for u in &named {
                let (unitig, fwd) = parse_oriented(u)?;
                let Some(pieces) = unitig_pieces.get(&unitig) else { continue };
                let mut steps = Vec::new();
                if fwd {
                    for (p, pf) in pieces {
                        steps.extend(step_of(p, *pf));
                    }
                } else {
                    for (p, pf) in pieces.iter().rev() {
                        steps.extend(step_of(p, !*pf));
                    }
                }
                parts.push(merge_steps(steps));
            }
            if parts.len() < 3 {
                continue;
            }
            let (Some(l), Some(r)) = (parts[0].last().cloned(), parts[parts.len() - 1].first().cloned()) else { continue };
            let rep: Vec<Step> = merge_steps(parts[1..parts.len() - 1].concat());
            let fields = ["repeat_bp", "reads_crossing", "pairs_crossing", "crossing_fragments", "expected", "p_at_most_observed", "class", "panel_takes"]
                .iter()
                .map(|c| pt.get(row, c).map(str::to_string))
                .collect::<Result<Vec<_>>>()?;
            pairings.entry(junction_key(&l, &rep, &r)).or_insert(fields.join("\t"));
        }
    }
    // presence
    let pt = Table::read(&od.join("presence.tsv"))?;
    let mut presence = FxHashMap::default();
    for row in &pt.rows {
        let seg = pt.get(row, "segment")?.to_string();
        // the run's panel must be this panel: every segment at its length
        let len: usize = pt.get(row, "length")?.parse()?;
        match panel.seqs.get(&seg) {
            Some(q) if q.len() == len => {}
            _ => bail!("{sample}: {o} segment {seg} ({len} bp) is not in the panel graph at that length; the run used another panel"),
        }
        let depth = pt.get(row, "depth")?.parse::<f64>().ok();
        let kmers: usize = pt.get(row, "kmers")?.parse().unwrap_or(0);
        let call = match pt.get(row, "call")? {
            c @ ("present" | "absent") if kmers >= MIN_OWN_KMERS => c,
            _ => "undetermined",
        };
        presence.insert(seg, (depth, call.to_string()));
    }
    // genotypes
    let gt = Table::read(&od.join("genotypes.tsv"))?;
    let mut calls: BTreeMap<u64, (Option<String>, Vec<String>)> = BTreeMap::new();
    for row in &gt.rows {
        let site: u64 = gt.get(row, "site")?.parse()?;
        let allele = gt.get(row, "allele")?.to_string();
        let entry = calls.entry(site).or_insert((None, Vec::new()));
        match gt.get(row, "call")? {
            "major" => entry.0 = Some(allele),
            "minor" => entry.1.push(allele),
            _ => {}
        }
    }
    let missing_pieces = missing.len();
    let relatives = closest_relatives(dir, o)?;
    Ok(Some(SampleOrganelle { walks, presence, new_segments, calls, junctions, pairings, missing_pieces, unmapped_unitigs, has_junctions, relatives }))
}

/// A junction (left anchor step, repeat steps, right anchor step) and its reverse are one.
pub fn junction_key(l: &Step, rep: &[Step], r: &Step) -> String {
    let fwd: Vec<String> = std::iter::once(l.text())
        .chain(std::iter::once("|".into()))
        .chain(rep.iter().map(Step::text))
        .chain(std::iter::once("|".into()))
        .chain(std::iter::once(r.text()))
        .collect();
    let rev: Vec<String> = std::iter::once(r.flip().text())
        .chain(std::iter::once("|".into()))
        .chain(rep.iter().rev().map(|s| s.flip().text()))
        .chain(std::iter::once("|".into()))
        .chain(std::iter::once(l.flip().text()))
        .collect();
    let f = fwd.join(",").replace(",|,", " | ");
    let r = rev.join(",").replace(",|,", " | ");
    if f <= r { f } else { r }
}

/// Site metadata from the first sample's genotypes.tsv rows: anchors, alleles, their panel runs.
#[derive(Default, Clone)]
pub struct Site {
    pub left: String,
    pub right: String,
    pub alleles: Vec<String>,
    pub runs: BTreeMap<String, Vec<String>>,
}

pub fn read_sites(path: &Path) -> Result<BTreeMap<u64, Site>> {
    let gt = Table::read(path)?;
    let mut sites: BTreeMap<u64, Site> = BTreeMap::new();
    for row in &gt.rows {
        let site: u64 = gt.get(row, "site")?.parse()?;
        let s = sites.entry(site).or_default();
        s.left = gt.get(row, "left")?.to_string();
        s.right = gt.get(row, "right")?.to_string();
        let allele = gt.get(row, "allele")?.to_string();
        if !s.alleles.contains(&allele) {
            s.alleles.push(allele.clone());
        }
        let runs: Vec<String> = gt.get(row, "reference_runs")?.split(',').filter(|r| !r.is_empty()).map(str::to_string).collect();
        s.runs.insert(allele, runs);
    }
    Ok(sites)
}

/// Nearest panel runs of one sample by its major alleles: (run, distance, sites compared),
/// nearest first. A run is compared at the sites where some allele lists it.
pub fn nearest_runs(calls: &BTreeMap<u64, (Option<String>, Vec<String>)>, sites: &BTreeMap<u64, Site>) -> Vec<(String, f64, usize)> {
    let mut tally: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    for (site, (major, _)) in calls {
        let (Some(major), Some(s)) = (major, sites.get(site)) else { continue };
        for (allele, runs) in &s.runs {
            for run in runs {
                let t = tally.entry(run.clone()).or_insert((0, 0));
                t.1 += 1;
                if allele == major {
                    t.0 += 1;
                }
            }
        }
    }
    let mut out: Vec<(String, f64, usize)> = tally
        .into_iter()
        .filter(|(_, (_, n))| *n > 0)
        .map(|(r, (agree, n))| (r, 1.0 - agree as f64 / n as f64, n))
        .collect();
    out.sort_by(|a, b| a.1.total_cmp(&b.1).then(b.2.cmp(&a.2)).then(a.0.cmp(&b.0)));
    out
}

/// Distance between two samples' major alleles at the sites both call.
pub fn sample_distance(a: &BTreeMap<u64, (Option<String>, Vec<String>)>, b: &BTreeMap<u64, (Option<String>, Vec<String>)>) -> (f64, usize) {
    let (mut n, mut diff) = (0usize, 0usize);
    for (site, (ma, _)) in a {
        if let (Some(ma), Some((Some(mb), _))) = (ma, b.get(site)) {
            n += 1;
            if ma != mb {
                diff += 1;
            }
        }
    }
    if n == 0 { (f64::NAN, 0) } else { (diff as f64 / n as f64, n) }
}

/// Sequence of an allele walk (e.g. "3+" or "71-,239-") in the panel.
fn walk_seq(panel: &Panel, walk: &str) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    // "-" is the empty allele (a deletion between the anchors)
    for s in walk.split(',').filter(|s| !s.is_empty() && *s != "-") {
        let (name, fwd) = parse_oriented(s)?;
        let seq = panel.seqs.get(&name).ok_or_else(|| anyhow!("allele segment {name} not in the panel"))?;
        if fwd { out.extend_from_slice(seq) } else { out.extend(revcomp(seq)) }
    }
    Ok(out)
}

/// A VCF position for a site: on a panel path holding its left anchor, the 1-based position of
/// the anchor's last base, that base, and the path name.
fn site_position(panel: &Panel, site: &Site) -> Result<Option<(String, usize, u8, Vec<u8>)>> {
    let (left, lf) = parse_oriented(&site.left)?;
    // the first allele's first run decides the contig
    for runs in site.runs.values() {
        for run in runs {
            let Some((pname, steps)) = panel.paths.iter().find(|(n, _)| n.split('#').next() == Some(run.as_str())) else { continue };
            let mut pos = 0usize;
            for (i, s) in steps.iter().enumerate() {
                let len = panel.seqs.get(&s.seg).map_or(0, Vec::len);
                if s.seg == left && s.fwd == lf {
                    let seq = walk_seq(panel, &s.text())?;
                    let anchor = *seq.last().unwrap_or(&b'N');
                    // the run's own allele: the steps between the anchors
                    let (right, rf) = parse_oriented(&site.right)?;
                    let mut own = Vec::new();
                    for t in &steps[i + 1..] {
                        if t.seg == right && t.fwd == rf {
                            break;
                        }
                        own.push(t.text());
                    }
                    let own_seq = walk_seq(panel, &own.join(","))?;
                    return Ok(Some((pname.clone(), pos + len, anchor, own_seq)));
                }
                pos += len;
            }
        }
    }
    Ok(None)
}

pub struct CohortArgs {
    pub samples: Vec<(String, PathBuf)>,
    pub out: PathBuf,
    pub organelles: Vec<String>,
    pub references: FxHashMap<String, PathBuf>,
    /// sample -> (median read length, data_limited); samples absent here count.
    pub read_lengths: FxHashMap<String, (Option<f64>, bool)>,
    /// Samples whose median read length is below this (or that are data_limited) keep their
    /// own columns and paths but do not count in any frequency.
    pub min_read_length: f64,
}

/// Reads the samples TSV: sample<TAB>pan reconstruct output directory.
pub fn read_samples(path: &Path) -> Result<Vec<(String, PathBuf)>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut out = Vec::new();
    for line in text.lines().filter(|l| !l.trim().is_empty() && !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 2 {
            bail!("{}: expected sample<TAB>directory, got {line:?}", path.display());
        }
        if f[0] == "sample" {
            continue;
        }
        out.push((f[0].to_string(), PathBuf::from(f[1])));
    }
    if out.is_empty() {
        bail!("{}: no samples", path.display());
    }
    Ok(out)
}

fn reference_of(dir: &Path, o: &str) -> Result<Option<PathBuf>> {
    let p = dir.join("reference.json");
    if !p.exists() {
        return Ok(None);
    }
    let v: serde_json::Value = serde_json::from_str(&fs::read_to_string(&p)?)?;
    Ok(v[o]["gfa"].as_str().map(PathBuf::from))
}

pub fn run(args: CohortArgs) -> Result<serde_json::Value> {
    let t0 = std::time::Instant::now();
    fs::create_dir_all(&args.out)?;
    let mut report = serde_json::json!({"tool": "ovasm pan cohort", "version": env!("CARGO_PKG_VERSION"),
        "samples": args.samples.iter().map(|(s, d)| serde_json::json!({"sample": s, "dir": d})).collect::<Vec<_>>(),
        "organelles": {}});
    // short-read samples: their own calls stand, but they are left out of every frequency
    let excluded: FxHashSet<String> = args
        .samples
        .iter()
        .filter(|(s, _)| args.read_lengths.get(s).is_some_and(|(len, lim)| *lim || len.is_some_and(|l| l < args.min_read_length)))
        .map(|(s, _)| s.clone())
        .collect();
    report["min_read_length"] = serde_json::json!(args.min_read_length);
    report["excluded_from_frequencies"] = serde_json::json!(args
        .samples
        .iter()
        .filter(|(s, _)| excluded.contains(s))
        .map(|(s, _)| serde_json::json!({"sample": s, "read_length_median": args.read_lengths[s].0, "data_limited": args.read_lengths[s].1}))
        .collect::<Vec<_>>());
    report["samples_without_read_length"] = serde_json::json!(args.samples.iter().filter(|(s, _)| !args.read_lengths.contains_key(s)).count());
    for o in &args.organelles {
        let reference = match args.references.get(o) {
            Some(p) => p.clone(),
            None => {
                let mut refs: Vec<PathBuf> = Vec::new();
                for (s, d) in &args.samples {
                    if let Some(r) = reference_of(d, o)? {
                        if !refs.contains(&r) {
                            refs.push(r);
                        }
                    } else if d.join(o).exists() {
                        bail!("{s}: no {o} graph in reference.json; give --reference {o}=<gfa>");
                    }
                }
                match refs.len() {
                    0 => continue,
                    1 => refs.remove(0),
                    _ => bail!("{o}: the samples were reconstructed on different panel graphs ({}); give --reference {o}=<gfa> only if they are the same graph", refs.iter().map(|p| p.display().to_string()).collect::<Vec<_>>().join(", ")),
                }
            }
        };
        let summary = run_organelle(o, &reference, &args.samples, &args.out.join(o), &excluded, &args.read_lengths)?;
        report["organelles"][o] = summary;
    }
    report["elapsed_seconds"] = serde_json::json!(t0.elapsed().as_secs_f64());
    report["peak_rss_mb"] = serde_json::json!(peak_rss_mb());
    fs::write(args.out.join("cohort.json"), serde_json::to_string_pretty(&report)? + "\n")?;
    Ok(report)
}

fn peak_rss_mb() -> Option<f64> {
    let s = fs::read_to_string("/proc/self/status").ok()?;
    let line = s.lines().find(|l| l.starts_with("VmHWM:"))?;
    let kb: f64 = line.split_whitespace().nth(1)?.parse().ok()?;
    Some(kb / 1024.0)
}

fn run_organelle(
    o: &str,
    reference: &Path,
    samples: &[(String, PathBuf)],
    out: &Path,
    excluded: &FxHashSet<String>,
    read_lengths: &FxHashMap<String, (Option<f64>, bool)>,
) -> Result<serde_json::Value> {
    let counted = |s: &str| !excluded.contains(s);
    fs::create_dir_all(out)?;
    let panel = Panel::read(reference)?;
    let mut new = NewSegments::default();
    let mut runs: Vec<(String, SampleOrganelle)> = Vec::new();
    let mut sites: BTreeMap<u64, Site> = BTreeMap::new();
    let mut skipped = Vec::new();
    for (s, d) in samples {
        match read_sample(s, d, o, &panel, &mut new)? {
            Some(r) => {
                if sites.is_empty() {
                    sites = read_sites(&d.join(o).join("genotypes.tsv"))?;
                }
                runs.push((s.clone(), r));
            }
            None => skipped.push(s.clone()),
        }
    }
    // cohort.gfa
    let mut links: BTreeMap<(Step, Step), (bool, FxHashSet<String>)> = BTreeMap::new();
    for l in &panel.links {
        links.entry(l.clone()).or_insert((true, FxHashSet::default()));
    }
    for (s, r) in &runs {
        for (_, w) in &r.walks {
            for pair in w.windows(2) {
                let key = canonical_link(&pair[0], &pair[1]);
                links.entry(key).or_insert((false, FxHashSet::default())).1.insert(s.clone());
            }
        }
    }
    {
        let mut g = BufWriter::new(fs::File::create(out.join("cohort.gfa"))?);
        writeln!(g, "H\tVN:Z:1.0")?;
        for name in &panel.names {
            let n = runs.iter().filter(|(s, r)| counted(s) && r.presence.get(name).is_some_and(|p| p.1 == "present")).count();
            writeln!(g, "S\t{name}\t{}\tLN:i:{}\tNS:i:{n}\tSO:Z:panel", String::from_utf8_lossy(&panel.seqs[name]), panel.seqs[name].len())?;
        }
        for (id, seq) in &new.seqs {
            let n = new.samples.get(id).map_or(0, |s| s.iter().filter(|x| counted(x)).count());
            writeln!(g, "S\t{id}\t{}\tLN:i:{}\tNS:i:{n}\tSO:Z:new", String::from_utf8_lossy(seq), seq.len())?;
        }
        for ((a, b), (in_panel, who)) in &links {
            let src = if *in_panel { "panel" } else { "reads" };
            writeln!(g, "L\t{}\t{}\t{}\t{}\t0M\tSO:Z:{src}\tRS:i:{}", a.seg, if a.fwd { '+' } else { '-' }, b.seg, if b.fwd { '+' } else { '-' }, who.iter().filter(|x| counted(x)).count())?;
        }
        for (s, r) in &runs {
            for (name, w) in &r.walks {
                if w.is_empty() {
                    continue;
                }
                let steps: Vec<String> = w.iter().map(Step::text).collect();
                let overlaps = vec!["0M"; w.len().saturating_sub(1)].join(",");
                writeln!(g, "P\t{s}#1#{name}\t{}\t{}", steps.join(","), if overlaps.is_empty() { "*".into() } else { overlaps })?;
            }
        }
    }
    // presence
    {
        let mut d = BufWriter::new(fs::File::create(out.join("presence.tsv"))?);
        let mut c = BufWriter::new(fs::File::create(out.join("presence.calls.tsv"))?);
        let head = std::iter::once("segment".to_string()).chain(std::iter::once("length".into())).chain(std::iter::once("source".into())).chain(runs.iter().map(|(s, _)| s.clone())).collect::<Vec<_>>().join("\t");
        writeln!(d, "{head}")?;
        writeln!(c, "{head}")?;
        for name in &panel.names {
            let len = panel.seqs[name].len();
            let (mut dv, mut cv) = (vec![name.clone(), len.to_string(), "panel".into()], vec![name.clone(), len.to_string(), "panel".into()]);
            for (_, r) in &runs {
                match r.presence.get(name) {
                    Some((depth, call)) => {
                        dv.push(match (call.as_str(), depth) {
                            ("absent", _) => "0".into(),
                            ("present", Some(x)) => format!("{x}"),
                            _ => "NA".into(),
                        });
                        cv.push(call.clone());
                    }
                    None => {
                        dv.push("NA".into());
                        cv.push("NA".into());
                    }
                }
            }
            writeln!(d, "{}", dv.join("\t"))?;
            writeln!(c, "{}", cv.join("\t"))?;
        }
        for (id, seq) in &new.seqs {
            let (mut dv, mut cv) = (vec![id.clone(), seq.len().to_string(), "new".into()], vec![id.clone(), seq.len().to_string(), "new".into()]);
            for (_, r) in &runs {
                let has = r.new_segments.contains(id);
                dv.push(if has { "1".into() } else { "0".into() });
                cv.push(if has { "present".into() } else { "absent".into() });
            }
            writeln!(d, "{}", dv.join("\t"))?;
            writeln!(c, "{}", cv.join("\t"))?;
        }
    }
    // genotypes
    let mut vcf_sites = 0usize;
    {
        let mut t = BufWriter::new(fs::File::create(out.join("genotypes.tsv"))?);
        writeln!(t, "site\tleft\tright\talleles\t{}", runs.iter().map(|(s, _)| s.as_str()).collect::<Vec<_>>().join("\t"))?;
        let mut v = BufWriter::new(fs::File::create(out.join("genotypes.vcf"))?);
        writeln!(v, "##fileformat=VCFv4.2")?;
        writeln!(v, "##source=ovasm pan cohort {}", env!("CARGO_PKG_VERSION"))?;
        writeln!(v, "##reference={}", reference.display())?;
        for (pname, steps) in &panel.paths {
            let len: usize = steps.iter().map(|s| panel.seqs.get(&s.seg).map_or(0, Vec::len)).sum();
            writeln!(v, "##contig=<ID={pname},length={len}>")?;
        }
        writeln!(v, "##INFO=<ID=SITE,Number=1,Type=Integer,Description=\"Site of the panel's genotypes.tsv\">")?;
        writeln!(v, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Share of called samples with each ALT allele as major (samples counted in frequencies)\">")?;
        writeln!(v, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Samples counted in frequencies with a major call\">")?;
        writeln!(v, "##INFO=<ID=WALK,Number=.,Type=String,Description=\"Panel segment walk of REF then each ALT\">")?;
        writeln!(v, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Major allele (haploid)\">")?;
        writeln!(v, "##FORMAT=<ID=MN,Number=1,Type=String,Description=\"Minor alleles seen (indices)\">")?;
        writeln!(v, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}", runs.iter().map(|(s, _)| s.as_str()).collect::<Vec<_>>().join("\t"))?;
        for (id, site) in &sites {
            let idx = |a: &str| site.alleles.iter().position(|x| x == a);
            let cols: Vec<String> = runs
                .iter()
                .map(|(_, r)| match r.calls.get(id) {
                    Some((Some(m), minors)) => {
                        let mut s = idx(m).map_or(".".into(), |i| i.to_string());
                        if !minors.is_empty() {
                            s.push('*');
                        }
                        s
                    }
                    _ => ".".into(),
                })
                .collect();
            writeln!(t, "{id}\t{}\t{}\t{}\t{}", site.left, site.right, site.alleles.join(";"), cols.join("\t"))?;
            // VCF: the alleles' sequences after the left anchor's last base
            let Some((contig, pos, anchor, own)) = site_position(&panel, site)? else { continue };
            let seqs: Vec<Vec<u8>> = site.alleles.iter().map(|a| walk_seq(&panel, a)).collect::<Result<_>>()?;
            let ref_i = seqs.iter().position(|s| *s == own).unwrap_or(0);
            let mut order: Vec<usize> = vec![ref_i];
            order.extend((0..seqs.len()).filter(|&i| i != ref_i));
            if order.len() < 2 {
                continue;
            }
            let allele_text = |i: usize| {
                let mut s = vec![anchor];
                s.extend_from_slice(&seqs[i]);
                String::from_utf8_lossy(&s).to_string()
            };
            let mut counts = vec![0usize; order.len()];
            let mut gts = Vec::new();
            for (s, r) in &runs {
                match r.calls.get(id) {
                    Some((Some(m), minors)) => {
                        let i = idx(m).and_then(|i| order.iter().position(|&x| x == i));
                        if let Some(i) = i.filter(|_| counted(s)) {
                            counts[i] += 1;
                        }
                        let mn: Vec<String> = minors.iter().filter_map(|a| idx(a).and_then(|i| order.iter().position(|&x| x == i))).map(|i| i.to_string()).collect();
                        gts.push(format!("{}:{}", i.map_or(".".into(), |i| i.to_string()), if mn.is_empty() { ".".into() } else { mn.join(",") }));
                    }
                    _ => gts.push(".:.".into()),
                }
            }
            let ns: usize = counts.iter().sum();
            let af: Vec<String> = counts[1..].iter().map(|&c| if ns == 0 { "0".into() } else { format!("{:.4}", c as f64 / ns as f64) }).collect();
            let walks: Vec<String> = order.iter().map(|&i| site.alleles[i].replace(',', "/")).collect();
            writeln!(
                v,
                "{contig}\t{pos}\tsite{id}\t{}\t{}\t.\tPASS\tSITE={id};AF={};NS={ns};WALK={}\tGT:MN\t{}",
                allele_text(order[0]),
                order[1..].iter().map(|&i| allele_text(i)).collect::<Vec<_>>().join(","),
                af.join(","),
                walks.join(","),
                gts.join("\t")
            )?;
            vcf_sites += 1;
        }
    }
    // junctions
    let mut jcount: BTreeMap<String, BTreeMap<&str, Vec<String>>> = BTreeMap::new();
    for (s, r) in runs.iter().filter(|(s, _)| counted(s)) {
        for (key, source, class) in &r.junctions {
            let kind = match (source.as_str(), class.as_str()) {
                ("contradicted_by_reads", _) | ("reads", "unsupported") => "reads_unsupported",
                ("reads", _) => "reads_supported",
                ("panel", _) => "panel_inferred",
                _ => "undetermined",
            };
            jcount.entry(key.clone()).or_default().entry(kind).or_default().push(s.clone());
        }
    }
    {
        let n = runs.iter().filter(|(s, _)| counted(s)).count().max(1) as f64;
        let mut j = BufWriter::new(fs::File::create(out.join("junctions.tsv"))?);
        writeln!(j, "left\trepeat\tright\treads_supported\treads_unsupported\tpanel_inferred\tundetermined\treads_frequency\tpanel_frequency\tsamples_reads\tsamples_panel")?;
        for (key, kinds) in &jcount {
            let parts: Vec<&str> = key.split(" | ").collect();
            let c = |k: &str| kinds.get(k).map_or(0, Vec::len);
            let names = |k: &str| kinds.get(k).map_or("-".to_string(), |v| v.join(","));
            writeln!(
                j,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}",
                parts.first().unwrap_or(&""),
                parts.get(1).unwrap_or(&""),
                parts.get(2).unwrap_or(&""),
                c("reads_supported"),
                c("reads_unsupported"),
                c("panel_inferred"),
                c("undetermined"),
                c("reads_supported") as f64 / n,
                c("panel_inferred") as f64 / n,
                names("reads_supported"),
                names("panel_inferred")
            )?;
        }
    }
    // pairings: per sample, every pairing its graph allowed, in panel terms
    {
        let mut f = BufWriter::new(fs::File::create(out.join("pairings.tsv"))?);
        writeln!(f, "left\trepeat\tright\tsample\tin_frequencies\trepeat_bp\treads_crossing\tpairs_crossing\tcrossing_fragments\texpected\tp_at_most_observed\tclass\tpanel_takes")?;
        let mut rows: Vec<(&String, &String, &String)> = runs.iter().flat_map(|(s, r)| r.pairings.iter().map(move |(k, v)| (k, s, v))).collect();
        rows.sort();
        for (key, s, v) in rows {
            let parts: Vec<&str> = key.split(" | ").collect();
            writeln!(f, "{}\t{}\t{}\t{}\t{}\t{}", parts.first().unwrap_or(&""), parts.get(1).unwrap_or(&""), parts.get(2).unwrap_or(&""), s, if counted(s) { "yes" } else { "no" }, v)?;
        }
    }
    // nearest
    {
        let mut f = BufWriter::new(fs::File::create(out.join("nearest.tsv"))?);
        writeln!(f, "sample\tnearest_panel_run\tdistance\tsites\tsecond_panel_run\tsecond_distance\tsecond_sites\tnearest_cohort_sample\tcohort_distance\tcohort_sites\tclosest_relatives\tread_length_median\tin_frequencies")?;
        for (i, (s, r)) in runs.iter().enumerate() {
            let near = nearest_runs(&r.calls, &sites);
            let fmt = |k: usize| near.get(k).map_or(("-".to_string(), "NA".to_string(), "0".to_string()), |(n, d, c)| (n.clone(), format!("{d:.5}"), c.to_string()));
            let (a, ad, an) = fmt(0);
            let (b, bd, bn) = fmt(1);
            let mut best: Option<(String, f64, usize)> = None;
            for (j, (t, q)) in runs.iter().enumerate() {
                if i == j {
                    continue;
                }
                let (d, n) = sample_distance(&r.calls, &q.calls);
                if n > 0 && best.as_ref().map_or(true, |b| d < b.1) {
                    best = Some((t.clone(), d, n));
                }
            }
            let (c, cd, cn) = best.map_or(("-".to_string(), "NA".to_string(), "0".to_string()), |(t, d, n)| (t, format!("{d:.5}"), n.to_string()));
            // the run's own posterior over panel members, by name
            let rel = if r.relatives.is_empty() {
                "-".to_string()
            } else {
                r.relatives.iter().map(|(n, w)| format!("{n}={w:.3e}")).collect::<Vec<_>>().join(",")
            };
            let len = read_lengths.get(s).and_then(|x| x.0).map_or("NA".to_string(), |l| format!("{l}"));
            let inf = if counted(s) { "yes" } else { "no" };
            writeln!(f, "{s}\t{a}\t{ad}\t{an}\t{b}\t{bd}\t{bn}\t{c}\t{cd}\t{cn}\t{rel}\t{len}\t{inf}")?;
        }
    }
    let paths: usize = runs.iter().map(|(_, r)| r.walks.len()).sum();
    Ok(serde_json::json!({
        "reference": reference, "samples": runs.len(), "samples_in_frequencies": runs.iter().filter(|(s, _)| counted(s)).count(),
        "skipped_without_this_organelle": skipped,
        "panel_segments": panel.names.len(), "new_segments": new.seqs.len(),
        "new_segment_bases": new.seqs.iter().map(|(_, s)| s.len()).sum::<usize>(),
        "links": links.len(), "read_links_not_in_panel": links.values().filter(|(p, _)| !p).count(),
        "paths": paths, "sites": sites.len(),
        "samples_without_molecules": runs.iter().filter(|(_, r)| r.walks.is_empty()).map(|(s, _)| s.clone()).collect::<Vec<_>>(),
        "samples_without_junctions_tsv": runs.iter().filter(|(_, r)| !r.has_junctions).map(|(s, _)| s.clone()).collect::<Vec<_>>(),
        "unmapped_unitigs_taken_as_new": runs.iter().filter(|(_, r)| r.unmapped_unitigs > 0).map(|(s, r)| serde_json::json!({"sample": s, "unitigs": r.unmapped_unitigs})).collect::<Vec<_>>(),
        "pieces_without_sequence": runs.iter().filter(|(_, r)| r.missing_pieces > 0).map(|(s, r)| serde_json::json!({"sample": s, "pieces": r.missing_pieces})).collect::<Vec<_>>(), "vcf_sites": vcf_sites, "junctions": jcount.len(),
        "outputs": ["cohort.gfa", "presence.tsv", "presence.calls.tsv", "genotypes.tsv", "genotypes.vcf", "junctions.tsv", "pairings.tsv", "nearest.tsv"],
        "new_segment_merge": "exact sequence in either orientation; near-identical new segments stay apart",
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn st(s: &str) -> Step {
        let (seg, fwd) = parse_oriented(s).unwrap();
        Step { seg, fwd }
    }

    #[test]
    fn walks_are_cut_where_the_panel_decided_a_join_either_way_round() {
        let w: Vec<String> = ["u1+", "u2+", "u3+", "u4+", "u5+"].iter().map(|s| s.to_string()).collect();
        let cut = cut_walk(&w, &[("u2+".into(), vec!["u3+".into()], "u4+".into())]);
        assert_eq!(cut, vec![vec!["u1+", "u2+", "u3+"], vec!["u4+", "u5+"]]);
        // the same join written from the other strand
        let cut = cut_walk(&w, &[("u4-".into(), vec!["u3-".into()], "u2-".into())]);
        assert_eq!(cut, vec![vec!["u1+", "u2+"], vec!["u3+", "u4+", "u5+"]]);
        assert_eq!(cut_walk(&w, &[]), vec![w.clone()]);
    }

    #[test]
    fn new_segments_merge_in_either_orientation() {
        let mut n = NewSegments::default();
        let a = n.step(b"AACG", "s1");
        let b = n.step(b"CGTT", "s2");
        assert_eq!(a.seg, b.seg);
        assert_ne!(a.fwd, b.fwd);
        assert_eq!(n.seqs.len(), 1);
        assert_eq!(n.samples[&a.seg].len(), 2);
    }

    #[test]
    fn junction_keys_and_links_are_strand_independent() {
        let k1 = junction_key(&st("1+"), &[st("2+")], &st("3+"));
        let k2 = junction_key(&st("3-"), &[st("2-")], &st("1-"));
        assert_eq!(k1, k2);
        assert_eq!(canonical_link(&st("1+"), &st("2+")), canonical_link(&st("2-"), &st("1-")));
    }

    #[test]
    fn split_panel_segments_merge_back_into_one_step() {
        let s = merge_steps(vec![st("5+"), st("5+"), st("6-"), st("n1+"), st("n1+")]);
        assert_eq!(s, vec![st("5+"), st("6-"), st("n1+"), st("n1+")]);
    }

    #[test]
    fn nearest_runs_count_agreement_at_the_sites_each_run_has() {
        let mut sites = BTreeMap::new();
        let mut s = Site { left: "1+".into(), right: "4+".into(), alleles: vec!["2+".into(), "3+".into()], runs: BTreeMap::new() };
        s.runs.insert("2+".into(), vec!["A".into(), "B".into()]);
        s.runs.insert("3+".into(), vec!["C".into()]);
        sites.insert(0, s.clone());
        s.runs.insert("2+".into(), vec!["A".into()]);
        s.runs.insert("3+".into(), vec!["B".into(), "C".into()]);
        sites.insert(1, s);
        let mut calls = BTreeMap::new();
        calls.insert(0, (Some("2+".to_string()), vec![]));
        calls.insert(1, (Some("2+".to_string()), vec![]));
        let near = nearest_runs(&calls, &sites);
        assert_eq!(near[0].0, "A");
        assert_eq!(near[0].1, 0.0);
        assert_eq!(near[1].0, "B");
        assert!((near[1].1 - 0.5).abs() < 1e-9);
        assert_eq!(near[2].1, 1.0);
    }

    #[test]
    fn a_toy_cohort_writes_every_table() {
        let dir = std::env::temp_dir().join(format!("ovasm-cohort-test-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        let panel = "H\tVN:Z:1.0\nS\t1\tAAAAC\nS\t2\tCCCCG\nS\t3\tGGGGT\nS\t4\tTTTTA\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nP\tR1#1#mt\t1+,2+,4+\t0M,0M\nP\tR2#1#mt\t1+,3+,4+\t0M,0M\n";
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join("panel.gfa"), panel).unwrap();
        for (s, allele, novel) in [("s1", "2+", "ACGTACGTTT"), ("s2", "3+", "AAACGTACGT")] {
            let d = dir.join(s);
            let o = d.join("mitochondrion");
            fs::create_dir_all(&o).unwrap();
            fs::write(d.join("reference.json"), format!("{{\"mitochondrion\":{{\"gfa\":\"{}\"}}}}", dir.join("panel.gfa").display())).unwrap();
            fs::write(o.join("individual.gfa"), format!("H\tVN:Z:1.0\nS\t1\tAAAAC\nS\t2\t{}\nS\t3\tTTTTA\nS\t4\t{novel}\n", if allele == "2+" { "CCCCG" } else { "GGGGT" })).unwrap();
            let refn = if allele == "2+" { 1 } else { 2 };
            fs::write(o.join("individual.segments.tsv"), format!("sample\tsegment\tlength\tpieces\n{s}_reference_supported\tref.0.0.5\t5\t1+\n{s}_reference_supported\tref.{refn}.0.5\t5\t2+\n{s}_reference_supported\tref.3.0.5\t5\t3+\n{s}_denovo\t{s}_denovo.u1\t15\t1+,2+,3+\n{s}_denovo\t{s}_denovo.u2\t10\t4+\n")).unwrap();
            fs::write(o.join("junctions.tsv"), "molecule\tleft\trepeat\tright\tclass\tsource\nm1\tu1+\tu2+\tu1+\tundetermined_long_repeat\tpanel\n").unwrap();
            // a bridge and its reverse: one pairing
            fs::write(o.join("pairings.tsv"), "left\trepeat\tright\trepeat_bp\treads_crossing\tpairs_crossing\tcrossing_fragments\texpected\tp_at_most_observed\tclass\tpanel_probability\truns_holding\truns_informative\tpanel_takes\nu1+\tu2+\tu1+\t10\t1\t2\t3\t4.00\t4.3e-1\tsupported\t0.5\t1\t2\tyes\nu1-\tu2-\tu1-\t10\t1\t2\t3\t4.00\t4.3e-1\tsupported\t0.5\t1\t2\tyes\n").unwrap();
            fs::write(o.join("presence.tsv"), format!("segment\tlength\tkmers\tdepth\trelative\tcall\n1\t5\t30\t30\t1\tpresent\n2\t5\t30\t{}\t1\t{}\n3\t5\t30\t{}\t1\t{}\n4\t5\t5\t30\t1\tpresent\n", if refn == 1 { "30" } else { "0" }, if refn == 1 { "present" } else { "absent" }, if refn == 2 { "30" } else { "0" }, if refn == 2 { "present" } else { "absent" })).unwrap();
            fs::write(o.join("genotypes.tsv"), format!("site\tleft\tright\tallele\tallele_bp\treference_runs\tkmers\tdepth\tshare\tcall\n0\t1+\t4+\t2+\t5\tR1\t1\t30\t1\t{}\n0\t1+\t4+\t3+\t5\tR2\t1\t30\t1\t{}\n", if refn == 1 { "major" } else { "absent" }, if refn == 2 { "major" } else { "absent" })).unwrap();
            fs::write(d.join("mitochondrion.fasta"), ">m1 circular=false path=u1+,u2+,u1+ configuration=1 of=2\nN\n>m2 circular=false path=u1+ configuration=2 of=2\nN\n").unwrap();
        }
        fs::write(dir.join("samples.tsv"), format!("s1\t{}\ns2\t{}\n", dir.join("s1").display(), dir.join("s2").display())).unwrap();
        let samples = read_samples(&dir.join("samples.tsv")).unwrap();
        fs::write(dir.join("s1/reconstruction.json"), r#"{"organelles":{"mitochondrion":{"final":{"pan_posterior":{"closest_relatives":[{"run":"R1","weight":0.9},{"run":"R2","weight":0.1}]}}}}}"#).unwrap();
        let args = |out: &str, read_lengths: FxHashMap<String, (Option<f64>, bool)>| CohortArgs { samples: samples.clone(), out: dir.join(out), organelles: vec!["mitochondrion".into()], references: FxHashMap::default(), read_lengths, min_read_length: MIN_READ_LENGTH };
        let rep = run(args("out", FxHashMap::default())).unwrap();
        let o = &rep["organelles"]["mitochondrion"];
        assert_eq!(o["samples"], 2);
        // the two samples' new segment is one sequence read from either strand
        assert_eq!(o["new_segments"], 1);
        let gfa = fs::read_to_string(dir.join("out/mitochondrion/cohort.gfa")).unwrap();
        // each sample's molecule is cut where the panel decided its join: two path pieces
        assert_eq!(gfa.lines().filter(|l| l.starts_with("P\ts1#")).count(), 2);
        assert!(gfa.contains("P\ts1#1#m1.r1\t1+,2+,4+,n"));
        let nearest = fs::read_to_string(dir.join("out/mitochondrion/nearest.tsv")).unwrap();
        assert!(nearest.contains("s1\tR1\t0.00000"));
        assert!(nearest.contains("s2\tR2\t0.00000"));
        let j = fs::read_to_string(dir.join("out/mitochondrion/junctions.tsv")).unwrap();
        // one panel-decided junction a sample, through the new segment read from either strand
        let rows: Vec<&str> = j.lines().skip(1).collect();
        assert_eq!(rows.len(), 2);
        assert!(rows.iter().all(|r| r.contains("\t0\t0\t1\t0\t")), "{j}");
        let pr = fs::read_to_string(dir.join("out/mitochondrion/pairings.tsv")).unwrap();
        // one row per sample, in panel segments, with the run's counts
        let rows: Vec<&str> = pr.lines().skip(1).collect();
        assert_eq!(rows.len(), 2, "{pr}");
        assert!(rows.iter().any(|r| r.contains("\ts1\tyes\t10\t1\t2\t3\t4.00\t4.3e-1\tsupported\tyes")), "{pr}");
        let vcf = fs::read_to_string(dir.join("out/mitochondrion/genotypes.vcf")).unwrap();
        let rec = vcf.lines().find(|l| !l.starts_with('#')).unwrap();
        assert!(rec.starts_with("R1#1#mt\t5\tsite0\tCCCCCG\tCGGGGT"), "{rec}");
        assert!(rec.ends_with("\t0:.\t1:."), "{rec}");
        assert!(rec.contains(";AF=0.5000;NS=2;"), "{rec}");
        // member names with their weights; segment 4 has too few own k-mers to call
        assert!(nearest.contains("\tR1=9.000e-1,R2=1.000e-1\tNA\tyes"), "{nearest}");
        let calls = fs::read_to_string(dir.join("out/mitochondrion/presence.calls.tsv")).unwrap();
        assert!(calls.lines().any(|l| l.starts_with("4\t5\tpanel\tundetermined\tundetermined")), "{calls}");
        assert!(calls.lines().any(|l| l.starts_with("2\t5\tpanel\tpresent\tabsent")), "{calls}");
        // a short-read sample keeps its own columns but leaves every frequency
        let mut lens = FxHashMap::default();
        lens.insert("s2".to_string(), (Some(50.0), true));
        let rep = run(args("short", lens)).unwrap();
        assert_eq!(rep["organelles"]["mitochondrion"]["samples_in_frequencies"], 1);
        assert_eq!(rep["excluded_from_frequencies"][0]["sample"], "s2");
        let vcf = fs::read_to_string(dir.join("short/mitochondrion/genotypes.vcf")).unwrap();
        let rec = vcf.lines().find(|l| !l.starts_with('#')).unwrap();
        assert!(rec.contains(";AF=0.0000;NS=1;") && rec.ends_with("\t0:.\t1:."), "{rec}");
        let j = fs::read_to_string(dir.join("short/mitochondrion/junctions.tsv")).unwrap();
        assert_eq!(j.lines().skip(1).count(), 1, "{j}");
        assert!(fs::read_to_string(dir.join("short/mitochondrion/nearest.tsv")).unwrap().contains("\t50\tno"));
        // a run without molecules (no mitochondrion.fasta, no junctions.tsv) still counts by its
        // presence and genotypes
        fs::remove_file(dir.join("s2/mitochondrion.fasta")).unwrap();
        fs::remove_file(dir.join("s2/mitochondrion/junctions.tsv")).unwrap();
        let rep = run(args("bare", FxHashMap::default())).unwrap();
        let o = &rep["organelles"]["mitochondrion"];
        assert_eq!(o["samples"], 2);
        assert_eq!(o["samples_without_molecules"][0], "s2");
        assert_eq!(o["samples_without_junctions_tsv"][0], "s2");
        let gfa = fs::read_to_string(dir.join("bare/mitochondrion/cohort.gfa")).unwrap();
        assert!(!gfa.lines().any(|l| l.starts_with("P\ts2#")));
        let _ = fs::remove_dir_all(&dir);
    }
}
