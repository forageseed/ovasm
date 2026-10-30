//! Reference-graph assisted short-read reconstruction. Reference bases and edges enter
//! the target graph only where complete target reads match them exactly. This is local
//! sequence/adjacency evidence, not long-range phasing or a certificate of origin.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as FmtWrite;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde_json::{json, Value};

use crate::bait::BaitSet;
use crate::evidence::{flip, parse_gfa, revcomp, Graph, Oriented};
use crate::kmer::CanonicalKmers;
use crate::{close, frontier, pan, pipeline, recruit};

const K: usize = 21; // The existing Illumina recruitment preset.
type Edge = (Oriented, Oriented);

fn edge(a: Oriented, b: Oriented) -> Edge {
    (a, b).min((flip(b), flip(a)))
}

struct Walk {
    name: String,
    steps: Vec<Oriented>,
    starts: Vec<usize>,
    sequence: Vec<u8>,
    circular: bool,
}

struct Reference {
    graph: Graph,
    walks: Vec<Walk>,
    sources: BTreeSet<String>,
    input: PathBuf,
}

impl Reference {
    fn read(input: &Path) -> Result<Self> {
        let graph = parse_gfa(input)?;
        if graph.overlaps.iter().any(|&n| n != 0) {
            bail!("reference must be a blunt OV-GFA/pangenome graph (0M links)");
        }
        let names: BTreeMap<_, _> = graph
            .names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.as_str(), i))
            .collect();
        if names.len() != graph.names.len() {
            bail!("duplicate reference segment names");
        }
        let edges: BTreeSet<_> = graph.links.iter().map(|&(a, b)| edge(a, b)).collect();
        let mut walks = Vec::new();
        let mut sources = BTreeSet::new();
        let mut path_names = BTreeSet::new();
        for line in fs::read_to_string(input)?.lines() {
            let f: Vec<_> = line.split('\t').collect();
            if f[0] != "P" {
                continue;
            }
            if f.len() < 4 || !path_names.insert(f[1].to_string()) {
                bail!("invalid or duplicate reference path");
            }
            let mut steps = Vec::new();
            let mut starts = Vec::new();
            let mut sequence = Vec::new();
            for token in f[2].split(',') {
                let (name, sign) =
                    token.split_at(token.len().checked_sub(1).context("empty path step")?);
                if sign != "+" && sign != "-" {
                    bail!("invalid path orientation");
                }
                let id = *names
                    .get(name)
                    .with_context(|| format!("path references unknown node {name}"))?;
                let step = (id, sign == "+");
                if let Some(&previous) = steps.last() {
                    if !edges.contains(&edge(previous, step)) {
                        bail!("reference path {} crosses an absent edge", f[1]);
                    }
                }
                starts.push(sequence.len());
                sequence.extend(graph.oriented_seq(step));
                steps.push(step);
            }
            let circular = f[4..].contains(&"TP:Z:circular");
            if sequence.is_empty() {
                bail!("empty reference path");
            }
            if circular && !edges.contains(&edge(*steps.last().unwrap(), steps[0])) {
                bail!("circular reference path {} lacks its closing edge", f[1]);
            }
            sources.insert(f[1].split('#').next().unwrap().to_string());
            walks.push(Walk {
                name: f[1].into(),
                steps,
                starts,
                sequence,
                circular,
            });
        }
        if walks.is_empty() {
            bail!("reference graph needs P walks, including all source configurations");
        }
        Ok(Self {
            graph,
            walks,
            sources,
            input: input.into(),
        })
    }

    /// The baited sequences: every P walk (a circular one with its boundary words of length
    /// `k`) and every node, with their FASTA names.
    fn bait_sequences(&self, k: usize) -> Vec<(String, Vec<u8>)> {
        let mut out = Vec::new();
        for (i, w) in self.walks.iter().enumerate() {
            let mut seq = w.sequence.clone();
            if w.circular {
                seq.extend(w.sequence.iter().cycle().take(k - 1));
            }
            out.push((format!("{i} source_path={}", w.name), seq));
        }
        // Also retain graph sequence not represented by the finite source paths.
        for (i, seq) in self.graph.seqs.iter().enumerate() {
            out.push((format!("node_{}", self.graph.names[i]), seq.clone()));
        }
        out
    }

    fn stored_bases(&self) -> u64 {
        self.graph.seqs.iter().map(|s| s.len() as u64).sum()
    }

    fn add_bait(&self, bait: &mut BaitSet, name: &str, out: &Path) -> Result<Value> {
        let t = bait.target(name)?;
        let mut fasta = String::new();
        for (name, seq) in self.bait_sequences(K) {
            bait.add_sequence(t, &seq);
            writeln!(fasta, ">{name}\n{}", std::str::from_utf8(&seq)?)?;
        }
        fs::write(out, fasta)?;
        let basis = self.stored_bases();
        if basis == 0 {
            bail!("empty reference graph");
        }
        // Resource normalization uses each stored graph base once, not the inflated
        // sum of eight genomes / alternative configurations. It is not genome size.
        bait.targets[t].seed_length = basis;
        bait.targets[t]
            .seed_files
            .push(self.input.display().to_string());
        Ok(
            json!({"gfa": self.input, "source_paths": self.walks.len(), "source_run_names": self.sources,
            "bait_fasta": out, "sampling_basis": "stored_graph_bases_once", "sampling_basis_bases": basis,
            "sampling_basis_is_physical_genome_size": false}),
        )
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
struct Footprint {
    intervals: Vec<(usize, usize, usize)>, // node and half-open forward coordinates
    edges: BTreeSet<Edge>,
}

fn union(mut ranges: Vec<(usize, usize, usize)>) -> Vec<(usize, usize, usize)> {
    ranges.sort_unstable();
    let mut merged: Vec<(usize, usize, usize)> = Vec::new();
    for (s, a, b) in ranges {
        if let Some(last) = merged.last_mut() {
            if last.0 == s && a <= last.2 {
                last.2 = last.2.max(b);
                continue;
            }
        }
        merged.push((s, a, b));
    }
    merged
}

fn intersect(a: &Footprint, b: &Footprint) -> Footprint {
    let mut ranges = Vec::new();
    let (mut i, mut j) = (0, 0);
    while i < a.intervals.len() && j < b.intervals.len() {
        let x = a.intervals[i];
        let y = b.intervals[j];
        if x.0 == y.0 && x.1.max(y.1) < x.2.min(y.2) {
            ranges.push((x.0, x.1.max(y.1), x.2.min(y.2)));
        }
        if (x.0, x.2) <= (y.0, y.2) {
            i += 1;
        } else {
            j += 1;
        }
    }
    Footprint {
        intervals: union(ranges),
        edges: a.edges.intersection(&b.edges).copied().collect(),
    }
}

impl Walk {
    fn matches(&self, start: usize, seq: &[u8]) -> bool {
        if !self.circular && start + seq.len() > self.sequence.len() {
            return false;
        }
        seq.iter()
            .enumerate()
            .all(|(i, &b)| self.sequence[(start + i) % self.sequence.len()] == b)
    }

    fn footprint(&self, mut start: usize, mut length: usize, graph: &Graph) -> Footprint {
        let mut fp = Footprint::default();
        while length > 0 {
            let i = self.starts.partition_point(|&p| p <= start) - 1;
            let step = self.steps[i];
            let offset = start - self.starts[i];
            let n = length.min(graph.seqs[step.0].len() - offset);
            let (a, b) = if step.1 {
                (offset, offset + n)
            } else {
                let end = graph.seqs[step.0].len() - offset;
                (end - n, end)
            };
            fp.intervals.push((step.0, a, b));
            length -= n;
            start += n;
            if length > 0 {
                let next = (i + 1) % self.steps.len();
                fp.edges.insert(edge(step, self.steps[next]));
                start %= self.sequence.len();
            }
        }
        fp.intervals = union(fp.intervals);
        fp
    }
}

struct Mapper<'a> {
    reference: &'a Reference,
    index: FxHashMap<u64, Vec<(usize, usize)>>,
}

impl<'a> Mapper<'a> {
    fn new(reference: &'a Reference) -> Self {
        let mut index: FxHashMap<u64, Vec<(usize, usize)>> = FxHashMap::default();
        for (i, w) in reference.walks.iter().enumerate() {
            let mut seq = w.sequence.clone();
            if w.circular {
                seq.extend(w.sequence.iter().cycle().take(K - 1));
            }
            for (p, key) in CanonicalKmers::new(&seq, K) {
                index.entry(key).or_default().push((i, p));
            }
        }
        Self { reference, index }
    }

    fn map(&self, seq: &[u8]) -> BTreeSet<Footprint> {
        let seq = seq.to_ascii_uppercase();
        let mut rarest = None;
        for (p, key) in CanonicalKmers::new(&seq, K) {
            let Some(hits) = self.index.get(&key) else {
                return BTreeSet::new();
            };
            if rarest
                .as_ref()
                .map_or(true, |(_, prev): &(usize, &Vec<(usize, usize)>)| {
                    hits.len() < prev.len()
                })
            {
                rarest = Some((p, hits));
            }
        }
        let Some((p, hits)) = rarest else {
            return BTreeSet::new();
        };
        let rc = revcomp(&seq);
        let mut seen = BTreeSet::new();
        let mut result = BTreeSet::new();
        for &(i, hit) in hits {
            let w = &self.reference.walks[i];
            for (reverse, (q, s)) in [(p, seq.as_slice()), (seq.len() - K - p, rc.as_slice())]
                .into_iter()
                .enumerate()
            {
                let start = if hit >= q {
                    hit - q
                } else if w.circular {
                    (hit + w.sequence.len() - q % w.sequence.len()) % w.sequence.len()
                } else {
                    continue;
                };
                if seen.insert((i, start, reverse)) && w.matches(start, s) {
                    result.insert(w.footprint(start, seq.len(), &self.reference.graph));
                }
            }
        }
        result
    }
}

fn count_ranges(coverage: &mut [Vec<u32>], fp: &Footprint) {
    for &(s, a, b) in &fp.intervals {
        for c in &mut coverage[s][a..b] {
            *c = c.saturating_add(1);
        }
    }
}

fn covered_ranges(coverage: &[u32]) -> Vec<(usize, usize)> {
    let mut result = Vec::new();
    let mut start = None;
    for (i, &c) in coverage.iter().chain(std::iter::once(&0)).enumerate() {
        if c > 0 && start.is_none() {
            start = Some(i);
        }
        if c == 0 {
            if let Some(a) = start.take() {
                result.push((a, i));
            }
        }
    }
    result
}

/// Each read record contributes once to a graph coordinate, even when dozens of
/// source configuration paths spell the same placement. Ambiguous alternatives
/// contribute to possible coverage; only their common intersection is emitted.
fn project(reference: &Reference, reads: &[PathBuf], sample: &str, out: &Path) -> Result<Value> {
    fs::create_dir_all(out)?;
    let mapper = Mapper::new(reference);
    let empty = || {
        reference
            .graph
            .seqs
            .iter()
            .map(|s| vec![0u32; s.len()])
            .collect::<Vec<_>>()
    };
    let mut certain = empty();
    let mut possible = empty();
    let mut edges: BTreeMap<Edge, u64> = BTreeMap::new();
    let (mut records, mut bases, mut matched, mut unique, mut ambiguous) =
        (0u64, 0u64, 0u64, 0u64, 0u64);
    let (batches, reader) = recruit::read_batches(reads);
    for batch in batches {
        let results: Vec<_> = batch
            .par_iter()
            .map(|r| (r.seq.len(), mapper.map(&r.seq)))
            .collect();
        for (len, placements) in results {
            records += 1;
            bases += len as u64;
            if placements.is_empty() {
                continue;
            }
            matched += 1;
            if placements.len() == 1 {
                unique += 1;
            } else {
                ambiguous += 1;
            }
            let mut common = placements.first().unwrap().clone();
            let mut all = Vec::new();
            for p in &placements {
                common = intersect(&common, p);
                all.extend_from_slice(&p.intervals);
            }
            count_ranges(
                &mut possible,
                &Footprint {
                    intervals: union(all),
                    edges: BTreeSet::new(),
                },
            );
            count_ranges(&mut certain, &common);
            for e in common.edges {
                *edges.entry(e).or_default() += 1;
            }
        }
    }
    reader
        .join()
        .map_err(|_| anyhow::anyhow!("reference projection reader panicked"))??;
    let mut gfa =
        format!("H\tVN:Z:1.0\tov:Z:organelleverse.graph.v1\tsn:Z:{sample}_reference_supported\n");
    let mut report =
        String::from("reference_node\tlength\tcertain_covered_bases\tpossible_covered_bases\n");
    let mut intervals = String::from("reference_node\tstart\tend\tstatus\n");
    let mut pieces = vec![Vec::new(); reference.graph.names.len()];
    let (mut certain_bp, mut possible_bp) = (0u64, 0u64);
    for (s, seq) in reference.graph.seqs.iter().enumerate() {
        let c = certain[s].iter().filter(|&&n| n > 0).count();
        let p = possible[s].iter().filter(|&&n| n > 0).count();
        certain_bp += c as u64;
        possible_bp += p as u64;
        writeln!(
            report,
            "{}\t{}\t{c}\t{p}",
            reference.graph.names[s],
            seq.len()
        )?;
        for (a, b) in covered_ranges(&certain[s]) {
            let name = format!("ref.{s}.{a}.{b}");
            writeln!(
                gfa,
                "S\t{name}\t{}\tLN:i:{}\tRC:i:{}",
                std::str::from_utf8(&seq[a..b])?,
                b - a,
                certain[s][a..b].iter().min().unwrap()
            )?;
            pieces[s].push((a, b, name));
        }
        // A run-length report includes both zero coverage and placements that cannot
        // distinguish alternative reference coordinates. Neither enters the target.
        let status = |i| {
            if certain[s][i] > 0 {
                "certain_exact"
            } else if possible[s][i] > 0 {
                "ambiguous_only"
            } else {
                "unsupported"
            }
        };
        let mut a = 0;
        while a < seq.len() {
            let state = status(a);
            let mut b = a + 1;
            while b < seq.len() && status(b) == state {
                b += 1;
            }
            writeln!(intervals, "{}\t{a}\t{b}\t{state}", reference.graph.names[s])?;
            a = b;
        }
    }
    let boundary = |o: Oriented, exit: bool| -> Option<&str> {
        let terminal_end = o.1 == exit;
        pieces[o.0]
            .iter()
            .find(|(a, b, _)| {
                if terminal_end {
                    *b == reference.graph.seqs[o.0].len()
                } else {
                    *a == 0
                }
            })
            .map(|(_, _, n)| n.as_str())
    };
    let mut kept_edges = 0;
    let mut edge_report = String::from("from\tto\tcertain_exact_read_records\temitted\n");
    for &(a, b) in &reference.graph.links {
        let n = edges.get(&edge(a, b)).copied().unwrap_or(0);
        let ends = boundary(a, true).zip(boundary(b, false));
        let kept = n > 0 && ends.is_some();
        if kept {
            let (x, y) = ends.unwrap();
            writeln!(
                gfa,
                "L\t{x}\t{}\t{y}\t{}\t0M\tRC:i:{n}",
                if a.1 { '+' } else { '-' },
                if b.1 { '+' } else { '-' }
            )?;
            kept_edges += 1;
        }
        writeln!(
            edge_report,
            "{}\t{}\t{n}\t{kept}",
            reference.graph.label(a),
            reference.graph.label(b)
        )?;
    }
    fs::write(out.join("observed-reference.gfa"), gfa)?;
    fs::write(out.join("node-coverage.tsv"), report)?;
    fs::write(out.join("reference-intervals.tsv"), intervals)?;
    fs::write(out.join("edge-support.tsv"), edge_report)?;
    let result = json!({"method":"complete_read_exact_match_to_all_reference_paths",
        "reference_gfa":reference.input, "reference_source_paths":reference.walks.len(), "reference_source_run_names":reference.sources,
        "input_read_files":reads, "records_scanned":records,"bases_scanned":bases,"exact_matching_records":matched,
        "unique_graph_placement_records":unique,"ambiguous_graph_placement_records":ambiguous,"unmatched_records":records-matched,
        "reference_stored_bases":reference.graph.seqs.iter().map(|s|s.len() as u64).sum::<u64>(),
        "certain_exact_reference_bases":certain_bp,"possible_exact_reference_bases":possible_bp,
        "reference_supported_links":kept_edges,"long_range_phase_identified":false,
        "nuclear_origin_certified":false,"physical_origin_certified":false,
        "limitations":["Exact complete-read matches exclude sequencing errors and new alleles; the separate de novo graph retains novel target sequence.",
            "Source P walks limit the searched configurations; a new target recombination can be absent from every source path.",
            "Local exact evidence does not establish origin, copy number, whole chromosome phase or physical topology."]});
    fs::write(
        out.join("projection.json"),
        serde_json::to_string_pretty(&result)?,
    )?;
    Ok(result)
}

/// Export maximal nonbranching oriented walks. Branches remain in the GFA;
/// a local cycle is labelled a graph cycle, not a physical circular chromosome.
fn unitigs(g: &Graph) -> Vec<(Vec<Oriented>, bool)> {
    let mut successors: BTreeMap<Oriented, BTreeSet<Oriented>> = BTreeMap::new();
    let mut predecessors: BTreeMap<Oriented, BTreeSet<Oriented>> = BTreeMap::new();
    for &(a, b) in &g.links {
        for (x, y) in [(a, b), (flip(b), flip(a))] {
            successors.entry(x).or_default().insert(y);
            predecessors.entry(y).or_default().insert(x);
        }
    }
    let degree =
        |map: &BTreeMap<Oriented, BTreeSet<Oriented>>, o| map.get(&o).map_or(0, BTreeSet::len);
    let mut used = BTreeSet::new();
    let mut output = Vec::new();
    let mut starts = Vec::new();
    for i in 0..g.names.len() {
        for o in [(i, true), (i, false)] {
            let extend_back = predecessors
                .get(&o)
                .filter(|p| p.len() == 1)
                .and_then(|p| p.first())
                .is_some_and(|&p| degree(&successors, p) == 1);
            if !extend_back {
                starts.push(o);
            }
        }
    }
    starts.extend((0..g.names.len()).map(|i| (i, true)));
    for first in starts {
        if used.contains(&first.0) {
            continue;
        }
        let mut walk = vec![first];
        used.insert(first.0);
        let mut circular = false;
        loop {
            let last = *walk.last().unwrap();
            let Some(next) = successors
                .get(&last)
                .filter(|p| p.len() == 1)
                .and_then(|p| p.first())
                .copied()
            else {
                break;
            };
            if degree(&predecessors, next) != 1 {
                break;
            }
            if next == first {
                circular = true;
                break;
            }
            if !used.insert(next.0) {
                break;
            }
            walk.push(next);
        }
        output.push((walk, circular));
    }
    output
}

/// Long node-internal k-mers cannot join a one-base allele between two already
/// shared flanks. Quotient only exact parallel nodes with identical sequence
/// and identical incoming/outgoing oriented neighbours. Distinct alleles,
/// isolated components and tandem self-loops remain distinct.
fn coalesce_parallel(file: &Path) -> Result<usize> {
    let old = parse_gfa(file)?;
    let text = fs::read_to_string(file)?;
    let reversed: Vec<_> = old.seqs.iter().map(|s| revcomp(s) < *s).collect();
    let seqs: Vec<_> = old
        .seqs
        .iter()
        .enumerate()
        .map(|(i, s)| if reversed[i] { revcomp(s) } else { s.clone() })
        .collect();
    let mut owners: Vec<_> = (0..seqs.len()).collect();
    let map = |o: Oriented, owners: &[usize]| (owners[o.0], o.1 != reversed[o.0]);
    loop {
        let mut incoming: BTreeMap<Oriented, BTreeSet<Oriented>> = BTreeMap::new();
        let mut outgoing: BTreeMap<Oriented, BTreeSet<Oriented>> = BTreeMap::new();
        for &(a, b) in &old.links {
            let a = map(a, &owners);
            let b = map(b, &owners);
            for (x, y) in [(a, b), (flip(b), flip(a))] {
                incoming.entry(y).or_default().insert(x);
                outgoing.entry(x).or_default().insert(y);
            }
        }
        let mut groups: BTreeMap<(Vec<u8>, Vec<Oriented>, Vec<Oriented>), usize> = BTreeMap::new();
        let mut replacement = owners.clone();
        let mut changed = false;
        for i in 0..seqs.len() {
            if owners[i] != i {
                continue;
            }
            let (Some(left), Some(right)) = (incoming.get(&(i, true)), outgoing.get(&(i, true)))
            else {
                continue;
            };
            if left.iter().chain(right).any(|o| o.0 == i) {
                continue;
            }
            let signature = (
                seqs[i].clone(),
                left.iter().copied().collect(),
                right.iter().copied().collect(),
            );
            let representative = *groups.entry(signature).or_insert(i);
            if representative != i {
                replacement[i] = representative;
                changed = true;
            }
        }
        if !changed {
            break;
        }
        for owner in &mut owners {
            *owner = replacement[*owner];
        }
    }
    let names: BTreeMap<_, _> = old
        .names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    let label = |o: Oriented| format!("{}{}", old.names[o.0], if o.1 { '+' } else { '-' });
    let mut output = String::new();
    let mut mapping = String::from("pan_node\tindividual_node\treversed\n");
    for (i, n) in old.names.iter().enumerate() {
        let o = map((i, true), &owners);
        if if o.1 {
            seqs[o.0].clone()
        } else {
            revcomp(&seqs[o.0])
        } != old.seqs[i]
        {
            bail!("parallel-node integration changed an input segment sequence");
        }
        writeln!(mapping, "{n}\t{}\t{}", old.names[o.0], !o.1)?;
        if owners[i] == i {
            writeln!(
                output,
                "S\t{n}\t{}\tLN:i:{}",
                std::str::from_utf8(&seqs[i])?,
                seqs[i].len()
            )?;
        }
    }
    let links: BTreeSet<_> = old
        .links
        .iter()
        .map(|&(a, b)| edge(map(a, &owners), map(b, &owners)))
        .collect();
    for &(a, b) in &links {
        writeln!(
            output,
            "L\t{}\t{}\t{}\t{}\t0M",
            old.names[a.0],
            if a.1 { '+' } else { '-' },
            old.names[b.0],
            if b.1 { '+' } else { '-' }
        )?;
    }
    for line in text.lines() {
        let mut f: Vec<_> = line.split('\t').map(str::to_string).collect();
        if f[0] == "H" {
            output = format!("{line}\n{output}");
        }
        if f[0] != "P" {
            continue;
        }
        let mut walk = Vec::new();
        let mut before = Vec::new();
        let mut after = Vec::new();
        for s in f[2].split(',') {
            let i = *names
                .get(&s[..s.len() - 1])
                .context("unknown integrated path node")?;
            let o = (i, s.ends_with('+'));
            let mapped = map(o, &owners);
            before.extend(old.oriented_seq(o));
            after.extend(if mapped.1 {
                seqs[mapped.0].clone()
            } else {
                revcomp(&seqs[mapped.0])
            });
            walk.push(mapped);
        }
        if before != after || walk.windows(2).any(|p| !links.contains(&edge(p[0], p[1]))) {
            bail!("parallel-node integration changed an input path or adjacency");
        }
        f[2] = walk.into_iter().map(label).collect::<Vec<_>>().join(",");
        writeln!(output, "{}", f.join("\t"))?;
    }
    fs::write(file, output)?;
    fs::write(file.with_extension("parallel-nodes.tsv"), mapping)?;
    Ok(owners
        .iter()
        .enumerate()
        .filter(|&(i, &owner)| i != owner)
        .count())
}

fn export_individual(
    observed: &Path,
    denovo: Option<&Path>,
    sample: &str,
    out: &Path,
) -> Result<Value> {
    let mut inputs = Vec::new();
    if !parse_gfa(observed)?.names.is_empty() {
        inputs.push(observed.into());
    }
    if let Some(p) = denovo {
        if !parse_gfa(p)?.names.is_empty() {
            inputs.push(p.into());
        }
    }
    if inputs.is_empty() {
        bail!("neither target reads nor de novo assembly support any sequence");
    }
    let file = out.join("individual.gfa");
    let merged = pan::run(
        &inputs,
        &pan::PanParams { k: 63 },
        &file,
        &out.join("integration.json"),
    )?;
    let text = fs::read_to_string(&file)?;
    let mut gfa = String::new();
    for line in text.lines() {
        if line.starts_with("H\t") {
            writeln!(
                gfa,
                "H\tVN:Z:1.0\tov:Z:organelleverse.individual.v1\tsn:Z:{sample}\tphase:Z:unresolved"
            )?;
        } else {
            // pan's ns counts input graph provenance, not biological sample counts.
            let fields: Vec<_> = line
                .split('\t')
                .filter(|s| !s.starts_with("ns:i:"))
                .collect();
            writeln!(gfa, "{}", fields.join("\t"))?;
        }
    }
    fs::write(&file, &gfa)?;
    let coalesced_parallel_nodes = coalesce_parallel(&file)?;
    gfa = fs::read_to_string(&file)?;
    let graph = parse_gfa(&file)?;
    let mut fasta = String::new();
    let mut lengths = Vec::new();
    for (i, (walk, circular)) in unitigs(&graph).iter().enumerate() {
        let seq: Vec<_> = walk.iter().flat_map(|&o| graph.oriented_seq(o)).collect();
        let name = format!("{sample}#1#contig.{}", i + 1);
        writeln!(
            fasta,
            ">{name} length={} graph_cycle={circular} physical_topology_identified=false\n{}",
            seq.len(),
            std::str::from_utf8(&seq)?
        )?;
        writeln!(
            gfa,
            "P\t{name}\t{}\t*\trole:Z:nonbranching_contig",
            walk.iter()
                .map(|&o| graph.label(o))
                .collect::<Vec<_>>()
                .join(",")
        )?;
        lengths.push(seq.len());
    }
    fs::write(&file, gfa)?;
    fs::write(out.join("individual.fasta"), fasta)?;
    Ok(
        json!({"individual_gfa":file,"individual_fasta":out.join("individual.fasta"),"contig_lengths":lengths,
        "segments":graph.names.len(),"links":graph.links.len(),"stored_bases":graph.seqs.iter().map(Vec::len).sum::<usize>(),
        "pre_context_merge_stored_bases":merged.pan_bases,"coalesced_parallel_nodes":coalesced_parallel_nodes,
        "verified_input_segments":merged.verified_segments,"verified_input_links":merged.verified_links,"verified_input_paths":merged.verified_paths,
        "reference_only_bases_emitted":0,"output_scope":"union_of_certain_read_compatible_reference_intervals_and_de_novo_target_graph",
        "complete_assembly":false,"long_range_phase_identified":false,"physical_molecule_count":null,"physical_topology_identified":false}),
    )
}

/// Ends of graph segments with no link: where the graph stops, not a molecule's chosen end.
fn dead_ends(g: &Graph) -> usize {
    let mut linked = BTreeSet::new();
    for &(a, b) in &g.links {
        // leaving `a` uses its end on the strand it is read; entering `b` its start
        linked.insert((a.0, a.1));
        linked.insert((b.0, !b.1));
    }
    (0..g.names.len())
        .map(|i| usize::from(!linked.contains(&(i, true))) + usize::from(!linked.contains(&(i, false))))
        .sum()
}

fn read_fasta(path: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    for line in fs::read_to_string(path)?.lines() {
        if let Some(h) = line.strip_prefix('>') {
            out.push((h.to_string(), Vec::new()));
        } else if let Some(last) = out.last_mut() {
            last.1.extend_from_slice(line.trim().as_bytes());
        }
    }
    Ok(out)
}

/// The organelle GFA with `molecule` closed: the P walk that spells it gets `fill` (a new
/// segment, unless empty, with `segment_tags`) and the link back to its first step, and is
/// marked circular, with `cj:Z:<how>`. None when no P walk spells the molecule.
fn close_gfa(
    text: &str,
    molecule: &[u8],
    fill: &[u8],
    how: &str,
    segment_tags: &str,
) -> Result<Option<String>> {
    let mut seqs: BTreeMap<&str, &[u8]> = BTreeMap::new();
    for line in text.lines() {
        let f: Vec<_> = line.split('\t').collect();
        if f[0] == "S" && f.len() > 2 {
            seqs.insert(f[1], f[2].as_bytes());
        }
    }
    let spell = |steps: &str| -> Option<Vec<u8>> {
        let mut s = Vec::new();
        for t in steps.split(',') {
            let (name, o) = t.split_at(t.len().checked_sub(1)?);
            let seq = seqs.get(name)?;
            s.extend(if o == "+" { seq.to_vec() } else { revcomp(seq) });
        }
        Some(s)
    };
    let Some(walk) = text.lines().find(|l| {
        let f: Vec<_> = l.split('\t').collect();
        f[0] == "P" && f.len() > 2 && spell(f[2]).as_deref() == Some(molecule)
    }) else {
        return Ok(None);
    };
    let f: Vec<_> = walk.split('\t').collect();
    let steps: Vec<&str> = f[2].split(',').collect();
    fn split(t: &str) -> (&str, &str) {
        t.split_at(t.len() - 1)
    }
    let ((first, fo), (last, lo)) = (split(steps[0]), split(steps[steps.len() - 1]));
    let prefix = first.rsplit_once('.').map_or(first, |p| p.0);
    let mut name = format!("{prefix}.{how}_join");
    while seqs.contains_key(name.as_str()) {
        name.push('_');
    }
    let mut new_steps: Vec<String> = steps.iter().map(|s| s.to_string()).collect();
    let mut added = String::new();
    if fill.is_empty() {
        writeln!(added, "L\t{last}\t{lo}\t{first}\t{fo}\t0M")?;
    } else {
        writeln!(
            added,
            "S\t{name}\t{}\tLN:i:{}{segment_tags}\nL\t{last}\t{lo}\t{name}\t+\t0M\nL\t{name}\t+\t{first}\t{fo}\t0M",
            std::str::from_utf8(fill)?,
            fill.len()
        )?;
        new_steps.push(format!("{name}+"));
    }
    let overlaps = if f.len() > 3 && f[3] == "*" {
        "*".to_string()
    } else {
        vec!["0M"; new_steps.len() - 1].join(",")
    };
    let mut closed = format!("P\t{}\t{}\t{overlaps}", f[1], new_steps.join(","));
    for tag in f.iter().skip(4).filter(|t| !t.starts_with("TP:Z:")) {
        write!(closed, "\t{tag}")?;
    }
    write!(closed, "\tTP:Z:circular\tcj:Z:{how}")?;
    let mut out = String::new();
    for line in text.lines() {
        writeln!(out, "{}", if line == walk { closed.as_str() } else { line })?;
    }
    out.push_str(&added);
    Ok(Some(out))
}

/// Most covers `linearize` reports per component when every configuration is wanted.
const MAX_CONFIGURATIONS: usize = 10_000;

/// A walk of the assembly graph as a key equal for the same molecule however it is written: a
/// linear walk or its reverse, a circular one from any step in either direction.
fn walk_key(steps: &[(String, bool)], circular: bool) -> Vec<(String, bool)> {
    let rev: Vec<(String, bool)> = steps.iter().rev().map(|(n, o)| (n.clone(), !o)).collect();
    if !circular {
        return steps.to_vec().min(rev);
    }
    let mut best: Option<Vec<(String, bool)>> = None;
    for w in [steps.to_vec(), rev] {
        for i in 0..w.len() {
            let r: Vec<_> = w[i..].iter().chain(&w[..i]).cloned().collect();
            if best.as_ref().map_or(true, |b| r < *b) {
                best = Some(r);
            }
        }
    }
    best.unwrap_or_default()
}

/// Every configuration of the organelle's graph its reads allow, when they do not decide one:
/// each cover of the anchor graph `linearize` finds (every pairing of every repeat's copies the
/// copy numbers allow), as the molecules it makes. Short reads cannot say which: a plant
/// mitochondrion's repeat longer than the fragments recombines, and its copies' pairings are
/// all as well supported (ZM4: 32 covers of five repeats, among them the held-out HiFi
/// assembly's two molecules). Written to `<organelle>.configurations.fasta` and as P walks of
/// `<organelle>.gfa` (`<sample>_configuration<N>#1#<n>`), all labelled unresolved; the
/// representative molecule stays one of them.
fn configurations(
    dir: &Path,
    out: &Path,
    organelle: &str,
    gfa_out: &Path,
    ref_walks: &[(Vec<u8>, bool)],
    ref_names: &[String],
) -> Result<(Value, Option<Representative>)> {
    let d = dir.join("de-novo");
    let rep = crate::linearize::run(
        &d.join("assembly.gfa"),
        &d.join("evidence.json"),
        crate::linearize::LinearizeParams {
            min_reads: 2,
            max_steps: 5_000_000,
            max_molecules: 0,
            max_alternatives: MAX_CONFIGURATIONS,
            optional_below: 0.5,
        },
        &d.join("configurations.best.fasta"),
        Some(&d.join("configurations.partial.fasta")),
        &d.join("configurations.json"),
    )?;
    let repeats: BTreeSet<String> = rep.repeats.iter().cloned().collect();
    // the component holding most of the organelle
    let Some(component) = rep.components.iter().max_by_key(|c| {
        c.best.as_ref().map_or(0, |b| b.molecules.iter().map(|m| m.length).sum::<usize>())
    }) else {
        return Ok((json!({"configurations": 0, "reason": "no component"}), None));
    };
    let covers: Vec<_> = component.best.iter().chain(&component.alternatives).collect();
    // segments of the organelle graph by the assembly segment they came from, in order
    let text = fs::read_to_string(gfa_out)?;
    let sample = text
        .lines()
        .next()
        .and_then(|h| h.split('\t').find_map(|t| t.strip_prefix("sn:Z:")))
        .unwrap_or("")
        .to_string();
    // The molecules are spelled on the assembly graph, each link's overlap taken off once. The
    // organelle graph is blunt: a segment with an overlapping link loses those bases for every
    // walk, and ZM4's 360 bp repeat (126 bp overlaps to two short tips) came out 126 bp short in
    // walks that pass it without the tips (145,947 bp against the HiFi molecule's 146,073).
    let asm = parse_gfa(&d.join("assembly.gfa"))?;
    let index: BTreeMap<&str, usize> = asm.names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let mut overlap: BTreeMap<Edge, usize> = BTreeMap::new();
    for (&(a, b), &o) in asm.links.iter().zip(&asm.overlaps) {
        overlap.insert((a, b), o);
        overlap.insert((flip(b), flip(a)), o);
    }
    // The organelle graph is the assembly graph made blunt by `unify`: a walk over assembly
    // segments is one over its pieces (`unify::piece_walk`), named as unify names them. A
    // repeat's bases shared with an overlapping neighbour are a piece owned by the neighbour:
    // ZM4's 360 bp repeat is `u13.1`, `u13.2` and `u12.2`.
    let blunt = crate::unify::blunt_with(&asm, |_, _| Ok(()))?;
    let mut owned = vec![0usize; asm.seqs.len()];
    let mut rank = Vec::with_capacity(blunt.pieces.len());
    for pc in &blunt.pieces {
        owned[pc.owner] += 1;
        rank.push(owned[pc.owner]);
    }
    let piece_name = |piece: usize| -> String {
        let o = blunt.pieces[piece].owner;
        if owned[o] == 1 {
            format!("{sample}.{}", asm.names[o])
        } else {
            format!("{sample}.{}.{}", asm.names[o], rank[piece])
        }
    };
    let spell = |steps: &[(String, bool)], circular: bool| -> Option<Vec<u8>> {
        let walk: Vec<Oriented> = steps
            .iter()
            .map(|(n, o)| index.get(n.as_str()).map(|&i| (i, *o)))
            .collect::<Option<_>>()?;
        let mut seq = asm.oriented_seq(walk[0]);
        for w in walk.windows(2) {
            let o = overlap.get(&(w[0], w[1])).copied().unwrap_or(0);
            seq.extend_from_slice(&asm.oriented_seq(w[1])[o..]);
        }
        if circular {
            let o = overlap.get(&(walk[walk.len() - 1], walk[0])).copied().unwrap_or(0);
            seq.truncate(seq.len() - o);
        }
        Some(seq)
    };
    let short = if organelle == "plastid" { "pt" } else { "mt" };
    let (mut seen, mut fasta, mut walks, mut listed) = (BTreeSet::new(), String::new(), String::new(), Vec::new());
    let mut skipped = 0usize;
    // (reads' score, molecules: (linearize molecule, steps, organelle-graph steps, sequence,
    // what the GFA walk spells))
    type Molecule<'a> = (&'a crate::linearize::Linearization, Vec<(String, bool)>, Vec<String>, Vec<u8>, usize);
    let mut found: Vec<(f64, Vec<Molecule>)> = Vec::new();
    'covers: for cover in covers {
        let mut key = Vec::new();
        let mut molecules = Vec::new();
        for m in &cover.molecules {
            let steps: Vec<(String, bool)> = m
                .path
                .split_whitespace()
                .map(|t| (t[..t.len() - 1].to_string(), t.ends_with('+')))
                .collect();
            key.push(walk_key(&steps, m.circular));
            // the organelle graph's steps (and what they spell) for it
            let walk: Option<Vec<Oriented>> = steps
                .iter()
                .map(|(n, o)| index.get(n.as_str()).map(|&i| (i, *o)))
                .collect();
            let Some(pw) = walk.and_then(|w| crate::unify::piece_walk(&asm, &blunt, &w, m.circular).ok()) else {
                skipped += 1;
                continue 'covers;
            };
            let gsteps: Vec<String> = pw
                .iter()
                .map(|&(p, f)| format!("{}{}", piece_name(p), if f { '+' } else { '-' }))
                .collect();
            let gseq: usize = pw.iter().map(|&(p, _)| blunt.piece_len(p)).sum();
            let Some(seq) = spell(&steps, m.circular) else {
                skipped += 1;
                continue 'covers;
            };
            molecules.push((m, steps, gsteps, seq, gseq));
        }
        key.sort();
        if seen.insert(key) {
            found.push((cover.score, molecules));
        }
    }
    // The reads leave the pairings open; the reference runs' assemblies show which pairings
    // individuals of the species have. Each configuration is ranked by how often the runs hold
    // the pairings it makes (`reference_prior`), first the one most like them.
    let pairings: Vec<Vec<Pairing>> = found
        .iter()
        .map(|(_, ms)| ms.iter().flat_map(|m| pairings_of(&m.1, m.0.circular, &repeats)).collect())
        .collect();
    let prior = reference_prior(&asm, &index, &overlap, &pairings, ref_walks, ref_names);
    let mut order: Vec<usize> = (0..found.len()).collect();
    order.sort_by(|&a, &b| {
        prior.scores[b]
            .total_cmp(&prior.scores[a])
            .then(found[b].0.total_cmp(&found[a].0))
            .then(a.cmp(&b))
    });
    // the circular representative: of the configurations whose molecules of at least
    // `pipeline::MIN_PIECE` are all circular, the one with fewest of them (a master circle when
    // there is one), the best ranked of those
    let mut circular: Option<Representative> = None;
    for (rank, &c) in order.iter().enumerate() {
        let (score, molecules) = &found[c];
        let n = rank + 1;
        let big: Vec<_> = molecules.iter().filter(|m| m.3.len() as u64 >= pipeline::MIN_PIECE).collect();
        if !big.is_empty()
            && big.iter().all(|m| m.0.circular)
            && circular.as_ref().map_or(true, |r| big.len() < r.pieces)
        {
            circular = Some(Representative {
                configuration: n,
                pieces: big.len(),
                molecules: molecules
                    .iter()
                    .enumerate()
                    .map(|(j, m)| (format!("{sample}_configuration{n}#1#{short}.{}", j + 1), m.3.clone(), m.0.circular))
                    .collect(),
            });
        }
        let runs = &prior.runs_with_all[c];
        let mut lengths = Vec::new();
        for (j, (m, _, gsteps, seq, gseq)) in molecules.iter().enumerate() {
            let name = format!("{sample}_configuration{n}#1#{short}.{}", j + 1);
            writeln!(
                fasta,
                ">{name} configuration={n} circular={} length={} path={} reads_decide=false reference_prior_ln={:.3} reference_runs_with_all_pairings={}",
                m.circular,
                seq.len(),
                m.path.replace(' ', ","),
                prior.scores[c],
                if runs.is_empty() { "none".to_string() } else { runs.join(",") }
            )?;
            for chunk in seq.chunks(80) {
                writeln!(fasta, "{}", std::str::from_utf8(chunk)?)?;
            }
            writeln!(
                walks,
                "P\t{name}\t{}\t*{}\tcf:i:{n}\tst:Z:unresolved_configuration\trp:f:{:.3}",
                gsteps.join(","),
                if m.circular { "\tTP:Z:circular" } else { "" },
                prior.scores[c]
            )?;
            // the GFA walk spells the blunt segments: shorter where a segment lost bases to an
            // overlapping link this molecule does not use
            lengths.push(json!({"length": seq.len(), "circular": m.circular, "gfa_walk_bp": gseq}));
        }
        let detail: Vec<Value> = pairings[c]
            .iter()
            .map(|p| {
                let (seen_in, informative) = prior.observed.get(p).cloned().unwrap_or_default();
                json!({"repeat": p.1 .0, "from": format!("{}{}", p.0 .0, if p.0 .1 {'+'} else {'-'}),
                    "to": format!("{}{}", p.2 .0, if p.2 .1 {'+'} else {'-'}),
                    "runs_with_it": seen_in.len(), "runs_with_its_flanks": informative.len(), "runs": seen_in})
            })
            .collect();
        listed.push(json!({"configuration": n, "reads_score_ln": score, "reference_prior_ln": prior.scores[c],
            "reference_runs_with_all_pairings": runs, "pairings": detail, "molecules": lengths}));
    }
    let path = out.join(format!("{organelle}.configurations.fasta"));
    fs::write(&path, fasta)?;
    fs::write(gfa_out, text + &walks)?;
    let best = found.iter().map(|f| f.0).fold(f64::NEG_INFINITY, f64::max);
    let tied = found.iter().filter(|f| best - f.0 < 1e-9).count();
    eprintln!(
        "[ovasm] reconstruct {organelle}: {} configurations ({tied} as well supported by the reads as the best), ranked by {} reference runs",
        listed.len(),
        prior.runs
    );
    Ok((
        json!({"configurations": listed.len(), "as_well_supported_by_reads_as_the_best": tied,
        "covers_searched": component.candidates, "search_truncated": component.search_truncated,
        "cover_size_limited": component.cover_size_limited, "skipped_covers": skipped,
        "reference_runs": prior.runs, "fasta": path, "list": listed,
        "circular_representative": circular.as_ref().map(|r| r.configuration),
        "basis": "every cover of the anchor graph linearize finds; the reads do not decide between \
                  them; ranked by how often the reference runs hold each repeat pairing (a prior \
                  from other individuals, not evidence of this one)"}),
        circular,
    ))
}

/// A configuration taken as the organelle's representative molecules (`configurations`):
/// its number, molecules of at least `pipeline::MIN_PIECE`, and every molecule as
/// (configuration walk name, sequence, circular).
struct Representative {
    configuration: usize,
    pieces: usize,
    molecules: Vec<(String, Vec<u8>, bool)>,
}

/// A repeat's copy entered from one step and left to another: (from, repeat, to), as the
/// lesser of it and its reverse.
type Pairing = ((String, bool), (String, bool), (String, bool));

fn pairings_of(steps: &[(String, bool)], circular: bool, repeats: &BTreeSet<String>) -> Vec<Pairing> {
    let n = steps.len();
    let flip = |s: &(String, bool)| (s.0.clone(), !s.1);
    (0..n)
        .filter(|&i| repeats.contains(&steps[i].0))
        .filter_map(|i| {
            let prev = if i > 0 { Some(&steps[i - 1]) } else { circular.then(|| &steps[n - 1]) }?;
            let next = if i + 1 < n { Some(&steps[i + 1]) } else { circular.then(|| &steps[0]) }?;
            let fwd = (prev.clone(), steps[i].clone(), next.clone());
            let rev = (flip(next), flip(&steps[i]), flip(prev));
            Some(fwd.min(rev))
        })
        .collect()
}

struct Prior {
    /// Reference runs (the first PanSN field of the walks).
    runs: usize,
    /// Each pairing: the runs that hold it, and the runs that hold its entry's flank at all.
    observed: BTreeMap<Pairing, (Vec<String>, Vec<String>)>,
    /// Each configuration: the sum of ln((runs with the pairing + 0.5) / (runs with its flank + 1))
    /// over its pairings.
    scores: Vec<f64>,
    /// Each configuration: the runs that hold every one of its pairings.
    runs_with_all: Vec<Vec<String>>,
}

/// Bases either side of a repeat that place a pairing on a reference walk.
const FLANK_K: usize = 31;

/// How often the reference walks hold each pairing: the entry's last `FLANK_K` bases followed,
/// at the distance the repeat's copy makes (within 5%, at least 20 bp: other individuals' copies differ
/// a little; 200 bp made both pairings of ZM4's 135 and 360 bp repeats "held by 8 of 8"), by the
/// exit's first `FLANK_K` bases, on either strand of any walk of a run.
fn reference_prior(
    asm: &Graph,
    index: &BTreeMap<&str, usize>,
    overlap: &BTreeMap<Edge, usize>,
    pairings: &[Vec<Pairing>],
    walks: &[(Vec<u8>, bool)],
    names: &[String],
) -> Prior {
    let run_of = |w: usize| names[w].split('#').next().unwrap_or("").to_string();
    let all_runs: BTreeSet<String> = (0..walks.len()).map(run_of).collect();
    let oriented = |s: &(String, bool)| index.get(s.0.as_str()).map(|&i| (i, s.1));
    // the probes: an entry's last bases and an exit's first, and the distance between them
    let distinct: BTreeSet<&Pairing> = pairings.iter().flatten().collect();
    let mut probes: BTreeMap<Vec<u8>, Vec<usize>> = BTreeMap::new();
    let mut plans: Vec<(&Pairing, usize, usize, usize, usize)> = Vec::new();
    for &p in &distinct {
        let (Some(a), Some(r), Some(b)) = (oriented(&p.0), oriented(&p.1), oriented(&p.2)) else {
            continue;
        };
        let (sa, sb) = (asm.oriented_seq(a), asm.oriented_seq(b));
        if sa.len() < FLANK_K || sb.len() < FLANK_K {
            continue;
        }
        let gap = asm.seqs[r.0].len() as i64
            - *overlap.get(&(a, r)).unwrap_or(&0) as i64
            - *overlap.get(&(r, b)).unwrap_or(&0) as i64;
        let (ia, ib) = (probes.len(), probes.len() + 1);
        probes.entry(sa[sa.len() - FLANK_K..].to_ascii_uppercase()).or_default().push(ia);
        probes.entry(sb[..FLANK_K].to_ascii_uppercase()).or_default().push(ib);
        plans.push((p, ia, ib, gap.max(0) as usize, asm.seqs[r.0].len()));
    }
    let longest = plans.iter().map(|p| p.3).max().unwrap_or(0) + 2 * FLANK_K + 400;
    // where each probe lies: (walk, reverse strand?, position)
    let mut at: BTreeMap<usize, Vec<(usize, bool, usize)>> = BTreeMap::new();
    for (w, (seq, circular)) in walks.iter().enumerate() {
        for reverse in [false, true] {
            let mut s = if reverse { revcomp(seq) } else { seq.to_ascii_uppercase() };
            if *circular {
                let wrap: Vec<u8> = s.iter().cycle().take(longest.min(s.len())).copied().collect();
                s.extend(wrap);
            }
            for (pos, x) in s.windows(FLANK_K).enumerate() {
                for &id in probes.get(x).map(Vec::as_slice).unwrap_or(&[]) {
                    at.entry(id).or_default().push((w, reverse, pos));
                }
            }
        }
    }
    let mut observed: BTreeMap<Pairing, (Vec<String>, Vec<String>)> = BTreeMap::new();
    for &(p, ia, ib, gap, len) in &plans {
        let tol = 20.max(len / 20) as i64;
        let (ea, eb) = (at.get(&ia).cloned().unwrap_or_default(), at.get(&ib).cloned().unwrap_or_default());
        let mut with: BTreeSet<String> = BTreeSet::new();
        for &(w, rv, pa) in &ea {
            if eb.iter().any(|&(w2, rv2, pb)| {
                w2 == w && rv2 == rv && ((pb as i64 - (pa + FLANK_K) as i64) - gap as i64).abs() <= tol
            }) {
                with.insert(run_of(w));
            }
        }
        let flanked: BTreeSet<String> = ea.iter().map(|&(w, _, _)| run_of(w)).collect();
        observed.insert(p.clone(), (with.into_iter().collect(), flanked.into_iter().collect()));
    }
    let mut scores = Vec::new();
    let mut runs_with_all = Vec::new();
    for ps in pairings {
        let mut s = 0.0;
        let mut all: BTreeSet<String> = all_runs.clone();
        for p in ps {
            let (with, flanked) = observed.get(p).cloned().unwrap_or_default();
            s += ((with.len() as f64 + 0.5) / (flanked.len() as f64 + 1.0)).ln();
            let with: BTreeSet<String> = with.into_iter().collect();
            all = all.intersection(&with).cloned().collect();
        }
        scores.push(s);
        runs_with_all.push(if ps.is_empty() { Vec::new() } else { all.into_iter().collect() });
    }
    Prior {
        runs: all_runs.len(),
        observed,
        scores,
        runs_with_all,
    }
}

/// The organelle's result as `ovasm run` gives it, at the top of the output: `<organelle>.gfa`
/// (the de novo organelle graph with its molecule paths) and `<organelle>.fasta` (its
/// molecules). A single linear molecule of the genome is first closed with read pairs (`close`)
/// when reads spell its join. It is complete as one circular molecule of at least
/// `pipeline::MIN_PIECE` from a graph without loose segment ends; otherwise every reason is
/// named. Whether the reads decide every repeat pairing is reported beside it (`decisive`):
/// short reads cannot decide a repeat longer than their fragments.
/// `individual.gfa` stays the graph for a species pan graph (reference intervals and de novo).
fn final_outputs(
    dir: &Path,
    out: &Path,
    organelle: &str,
    de_novo: &Value,
    reads: &[PathBuf],
    walks: &[(Vec<u8>, bool)],
    walk_names: &[String],
) -> Result<Value> {
    let graph = dir.join("de-novo/organelle.gfa");
    let fasta = dir.join("de-novo/molecules.fasta");
    let (gfa_out, fasta_out) = (
        out.join(format!("{organelle}.gfa")),
        out.join(format!("{organelle}.fasta")),
    );
    let mut gfa_text = fs::read_to_string(&graph)?;
    let mut molecules = if fasta.exists() {
        read_fasta(&fasta)?
    } else {
        Vec::new()
    };
    let mut circular: Vec<bool> = molecules
        .iter()
        .map(|(h, _)| h.contains("circular=true"))
        .collect();
    let pieces: Vec<usize> = (0..molecules.len())
        .filter(|&i| molecules[i].1.len() as u64 >= pipeline::MIN_PIECE)
        .collect();
    let mut closure = Value::Null;
    let mut reference_fill = Value::Null;
    let mut reference_bp = 0usize;
    if let [i] = pieces[..] {
        if !circular[i] {
            eprintln!("[ovasm] reconstruct {organelle}: closing the linear molecule with read pairs");
            let mut c = close::close(&molecules[i].1, reads)?;
            // the reads spell the join; else the pairs join the ends and a reference walk holds
            // sequence of their distance between them
            let mut join: Option<(Vec<u8>, &str, String, String)> = None;
            if c.fill.is_some() {
                join = Some((
                    c.sequence.clone(),
                    "read_pairs",
                    String::new(),
                    format!("closed_by=read_pairs join_bp={} joining_pairs={}", c.sequence.len(), c.joining_pairs),
                ));
            } else if let (true, Some(gap)) = (c.joining_pairs >= close::MIN_JOINING_PAIRS, c.gap_estimate) {
                if let Some(f) = close::reference_fill(&molecules[i].1, walks, gap) {
                    let source = &walk_names[f.walk];
                    join = Some((
                        f.sequence.clone(),
                        "reference_fill",
                        format!("\trf:Z:{source}"),
                        format!(
                            "closed_by=reference_fill reference_fill_bp={} reference_walk={source} reference_similarity={:.3} joining_pairs={} pair_gap_estimate={gap:.0}",
                            f.length, f.similarity, c.joining_pairs
                        ),
                    ));
                    reference_fill = json!({"walk": source, "reverse": f.reverse, "length": f.length,
                        "similarity": f.similarity, "candidates": f.candidates,
                        "runner_up": f.runner_up.map(|(w, s)| json!({"walk": walk_names[w], "similarity": s})),
                        "pair_gap_estimate": gap,
                        "basis": "the agreeing reference walk most like the individual (share of its sampled 31-mers); reference bases, not read evidence"});
                }
            }
            if let Some((fill, how, tags, note)) = join {
                match close_gfa(&gfa_text, &molecules[i].1, &fill, how, &tags)? {
                    Some(t) => {
                        gfa_text = t;
                        if how == "reference_fill" {
                            reference_bp = fill.len();
                        }
                        let (h, s) = &mut molecules[i];
                        s.extend_from_slice(&fill);
                        let renamed = h
                            .split(' ')
                            .map(|w| {
                                if w == "circular=false" {
                                    "circular=true".to_string()
                                } else if w.starts_with("length=") {
                                    format!("length={}", s.len())
                                } else {
                                    w.to_string()
                                }
                            })
                            .collect::<Vec<_>>()
                            .join(" ");
                        *h = format!("{renamed} {note}");
                        circular[i] = true;
                    }
                    None => {
                        c.fill = None;
                        reference_fill = Value::Null;
                        c.stopped = Some("no P walk of the organelle graph spells the molecule".into());
                    }
                }
            }
            closure = serde_json::to_value(&c)?;
        }
    }
    fs::write(&gfa_out, &gfa_text)?;
    if !molecules.is_empty() {
        let mut text = String::new();
        for (h, s) in &molecules {
            writeln!(text, ">{h}")?;
            for chunk in s.chunks(80) {
                writeln!(text, "{}", std::str::from_utf8(chunk)?)?;
            }
        }
        fs::write(&fasta_out, text)?;
    }
    let decisive = de_novo["decisive"].as_bool().unwrap_or(false);
    let loose = dead_ends(&parse_gfa(&gfa_out)?);
    let mut reasons = Vec::new();
    if pieces.len() != 1 {
        reasons.push(format!(
            "{} molecules of at least {} bp",
            pieces.len(),
            pipeline::MIN_PIECE
        ));
    }
    if pieces.iter().any(|&i| !circular[i]) {
        reasons.push("the molecule is linear".to_string());
    }
    if loose > 0 {
        reasons.push(format!("{loose} loose segment ends in the organelle graph"));
    }
    let status = match (reasons.is_empty(), reference_bp) {
        (false, _) => "incomplete",
        (true, 0) => "complete",
        // closed, but with bases a reference walk holds, not the reads
        (true, _) => "complete_with_reference_fill",
    };
    // the reads leave the configuration open: every one the graph allows
    let (configurations_found, chosen) = if status == "incomplete" && !decisive && !pieces.is_empty() {
        configurations(dir, out, organelle, &gfa_out, walks, walk_names)?
    } else {
        (Value::Null, None)
    };
    let mut status = status.to_string();
    let mut representative = Value::Null;
    let (mut lengths, mut circular) = (molecules.iter().map(|m| m.1.len()).collect::<Vec<_>>(), circular);
    // A configuration of circular molecules, when the reads leave one open, is the representative:
    // as `run` gives a circle for a molecule of the graph, not the linear walk that ends on short
    // tips beside a repeat (ZM4: two dead-end tips of 177 and 168 bp at its 360 bp repeat).
    if let Some(r) = &chosen {
        let total = configurations_found["configurations"].as_u64().unwrap_or(0);
        let mut text = String::new();
        for (name, seq, circ) in &r.molecules {
            let unplaced = (seq.len() as u64) < pipeline::MIN_PIECE;
            writeln!(
                text,
                ">{name} circular={circ} length={} configuration={} of={total} chosen_by=fewest_circular_molecules reads_decide=false{}",
                seq.len(),
                r.configuration,
                if unplaced { " unplaced=true" } else { "" }
            )?;
            for chunk in seq.chunks(80) {
                writeln!(text, "{}", std::str::from_utf8(chunk)?)?;
            }
        }
        fs::write(&fasta_out, text)?;
        lengths = r.molecules.iter().map(|m| m.1.len()).collect();
        circular = r.molecules.iter().map(|m| m.2).collect();
        let small = r.molecules.iter().filter(|m| (m.1.len() as u64) < pipeline::MIN_PIECE).count();
        reasons = vec![format!(
            "the reads leave the repeat pairings open: configuration {} of {total}, the fewest circular molecules ({}), is the representative; all are in {organelle}.configurations.fasta",
            r.configuration, r.pieces
        )];
        if small > 0 {
            reasons.push(format!("{small} piece(s) under {} bp not placed: the graph's short dead-end tips", pipeline::MIN_PIECE));
        }
        status = "circular_configuration_chosen".to_string();
        representative = json!({"configuration": r.configuration, "circular_molecules": r.pieces,
            "chosen_by": "fewest circular molecules, then the reference runs' ranking"});
    }
    // stretches whose length the reads do not fix (collapsed tandem arrays, microsatellites)
    let uncertain = if fasta_out.exists() {
        let final_molecules: Vec<(String, Vec<u8>, bool)> = read_fasta(&fasta_out)?
            .into_iter()
            .map(|(h, s)| (h.split(' ').next().unwrap_or("").to_string(), s, h.contains("circular=true")))
            .collect();
        let u = crate::uncertain::run(&final_molecules, reads, pipeline::MIN_PIECE as usize, dir)?;
        eprintln!("[ovasm] reconstruct {organelle}: {} stretch(es) of uncertain length", u["marked"]);
        u
    } else {
        Value::Null
    };
    Ok(json!({"gfa": gfa_out, "fasta": (!molecules.is_empty()).then_some(fasta_out),
        "molecules": lengths.len(), "lengths": lengths, "circular": circular,
        "pieces_of_at_least_1kb": pieces.len(), "closure": closure,
        "reference_fill": reference_fill, "reference_fill_bp": reference_bp,
        "configurations": configurations_found, "representative": representative,
        "uncertain_regions": uncertain,
        "graph_loose_ends": loose, "status": status,
        "reasons": reasons, "decisive": decisive,
        "repeats": if decisive {"every pairing decided by reads"} else {
            "some repeat pairings not decided by reads (a repeat longer than the fragments, or too few spanning reads): the molecule is one configuration the graph allows"}}))
}

/// The individual written back into the species graph: the reference graph and the
/// individual's organelle graph (`<organelle>.gfa`, its segments and links with only its
/// representative molecules' walks, PanSN haplotype 1, renamed to the sample; neither the
/// unresolved configurations nor other conformations) merged by `pan` into
/// `<organelle>.pan.gfa`. Bases the reference holds are shared; the rest is the individual's
/// new sequence, added as segments of its own.
fn write_back(
    reference: &Path,
    organelle_gfa: &Path,
    sample: &str,
    organelle: &str,
    configuration: Option<u64>,
    out: &Path,
    dir: &Path,
) -> Result<Value> {
    let text = fs::read_to_string(organelle_gfa)?;
    // the representative's walks: a chosen configuration's, else haplotype 1's
    let from = match configuration {
        Some(n) => format!("{sample}_denovo_configuration{n}#"),
        None => format!("{sample}_denovo#"),
    };
    let mut kept = format!("H\tVN:Z:1.0\tov:Z:organelleverse.graph.v1\tsn:Z:{sample}\n");
    let mut paths = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "S" | "L" => writeln!(kept, "{line}")?,
            "P" if f.len() > 2 => {
                if !f[1].starts_with(&from) || f[1].split('#').nth(1) != Some("1") {
                    continue;
                }
                let name = f[1].replacen(&from, &format!("{sample}#"), 1);
                writeln!(kept, "P\t{name}\t{}", f[2..].join("\t"))?;
                paths.push(name);
            }
            _ => {}
        }
    }
    let input = dir.join("individual.for-pan.gfa");
    fs::write(&input, kept)?;
    let pan_gfa = out.join(format!("{organelle}.pan.gfa"));
    let rep = pan::run(
        &[reference.to_path_buf(), input],
        &pan::PanParams { k: 63 },
        &pan_gfa,
        &dir.join("pan.json"),
    )?;
    let reference_bases: usize = parse_gfa(reference)?.seqs.iter().map(Vec::len).sum();
    Ok(json!({"pan_gfa": pan_gfa, "individual_paths": paths, "pan_bases": rep.pan_bases,
        "reference_bases": reference_bases,
        "individual_new_bases": rep.pan_bases.saturating_sub(reference_bases),
        "pieces": rep.pieces, "links": rep.links,
        "basis": "pan merge (k 63) of the reference graph and the individual's representative molecules"}))
}

pub struct Config {
    pub reads: Vec<PathBuf>,
    pub target_reads: BTreeMap<String, Vec<PathBuf>>,
    pub references: BTreeMap<String, PathBuf>,
    pub sample: String,
    pub out: PathBuf,
    pub target_depth: Option<f64>,
    pub salt: u64,
    /// Walk the individual's sequence that no reference path holds before recruitment
    /// (`frontier`); WGS input only.
    pub frontier: Option<frontier::Params>,
}

pub fn run(cfg: &Config) -> Result<Value> {
    if cfg.sample.is_empty() || cfg.sample.contains(['#', '\t', '\n', ' ']) {
        bail!("sample must be a nonempty PanSN field");
    }
    let targeted = !cfg.target_reads.is_empty();
    if targeted {
        if !cfg.reads.is_empty()
            || cfg.target_depth.is_some()
            || cfg.target_reads.len() != 2
            || ["mitochondrion", "plastid"]
                .iter()
                .any(|&o| cfg.target_reads.get(o).map_or(true, |p| p.is_empty()))
        {
            bail!("target mode needs both labelled organelle read sets, no raw WGS input or depth sampling");
        }
    } else if cfg.reads.is_empty() {
        bail!("give target WGS reads");
    }
    if let Some(d) = cfg.target_depth {
        if !d.is_finite() || d <= 0.0 {
            bail!("target depth must be positive and finite");
        }
    }
    if cfg.out.exists() && fs::read_dir(&cfg.out)?.next().is_some() {
        bail!("output directory is nonempty; preserve the previous run");
    }
    let start = Instant::now();
    // Validate all references before creating output or scanning large inputs.
    let mut references = BTreeMap::new();
    for o in ["mitochondrion", "plastid"] {
        references.insert(
            o,
            Reference::read(
                cfg.references
                    .get(o)
                    .with_context(|| format!("missing {o} reference graph"))?,
            )?,
        );
    }
    fs::create_dir_all(cfg.out.join("reference"))?;
    let mut bait = BaitSet::new(K);
    let mut reference_info = BTreeMap::new();
    for (&o, r) in &references {
        reference_info.insert(
            o,
            r.add_bait(
                &mut bait,
                o,
                &cfg.out.join("reference").join(format!("{o}.fasta")),
            )?,
        );
    }
    fs::write(
        cfg.out.join("reference.json"),
        serde_json::to_string_pretty(&reference_info)?,
    )?;
    // The individual's sequence that no reference path holds joins the bait: recruitment alone
    // reaches only a fifth of a read past the reference each pass (`frontier`).
    let frontier_info = match (&cfg.frontier, targeted) {
        (Some(p), false) => {
            eprintln!(
                "[ovasm] reconstruct: counting {}-mers of the WGS reads ({} GiB table)",
                frontier::K,
                (2u64 << p.slots_log2) >> 30
            );
            let counts = frontier::KmerCounts::from_reads(&cfg.reads, p.slots_log2)?;
            let mut targets = BTreeMap::new();
            for (&o, r) in &references {
                let sequences: Vec<Vec<u8>> = r
                    .bait_sequences(frontier::K)
                    .into_iter()
                    .map(|(_, s)| s)
                    .collect();
                let ext = frontier::extend(&counts, &sequences, r.stored_bases(), p);
                let t = bait.target(o)?;
                for s in &ext.sequences {
                    bait.add_sequence(t, s);
                }
                frontier::write_fasta(
                    &cfg.out.join("reference").join(format!("{o}.frontier.fasta")),
                    o,
                    &ext,
                )?;
                eprintln!(
                    "[ovasm] reconstruct {o}: depth {:.0}, {} ends, {} extended, {} new {}-mers in {} pieces{}",
                    ext.depth,
                    ext.ends,
                    ext.ends_extended,
                    ext.new_kmers,
                    frontier::K,
                    ext.unitigs,
                    if ext.budget_reached { " (budget reached)" } else { "" }
                );
                targets.insert(o, ext);
            }
            json!({"status":"run","params":p,"counted_reads":counts.reads,"counted_bases":counts.bases,
                "counting_seconds":counts.seconds,"table_bytes":2u64 << p.slots_log2,"targets":targets,
                "role":"bait_only_never_emitted_as_sequence"})
        }
        (Some(_), true) => json!({"status":"not_run_explicit_target_inputs"}),
        (None, _) => json!({"status":"disabled"}),
    };
    let recruitment = if targeted {
        json!({"status":"not_run_explicit_target_inputs","inputs":cfg.target_reads,
            "upstream_recruitment_certified":false})
    } else {
        serde_json::to_value(recruit::run(
            bait,
            recruit::RecruitConfig {
                replay_dir: None,
                inputs: cfg.reads.clone(),
                params: recruit::ClassifyParams {
                    min_distinct: 0.0,
                    min_cover: 0.8,
                    min_hit_frac: 0.2,
                    numt_min_block: usize::MAX,
                    block_gap: 40,
                },
                iterations: 20,
                extend_min_count: 3,
                extend_depth: 200.0,
                saturation: 0.05,
                target_depth: cfg.target_depth,
                salt: cfg.salt,
                out_dir: cfg.out.join("recruit"),
                bootstrap: false,
                extend_gate: false,
                max_read_depth: None,
                max_bases: None,
                // read pairs close the circles (`close`)
                mates_together: true,
                kmer_gate: false,
            },
        )?)?
    };
    let mut organelles = BTreeMap::new();
    for (&o, r) in &references {
        let dir = cfg.out.join(o);
        fs::create_dir_all(&dir)?;
        let reads = if targeted {
            cfg.target_reads[o].clone()
        } else {
            vec![cfg.out.join("recruit").join(format!("{o}.fastq"))]
        };
        eprintln!("[ovasm] reconstruct {o}: projecting exact target read contexts on all {} reference paths",r.walks.len());
        let projection = project(r, &reads, &cfg.sample, &dir)?;
        // the individual's allele at every known site of the reference graph
        eprintln!("[ovasm] reconstruct {o}: genotyping the reference graph's sites");
        let walk_steps: Vec<(String, Vec<Oriented>)> =
            r.walks.iter().map(|w| (w.name.clone(), w.steps.clone())).collect();
        let genotypes = crate::genotype::run(&r.graph, &walk_steps, &reads, &dir)?;
        eprintln!(
            "[ovasm] reconstruct {o}: {} sites, {} called, {} heteroplasmic; segments {}",
            genotypes["sites"], genotypes["called"], genotypes["heteroplasmic_sites"], genotypes["segments"]
        );
        let assembly = pipeline::run(&pipeline::RunConfig {
            reads: reads.clone(),
            organelle: if o == "mitochondrion" {
                pipeline::OrganelleChoice::Mitochondrion
            } else {
                pipeline::OrganelleChoice::Plastid
            },
            out: dir.join("de-novo"),
            read_type: pipeline::ReadType::Sr,
            read_set: pipeline::ReadSet::Target,
            seeds: BTreeMap::new(),
            builtin_seeds: BTreeSet::new(),
            recruit: pipeline::RecruitBy::Seeds,
            recruit_depth: None,
            recruit_cache: None,
            auto_recruit_cache: false,
            gated_recruit: false,
            recruit_max_depth: None,
            recruit_max_bases: None,
            recruit_iterations: None,
            read_context: None,
            single_seed: false,
            mode: pipeline::Mode::Fast,
            replicates: 3,
            max_replicates: 5,
            replicate_depth: None,
            salt: cfg.salt,
            sample: format!("{}_denovo", cfg.sample),
            component_panel: None,
        });
        let denovo_path = dir.join("de-novo/organelle.gfa");
        let integration = export_individual(
            &dir.join("observed-reference.gfa"),
            denovo_path.exists().then_some(denovo_path.as_path()),
            &cfg.sample,
            &dir,
        )?;
        let de_novo = match assembly {
            Ok(s) => s,
            Err(e) => json!({"error":format!("{e:#}")}),
        };
        let final_output = if de_novo.get("error").is_none() {
            let walks: Vec<(Vec<u8>, bool)> =
                r.walks.iter().map(|w| (w.sequence.clone(), w.circular)).collect();
            let names: Vec<String> = r.walks.iter().map(|w| w.name.clone()).collect();
            final_outputs(&dir, &cfg.out, o, &de_novo, &reads, &walks, &names)?
        } else {
            json!({"status":"failed","reasons":["de novo assembly failed"]})
        };
        let pan = if de_novo.get("error").is_none() {
            let chosen = final_output["representative"]["configuration"].as_u64();
            let p = write_back(&r.input, &cfg.out.join(format!("{o}.gfa")), &cfg.sample, o, chosen, &cfg.out, &dir)?;
            eprintln!(
                "[ovasm] reconstruct {o}: written back into the reference graph, {} new bases",
                p["individual_new_bases"]
            );
            p
        } else {
            Value::Null
        };
        eprintln!(
            "[ovasm] reconstruct {o}: {} ({})",
            final_output["status"].as_str().unwrap_or("?"),
            final_output["reasons"]
                .as_array()
                .map(|r| r.iter().filter_map(Value::as_str).collect::<Vec<_>>().join("; "))
                .unwrap_or_default()
        );
        organelles.insert(
            o,
            json!({"reference":reference_info[o],"projection":projection,"individual":integration,
            "de_novo":de_novo,"final":final_output,"genotypes":genotypes,"pan":pan}),
        );
    }
    let workflow_completed = organelles
        .values()
        .all(|s| s["de_novo"].get("error").is_none());
    let result = json!({"tool":"ovasm reconstruct","format":"organelleverse.reference-reconstruction.v1","sample":cfg.sample,
        "workflow_completed":workflow_completed,
        "inputs":cfg.reads,"target_read_inputs":cfg.target_reads,
        "input_scope":if targeted {"explicit_target_reads"} else {"whole_genome_recruitment"},
        "reference_graphs":cfg.references,"frontier":frontier_info,"recruitment":recruitment,"organelles":organelles,
        "elapsed_seconds":start.elapsed().as_secs_f64(),"peak_rss_mb":recruit::peak_rss_mb(),
        "strict_score":null,"relaxed_score":null,"independent_truth_evaluated":false,
        "physical_molecule_count":null,"physical_topology_identified":false});
    fs::write(
        cfg.out.join("reconstruction.json"),
        serde_json::to_string_pretty(&result)?,
    )?;
    if !workflow_completed {
        bail!("de novo target assembly failed; reconstruction.json preserves errors and partial graph evidence");
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn reference(text: &str) -> Reference {
        let p = std::env::temp_dir().join(format!(
            "ovasm-reconstruct-{}-{:?}.gfa",
            std::process::id(),
            std::thread::current().id()
        ));
        fs::write(&p, text).unwrap();
        let r = Reference::read(&p).unwrap();
        fs::remove_file(p).unwrap();
        r
    }
    fn random(n: usize) -> String {
        let mut x = 811u64;
        (0..n)
            .map(|_| {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                b"ACGT"[(x % 4) as usize] as char
            })
            .collect()
    }
    #[test]
    fn threads_across_tiny_variant_nodes_and_reverse_complement() {
        let s = random(120);
        let g=format!("S\ta\t{}\nS\tb\t{}\nS\tc\t{}\nL\ta\t+\tb\t+\t0M\nL\tb\t+\tc\t+\t0M\nP\tone#1#pt\ta+,b+,c+\t*\n",&s[..60],&s[60..61],&s[61..]);
        let r = reference(&g);
        let m = Mapper::new(&r);
        let seq = &s.as_bytes()[20..100];
        let hits = m.map(seq);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits, m.map(&revcomp(seq)));
        assert_eq!(hits.first().unwrap().edges.len(), 2);
    }
    #[test]
    fn duplicate_source_paths_do_not_multiply_read_support() {
        let s = random(90);
        let r = reference(&format!(
            "S\ta\t{s}\nP\tone#1#pt\ta+\t*\nP\ttwo#1#pt\ta+\t*\n"
        ));
        assert_eq!(r.sources.len(), 2);
        assert_eq!(Mapper::new(&r).map(&s.as_bytes()[10..70]).len(), 1);
    }
    #[test]
    fn ambiguous_repeat_coordinates_are_not_selected() {
        let s = random(90);
        let r = reference(&format!(
            "S\ta\t{s}\nS\tb\t{s}\nP\tone#1#pt\ta+\t*\nP\ttwo#1#pt\tb+\t*\n"
        ));
        let hits = Mapper::new(&r).map(&s.as_bytes()[10..70]);
        assert_eq!(hits.len(), 2);
        assert!(intersect(hits.first().unwrap(), hits.last().unwrap())
            .intervals
            .is_empty());
    }
    #[test]
    fn circular_boundary_and_unsupported_reference_bases() {
        let s = random(90);
        let r = reference(&format!(
            "S\ta\t{s}\nL\ta\t+\ta\t+\t0M\nP\tone#1#pt\ta+\t*\tTP:Z:circular\n"
        ));
        let seq = [&s.as_bytes()[70..], &s.as_bytes()[..30]].concat();
        let hits = Mapper::new(&r).map(&seq);
        let fp = hits.first().unwrap();
        assert_eq!(fp.intervals, vec![(0, 0, 30), (0, 70, 90)]);
        assert_eq!(fp.edges.len(), 1);
        let mut c = vec![vec![0; 90]];
        count_ranges(&mut c, fp);
        assert_eq!(covered_ranges(&c[0]), vec![(0, 30), (70, 90)]);
    }
    #[test]
    fn second_source_allele_is_used_and_novel_sequence_is_unmatched() {
        let s = random(90);
        let mut variant = s.as_bytes().to_vec();
        variant[45] = if variant[45] == b'A' { b'C' } else { b'A' };
        let r = reference(&format!(
            "S\ta\t{s}\nS\tb\t{}\nP\tone#1#pt\ta+\t*\nP\ttwo#1#pt\tb+\t*\n",
            std::str::from_utf8(&variant).unwrap()
        ));
        let hits = Mapper::new(&r).map(&variant[15..75]);
        assert_eq!(hits.len(), 1);
        assert_eq!(hits.first().unwrap().intervals[0].0, 1);
        variant[35] = if variant[35] == b'A' { b'C' } else { b'A' };
        assert!(Mapper::new(&r).map(&variant[15..75]).is_empty());
    }
    #[test]
    fn reverse_match_with_central_rarest_anchor_is_not_suppressed() {
        let s = random(81);
        let r = reference(&format!("S\ta\t{s}\nP\tone#1#pt\ta+\t*\n"));
        let m = Mapper::new(&r);
        assert_eq!(m.map(&revcomp(s.as_bytes())).len(), 1);
    }
    #[test]
    fn projection_excludes_uncovered_dna_and_integration_preserves_novel_dna() {
        let s = random(160);
        let r = reference(&format!(
            "S\ta\t{s}\nP\tone#1#pt\ta+\t*\nP\ttwo#1#pt\ta+\t*\n"
        ));
        let dir = std::env::temp_dir().join(format!("ovasm-project-test-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let reads = dir.join("reads.fa");
        fs::write(&reads, format!(">r\n{}\n", &s[30..100])).unwrap();
        let report = project(&r, &[reads], "target", &dir).unwrap();
        assert_eq!(report["exact_matching_records"], 1);
        assert_eq!(report["certain_exact_reference_bases"], 70);
        let observed = parse_gfa(&dir.join("observed-reference.gfa")).unwrap();
        assert_eq!(observed.seqs, vec![s.as_bytes()[30..100].to_vec()]);
        let novel = dir.join("novel.gfa");
        let n = revcomp(&s.as_bytes()[100..]);
        fs::write(
            &novel,
            format!(
                "H\tVN:Z:1.0\tsn:Z:target_denovo\nS\tn\t{}\nP\ttarget_denovo#1#pt\tn+\t*\n",
                std::str::from_utf8(&n).unwrap()
            ),
        )
        .unwrap();
        let integrated = export_individual(
            &dir.join("observed-reference.gfa"),
            Some(&novel),
            "target",
            &dir,
        )
        .unwrap();
        assert_eq!(integrated["stored_bases"], 130);
        assert_eq!(integrated["verified_input_paths"], 1);
        assert_eq!(integrated["complete_assembly"], false);
        fs::remove_dir_all(&dir).unwrap();
    }
    #[test]
    fn loose_ends_are_segment_ends_without_a_link() {
        let chain = reference("S\ta\tACGT\nS\tb\tGG\nL\ta\t+\tb\t-\t0M\nP\tone#1#pt\ta+,b-\t*\n");
        assert_eq!(dead_ends(&chain.graph), 2);
        let circle = reference(
            "S\ta\tACGT\nS\tb\tGG\nL\ta\t+\tb\t-\t0M\nL\tb\t-\ta\t+\t0M\nP\tone#1#pt\ta+,b-\t*\tTP:Z:circular\n",
        );
        assert_eq!(dead_ends(&circle.graph), 0);
    }

    #[test]
    fn a_closed_molecule_gets_its_join_segment_and_a_circular_walk() {
        let gfa = "H\tVN:Z:1.0\nS\ts.u0\tAAAACCCC\nS\ts.u1\tGGGTT\nL\ts.u0\t+\ts.u1\t-\t0M\nP\ts#1#pt.1\ts.u0+,s.u1-\t0M\tCH:i:1\n";
        let molecule = b"AAAACCCCAACCC";
        let closed = close_gfa(gfa, molecule, b"TGT", "read_pairs", "").unwrap().unwrap();
        assert!(closed.contains("S\ts.read_pairs_join\tTGT"));
        assert!(closed.contains("P\ts#1#pt.1\ts.u0+,s.u1-,s.read_pairs_join+\t0M,0M\tCH:i:1\tTP:Z:circular\tcj:Z:read_pairs"));
        let p = std::env::temp_dir().join(format!("ovasm-closed-{}.gfa", std::process::id()));
        fs::write(&p, &closed).unwrap();
        assert_eq!(dead_ends(&parse_gfa(&p).unwrap()), 0);
        fs::remove_file(p).unwrap();
        // no walk spells another sequence
        assert!(close_gfa(gfa, b"AAAACCCC", b"T", "read_pairs", "").unwrap().is_none());
    }

    #[test]
    fn reference_runs_rank_the_pairings_they_hold() {
        let (a, r, b, c) = (random(200), &random(400)[300..], &random(700)[400..], &random(1000)[700..]);
        let g = reference(&format!(
            "S\tA\t{a}\nS\tR\t{r}\nS\tB\t{b}\nS\tC\t{c}\nL\tA\t+\tR\t+\t0M\nL\tR\t+\tB\t+\t0M\nL\tR\t+\tC\t+\t0M\nP\tx#1#mt\tA+,R+,B+\t*\n"
        ));
        let index: BTreeMap<&str, usize> = g.graph.names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
        let s = |n: &str| (n.to_string(), true);
        let steps_b = vec![s("A"), s("R"), s("B")];
        let steps_c = vec![s("A"), s("R"), s("C")];
        let repeats: BTreeSet<String> = [String::from("R")].into();
        let pairings = vec![pairings_of(&steps_b, false, &repeats), pairings_of(&steps_c, false, &repeats)];
        assert_eq!(pairings[0].len(), 1);
        // two runs: one holds A R B, one A R C reversed; a third holds A R B again
        let walks = vec![
            ([a.as_bytes(), r.as_bytes(), b.as_bytes()].concat(), false),
            (revcomp(&[a.as_bytes(), r.as_bytes(), c.as_bytes()].concat()), false),
            ([a.as_bytes(), r.as_bytes(), b.as_bytes()].concat(), false),
        ];
        let names = vec!["one#1#mt".to_string(), "two#1#mt".to_string(), "three#1#mt".to_string()];
        let p = reference_prior(&g.graph, &index, &BTreeMap::new(), &pairings, &walks, &names);
        assert_eq!(p.runs, 3);
        assert_eq!(p.runs_with_all[0], vec!["one", "three"]);
        assert_eq!(p.runs_with_all[1], vec!["two"]);
        assert!(p.scores[0] > p.scores[1]);
    }

    #[test]
    fn a_molecule_has_one_key_however_it_is_written() {
        let s = |v: &[(&str, bool)]| v.iter().map(|(n, o)| (n.to_string(), *o)).collect::<Vec<_>>();
        let a = s(&[("u0", true), ("u9", true), ("u1", false)]);
        let rotated = s(&[("u9", true), ("u1", false), ("u0", true)]);
        let reversed = s(&[("u1", true), ("u9", false), ("u0", false)]);
        assert_eq!(walk_key(&a, true), walk_key(&rotated, true));
        assert_eq!(walk_key(&a, true), walk_key(&reversed, true));
        // a linear walk is the same backwards, not from another step
        assert_eq!(walk_key(&a, false), walk_key(&reversed, false));
        assert_ne!(walk_key(&a, false), walk_key(&rotated, false));
    }

    #[test]
    fn compaction_keeps_branches_and_all_bases() {
        let r=reference("S\ta\tAAAA\nS\tb\tC\nS\tc\tG\nS\td\tTTTT\nL\ta\t+\tb\t+\t0M\nL\ta\t+\tc\t+\t0M\nL\tb\t+\td\t+\t0M\nL\tc\t+\td\t+\t0M\nP\tone#1#pt\ta+,b+,d+\t*\nP\ttwo#1#pt\ta+,c+,d+\t*\n");
        let walks = unitigs(&r.graph);
        assert_eq!(walks.len(), 4);
        assert_eq!(walks.iter().map(|(w, _)| w.len()).sum::<usize>(), 4);
    }
    #[test]
    fn parallel_exact_alleles_merge_but_distinct_alleles_and_isolates_do_not() {
        let p = std::env::temp_dir().join(format!("ovasm-context-test-{}.gfa", std::process::id()));
        fs::write(&p,"H\tVN:Z:1.0\nS\tl\tAAAA\nS\tx\tC\nS\ty\tG\nS\ta\tA\nS\tr\tTTTT\nS\tisolated\tC\nL\tl\t+\tx\t+\t0M\nL\tx\t+\tr\t+\t0M\nL\tl\t+\ty\t-\t0M\nL\ty\t-\tr\t+\t0M\nL\tl\t+\ta\t+\t0M\nL\ta\t+\tr\t+\t0M\nP\tfirst#1#pt\tl+,x+,r+\t*\nP\tsecond#1#pt\tl+,y-,r+\t*\n").unwrap();
        assert_eq!(coalesce_parallel(&p).unwrap(), 1);
        let g = parse_gfa(&p).unwrap();
        assert_eq!(g.names.len(), 5);
        assert_eq!(g.links.len(), 4);
        assert!(g.names.iter().any(|n| n == "isolated"));
        let paths: Vec<_> = fs::read_to_string(&p)
            .unwrap()
            .lines()
            .filter(|s| s.starts_with("P\t"))
            .map(|s| s.split('\t').nth(2).unwrap().to_string())
            .collect();
        assert_eq!(paths[0], paths[1]);
        fs::remove_file(p.with_extension("parallel-nodes.tsv")).unwrap();
        fs::remove_file(p).unwrap();
    }
}
