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
use crate::{pan, pipeline, recruit};

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

    fn add_bait(&self, bait: &mut BaitSet, name: &str, out: &Path) -> Result<Value> {
        let t = bait.target(name)?;
        let mut fasta = String::new();
        for (i, w) in self.walks.iter().enumerate() {
            let mut seq = w.sequence.clone();
            if w.circular {
                seq.extend(w.sequence.iter().cycle().take(K - 1));
            }
            bait.add_sequence(t, &seq);
            writeln!(
                fasta,
                ">{i} source_path={}\n{}",
                w.name,
                std::str::from_utf8(&seq)?
            )?;
        }
        // Also retain graph sequence not represented by the finite source paths.
        for (i, seq) in self.graph.seqs.iter().enumerate() {
            bait.add_sequence(t, seq);
            writeln!(
                fasta,
                ">node_{}\n{}",
                self.graph.names[i],
                std::str::from_utf8(seq)?
            )?;
        }
        fs::write(out, fasta)?;
        let basis: u64 = self.graph.seqs.iter().map(|s| s.len() as u64).sum();
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

pub struct Config {
    pub reads: Vec<PathBuf>,
    pub target_reads: BTreeMap<String, Vec<PathBuf>>,
    pub references: BTreeMap<String, PathBuf>,
    pub sample: String,
    pub out: PathBuf,
    pub target_depth: Option<f64>,
    pub salt: u64,
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
        let assembly = pipeline::run(&pipeline::RunConfig {
            reads,
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
        organelles.insert(
            o,
            json!({"reference":reference_info[o],"projection":projection,"individual":integration,
            "de_novo":match assembly {Ok(s)=>s,Err(e)=>json!({"error":format!("{e:#}")})} }),
        );
    }
    let workflow_completed = organelles
        .values()
        .all(|s| s["de_novo"].get("error").is_none());
    let result = json!({"tool":"ovasm reconstruct","format":"organelleverse.reference-reconstruction.v1","sample":cfg.sample,
        "workflow_completed":workflow_completed,
        "inputs":cfg.reads,"target_read_inputs":cfg.target_reads,
        "input_scope":if targeted {"explicit_target_reads"} else {"whole_genome_recruitment"},
        "reference_graphs":cfg.references,"recruitment":recruitment,"organelles":organelles,
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
