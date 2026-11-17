//! Stage 3b: graph-to-graph pangenome (`ovasm pan merge`).
//!
//! Every sample's OV-GFA (its segments, links and configuration paths, isomers included) goes
//! into one graph, and bases the samples share become one base: a canonical k-mer (odd k, so
//! none is its own reverse complement) found in several places joins those places base by
//! base. The pieces are the non-branching runs of the joined base graph, cut at every sample
//! segment end: the base-identity blunting of `ovasm assemble unify`, across samples. Each sample
//! segment is then a walk over whole pieces, each sample link an adjacency of pieces and each
//! sample path (every configuration) a pangenome path. Linear pangenome builders take one
//! sequence per sample, so the isomers a sample's graph holds are lost before they start; here
//! none is. `run` rebuilds every sample segment, link and path from the pieces and fails unless
//! all of them match.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;

use crate::evidence::{parse_gfa, Graph, Oriented};
use crate::recruit::peak_rss_mb;
use crate::unify::{blunt_with, parse_walk, piece_walk, valid_pansn_field, StrandUf};

pub const FORMAT: &str = "organelleverse.pangraph.v1";

pub struct PanParams {
    /// Shared k-mer length that joins samples' bases (odd, at most 63).
    pub k: usize,
}

#[derive(Serialize)]
pub struct SampleSummary {
    pub sample: String,
    pub gfa: String,
    pub segments: usize,
    pub links: usize,
    pub paths: Vec<String>,
    pub bases: usize,
}

#[derive(Serialize)]
pub struct PanReport {
    pub format: &'static str,
    pub k: usize,
    pub samples: Vec<SampleSummary>,
    /// Bases of all sample segments together.
    pub sample_bases: usize,
    pub pieces: usize,
    pub pan_bases: usize,
    pub links: usize,
    /// Pieces (and their bases) in a segment of every sample.
    pub core_pieces: usize,
    pub core_bases: usize,
    /// Pieces by how many samples hold them: index = samples, value = (pieces, bases).
    pub pieces_by_samples: Vec<(usize, usize)>,
    /// Shared k-mer occurrences that joined bases, and joins refused because they would make a
    /// base its own reverse complement (left unjoined: less merging, nothing lost).
    pub kmer_joins: u64,
    pub join_conflicts: u64,
    pub verified_segments: usize,
    pub verified_links: usize,
    pub verified_paths: usize,
    pub elapsed_seconds: f64,
    pub peak_rss_mb: Option<f64>,
}

/// A sample path over pieces: (name, steps, circular, sample index).
type PanPath = (String, Vec<(usize, bool)>, bool, usize);

/// One input entry for `ovasm pan merge`: optional explicit sample name + GFA file or run directory.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PanInput {
    pub sample: Option<String>,
    pub path: PathBuf,
}

/// Parse a `.tsv` or `.txt` batch sample list (`sample<TAB>path`, `sample path`, or `path` per line).
pub fn parse_sample_list(list_path: &Path) -> Result<Vec<PanInput>> {
    let raw_text = fs::read_to_string(list_path)
        .with_context(|| format!("cannot read {}", list_path.display()))?;
    let text = raw_text.strip_prefix('\u{feff}').unwrap_or(&raw_text);
    let base_dir = list_path.parent().filter(|p| !p.as_os_str().is_empty());

    let resolve_rel = |raw: &str| -> PathBuf {
        let p = PathBuf::from(raw);
        if p.is_absolute() || p.exists() {
            p
        } else if let Some(b) = base_dir {
            let cand = b.join(raw);
            if cand.exists() { cand } else { p }
        } else {
            p
        }
    };
    let is_two_col_header = |c0: &str, c1: &str| -> bool {
        matches!(
            c0.to_ascii_lowercase().as_str(),
            "sample" | "sample_id" | "accession" | "id" | "name"
        ) && !resolve_rel(c1).exists()
    };
    let is_one_col_header = |c0: &str| -> bool {
        matches!(
            c0.to_ascii_lowercase().as_str(),
            "path" | "gfa" | "file" | "dir" | "directory"
        ) && !resolve_rel(c0).exists()
    };

    let mut out = Vec::new();
    let mut first_data_line = true;
    for raw_line in text.lines() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line.contains('\t') {
            let cols: Vec<&str> = line.split('\t').map(str::trim).collect();
            if cols.len() == 1 || cols.get(1).is_none_or(|c| c.is_empty()) {
                if first_data_line && is_one_col_header(cols[0]) {
                    first_data_line = false;
                    continue;
                }
                out.push(PanInput { sample: None, path: resolve_rel(cols[0]) });
            } else {
                if first_data_line && is_two_col_header(cols[0], cols[1]) {
                    first_data_line = false;
                    continue;
                }
                if cols[0].is_empty() {
                    bail!("{}: empty sample name in line {raw_line:?}", list_path.display());
                }
                out.push(PanInput {
                    sample: Some(cols[0].to_string()),
                    path: resolve_rel(cols[1]),
                });
            }
        } else if resolve_rel(line).exists() {
            out.push(PanInput { sample: None, path: resolve_rel(line) });
        } else {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() == 1 {
                if first_data_line && is_one_col_header(parts[0]) {
                    first_data_line = false;
                    continue;
                }
                out.push(PanInput { sample: None, path: resolve_rel(parts[0]) });
            } else {
                let rest = line[parts[0].len()..].trim();
                let raw_path = if resolve_rel(rest).exists() { rest } else { parts[1] };
                if first_data_line && is_two_col_header(parts[0], raw_path) {
                    first_data_line = false;
                    continue;
                }
                out.push(PanInput {
                    sample: Some(parts[0].to_string()),
                    path: resolve_rel(raw_path),
                });
            }
        }
        first_data_line = false;
    }
    if out.is_empty() {
        bail!("{}: no sample graphs listed", list_path.display());
    }
    Ok(out)
}

/// Resolve either a `.gfa` file or an `ovasm run` output directory to its `organelle.gfa`.
pub fn resolve_input_gfa(raw_path: &Path, organelle: Option<&str>) -> Result<PathBuf> {
    if raw_path.is_file() {
        return Ok(raw_path.to_path_buf());
    }
    if raw_path.is_dir() {
        if let Some(org) = organelle {
            let sub = raw_path.join(org).join("organelle.gfa");
            if sub.is_file() {
                return Ok(sub);
            }
            let direct = raw_path.join("organelle.gfa");
            if raw_path.file_name().and_then(|s| s.to_str()) == Some(org) && direct.is_file() {
                return Ok(direct);
            }
            bail!("{}: directory has no {org}/organelle.gfa", raw_path.display());
        }
        let direct = raw_path.join("organelle.gfa");
        if direct.is_file() {
            return Ok(direct);
        }
        if raw_path.join("mitochondrion/organelle.gfa").is_file()
            || raw_path.join("plastid/organelle.gfa").is_file()
        {
            bail!(
                "{}: run directory contains organelle subfolders; pass -G mt, -G pt or -G all (or point directly to organelle.gfa)",
                raw_path.display()
            );
        }
        bail!("{}: directory has no organelle.gfa", raw_path.display());
    }
    bail!("{}: no such GFA file or directory", raw_path.display())
}

/// One sample graph as read: its graph, P lines (name, step labels, circular) and name.
struct SampleGraph {
    sample: String,
    path: PathBuf,
    graph: Graph,
    paths: Vec<(String, String, bool)>,
    path_tags: BTreeMap<String, Vec<String>>,
}

fn rename_sample(sg: &mut SampleGraph, new_sample: &str) -> Result<()> {
    if !valid_pansn_field(new_sample) {
        bail!(
            "{}: sample {new_sample:?} must use only letters, digits, '_', '.' and '-'",
            sg.path.display()
        );
    }
    let old_sample = sg.sample.clone();
    let old_seg_prefix = format!("{old_sample}.");
    let new_seg_prefix = format!("{new_sample}.");
    if old_sample == new_sample
        && sg.graph.names.iter().all(|n| n.starts_with(&new_seg_prefix))
    {
        return Ok(());
    }
    let mut seg_map: FxHashMap<String, String> = FxHashMap::default();
    let mut seen_segs: FxHashSet<String> = FxHashSet::default();
    for name in &mut sg.graph.names {
        let old_name = name.clone();
        let new_name = if let Some(rest) = old_name.strip_prefix(&old_seg_prefix) {
            format!("{new_sample}.{rest}")
        } else if old_name.starts_with(&new_seg_prefix) {
            old_name.clone()
        } else {
            format!("{new_sample}.{old_name}")
        };
        if !seen_segs.insert(new_name.clone()) {
            bail!("{}: duplicate segment name {new_name} after renaming to {new_sample}", sg.path.display());
        }
        *name = new_name.clone();
        seg_map.insert(old_name, new_name);
    }

    let old_path_prefix = format!("{old_sample}#");
    let new_path_prefix = format!("{new_sample}#");
    let mut new_path_tags = BTreeMap::new();
    let mut seen_paths = FxHashSet::default();
    for (pname, labels, _) in &mut sg.paths {
        let old_pname = pname.clone();
        let new_pname = if let Some(rest) = old_pname.strip_prefix(&old_path_prefix) {
            format!("{new_sample}#{rest}")
        } else if old_pname.starts_with(&new_path_prefix) {
            old_pname.clone()
        } else {
            format!("{new_sample}#{old_pname}")
        };
        if !seen_paths.insert(new_pname.clone()) {
            bail!("{}: duplicate path name {new_pname} after renaming to {new_sample}", sg.path.display());
        }
        let new_labels = labels
            .split(' ')
            .map(|step| {
                let (seg, orient) = step.split_at(step.len() - 1);
                let mapped = seg_map.get(seg).map(String::as_str).unwrap_or(seg);
                format!("{mapped}{orient}")
            })
            .collect::<Vec<_>>()
            .join(" ");
        *pname = new_pname.clone();
        *labels = new_labels;
        if let Some(tags) = sg.path_tags.remove(&old_pname) {
            new_path_tags.insert(new_pname, tags);
        }
    }
    sg.path_tags = new_path_tags;
    sg.sample = new_sample.to_string();
    Ok(())
}

fn read_sample(path: &Path, sample_override: Option<&str>) -> Result<SampleGraph> {
    let text =
        fs::read_to_string(path).with_context(|| format!("cannot read {}", path.display()))?;
    let graph = parse_gfa(path)?;
    let filled_segments: FxHashSet<&str> = text.lines().filter_map(|line| {
        let f: Vec<_> = line.split('\t').collect();
        (f.first() == Some(&"S") && f.len() >= 3 && f[3..].iter().any(|t| t.starts_with("rf:"))).then(|| f[1])
    }).collect();
    // Segment tags are lost during graph resegmentation, so transfer uncertain
    // segment provenance to every input path traversing that segment.
    let uncertain_segments: FxHashSet<&str> = text.lines().filter_map(|line| {
        let f: Vec<_> = line.split('\t').collect();
        (f.first() == Some(&"S") && f.len() >= 3 && f[3..].iter().any(|t| t.starts_with("up:") || *t == "st:Z:unresolved_configuration")).then(|| f[1])
    }).collect();
    let mut sample = None;
    let mut paths = Vec::new();
    let mut path_tags = BTreeMap::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        match f.first() {
            Some(&"H") => {
                sample = f
                    .iter()
                    .find_map(|t| t.strip_prefix("sn:Z:"))
                    .map(str::to_string);
            }
            Some(&"P") if f.len() >= 3 => {
                let steps: Vec<&str> = f[2].split(',').collect();
                let walk = steps
                    .iter()
                    .map(|s| {
                        let (name, o) = s.split_at(s.len() - 1);
                        if o != "+" && o != "-" {
                            bail!("{}: path step {s:?} has no orientation", path.display());
                        }
                        Ok(format!("{name}{o}"))
                    })
                    .collect::<Result<Vec<_>>>()?
                    .join(" ");
                let circular = f[3..].contains(&"TP:Z:circular");
                let mut tags: Vec<String> = f.get(4..).unwrap_or(&[]).iter()
                    .filter(|t| !t.starts_with("TP:" )).map(|t| t.to_string()).collect();
                if steps.iter().any(|s| filled_segments.contains(&s[..s.len()-1])) && !tags.iter().any(|t| t.starts_with("rf:")) {
                    tags.push("rf:Z:contains_reference_fill".into());
                }
                if steps.iter().any(|s| uncertain_segments.contains(&s[..s.len()-1])) && !tags.iter().any(|t| t.starts_with("up:")) {
                    tags.push("up:Z:unresolved_segment".into());
                }
                path_tags.insert(f[1].to_string(), tags);
                paths.push((f[1].to_string(), walk, circular));
            }
            _ => {}
        }
    }
    let sample = match sample {
        Some(s) => s,
        None => path
            .file_name()
            .and_then(|n| n.to_str())
            .and_then(|n| n.split('.').next())
            .context("no sample name (H sn:Z:) and no usable file name")?
            .to_string(),
    };
    if sample_override.is_none() && !valid_pansn_field(&sample) {
        bail!(
            "{}: sample {sample:?} must use only letters, digits, '_', '.' and '-'",
            path.display()
        );
    }
    if graph.overlaps.iter().any(|&o| o != 0) {
        bail!(
            "{}: links must be blunt (0M); run `ovasm assemble unify` first",
            path.display()
        );
    }
    // Where each unresolved walk's own graph branches (`uj:Z:`, offsets in bases from the walk's
    // start, 0 for a circular walk's closing junction, `-` for none): the sample's reads did not
    // decide what follows what across a link that is not the only one at both of its ends, so
    // the walk is evidence of the sample's junctions only between these offsets
    // (`pan_posterior`: an Arabidopsis HiFi assembly unresolved at its two recombining repeats
    // still shows how every shorter repeat is joined). Only for the graph of one sample: a
    // merged graph branches wherever any of its samples does, and its walks keep the offsets
    // they came with.
    if paths.iter().all(|(name, _, _)| name.split('#').next() == Some(sample.as_str())) {
        let index: FxHashMap<&str, usize> = graph.names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
        // links leaving each oriented segment
        let mut leaving: FxHashMap<Oriented, usize> = FxHashMap::default();
        for &(a, b) in &graph.links {
            *leaving.entry(a).or_default() += 1;
            *leaving.entry((b.0, !b.1)).or_default() += 1;
        }
        let branches = |x: Oriented, y: Oriented| {
            leaving.get(&x).copied().unwrap_or(0) != 1 || leaving.get(&(y.0, !y.1)).copied().unwrap_or(0) != 1
        };
        for (name, labels, circular) in &paths {
            let Some(tags) = path_tags.get_mut(name) else { continue };
            let unresolved = tags.iter().any(|t| t == "st:Z:unresolved_configuration" || t.starts_with("up:"));
            if !unresolved || tags.iter().any(|t| t.starts_with("uj:")) {
                continue;
            }
            let walk: Vec<Oriented> = labels
                .split(' ')
                .filter_map(|l| {
                    let (n, o) = l.split_at(l.len() - 1);
                    index.get(n).map(|&i| (i, o == "+"))
                })
                .collect();
            if walk.is_empty() || walk.len() != labels.split(' ').count() {
                continue; // an unknown segment: `parse_walk` reports it
            }
            let (mut cuts, mut offset) = (Vec::new(), 0usize);
            if *circular && branches(walk[walk.len() - 1], walk[0]) {
                cuts.push(0);
            }
            for w in walk.windows(2) {
                offset += graph.seqs[w[0].0].len();
                if branches(w[0], w[1]) {
                    cuts.push(offset);
                }
            }
            tags.push(if cuts.is_empty() {
                "uj:Z:-".to_string()
            } else {
                format!("uj:Z:{}", cuts.iter().map(usize::to_string).collect::<Vec<_>>().join(","))
            });
        }
    }
    let mut sg = SampleGraph {
        sample,
        path: path.to_path_buf(),
        graph,
        paths,
        path_tags,
    };
    if let Some(new_sample) = sample_override {
        rename_sample(&mut sg, new_sample)?;
    }
    Ok(sg)
}

#[inline]
fn code(b: u8) -> Option<u128> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// `(start, canonical k-mer, forward?)` of every ACGT k-mer of `seq`, 2-bit packed (k <= 63).
fn kmers(seq: &[u8], k: usize) -> Vec<(usize, u128, bool)> {
    let mask: u128 = (1u128 << (2 * k)) - 1;
    let shift = 2 * (k as u32 - 1);
    let (mut fwd, mut rev, mut valid) = (0u128, 0u128, 0usize);
    let mut out = Vec::with_capacity(seq.len().saturating_sub(k - 1));
    for (i, &b) in seq.iter().enumerate() {
        match code(b) {
            Some(c) => {
                fwd = ((fwd << 2) | c) & mask;
                rev = (rev >> 2) | ((3 - c) << shift);
                valid += 1;
                if valid >= k {
                    // odd k: a k-mer never equals its reverse complement
                    out.push((i + 1 - k, fwd.min(rev), fwd < rev));
                }
            }
            None => {
                valid = 0;
                fwd = 0;
                rev = 0;
            }
        }
    }
    out
}

/// Join every occurrence of a shared canonical k-mer to its first occurrence, base by base.
/// Along a run of k-mers that match consecutively (same diagonal) only the new last base is
/// joined. Returns (k-mer occurrences joined, joins refused as self-complementary).
fn join_shared(g: &Graph, k: usize, uf: &mut StrandUf, base0: &[usize]) -> (u64, u64) {
    // canonical k-mer -> (segment, start, forward?) of its first occurrence
    let mut first: FxHashMap<u128, (u32, u32, bool)> = FxHashMap::default();
    let (mut joined, mut conflicts) = (0u64, 0u64);
    for (s, seq) in g.seqs.iter().enumerate() {
        // previous k-mer's partner: (segment, start, parity)
        let mut prev: Option<(usize, u32, usize, u8)> = None;
        for (i, km, fw) in kmers(seq, k) {
            let &mut (s0, i0, f0) = first.entry(km).or_insert((s as u32, i as u32, fw));
            if s0 as usize == s && i0 as usize == i {
                prev = None;
                continue;
            }
            joined += 1;
            let p = u8::from(fw != f0);
            // base i + t pairs with i0 + t (same strand) or i0 + k - 1 - t (opposite)
            let partner = |t: usize| -> usize {
                if p == 0 {
                    i0 as usize + t
                } else {
                    i0 as usize + k - 1 - t
                }
            };
            let continues = matches!(prev, Some((pi, ps0, pi0, pp))
                if pi + 1 == i && ps0 == s0 && pp == p
                    && if p == 0 { pi0 + 1 == i0 as usize } else { i0 as usize + 1 == pi0 });
            let ts: Vec<usize> = if continues {
                vec![k - 1]
            } else {
                (0..k).collect()
            };
            for t in ts {
                if !uf.union(base0[s] + i + t, base0[s0 as usize] + partner(t), p) {
                    conflicts += 1;
                }
            }
            prev = Some((i, s0, i0 as usize, p));
        }
    }
    (joined, conflicts)
}

pub fn run(inputs: &[PathBuf], p: &PanParams, out: &Path, out_json: &Path) -> Result<PanReport> {
    let named: Vec<PanInput> = inputs
        .iter()
        .map(|path| PanInput { sample: None, path: path.clone() })
        .collect();
    run_named(&named, None, p, out, out_json)
}

pub fn run_named(
    inputs: &[PanInput],
    organelle: Option<&str>,
    p: &PanParams,
    out: &Path,
    out_json: &Path,
) -> Result<PanReport> {
    let t0 = Instant::now();
    if p.k % 2 == 0 || !(15..=63).contains(&p.k) {
        bail!("k must be odd and in 15..=63 (got {})", p.k);
    }
    let samples = inputs
        .iter()
        .map(|entry| {
            let gfa = resolve_input_gfa(&entry.path, organelle)?;
            read_sample(&gfa, entry.sample.as_deref())
        })
        .collect::<Result<Vec<_>>>()?;
    let mut seen_samples = FxHashSet::default();
    for sg in &samples {
        if !seen_samples.insert(sg.sample.clone()) {
            bail!("sample {} is given twice", sg.sample);
        }
    }

    // One graph of every sample's segments and links.
    let mut g = Graph {
        names: Vec::new(),
        seqs: Vec::new(),
        links: Vec::new(),
        overlaps: Vec::new(),
    };
    let mut owner: Vec<usize> = Vec::new(); // segment -> sample
    let mut offset: Vec<usize> = Vec::new(); // sample -> first segment
    let mut names: FxHashSet<String> = FxHashSet::default();
    for (si, sg) in samples.iter().enumerate() {
        let off = g.names.len();
        offset.push(off);
        for (name, seq) in sg.graph.names.iter().zip(&sg.graph.seqs) {
            if !names.insert(name.clone()) {
                bail!(
                    "segment {name} occurs in two samples' graphs; OV-GFA names carry the sample"
                );
            }
            g.names.push(name.clone());
            g.seqs.push(seq.clone());
            owner.push(si);
        }
        for (&(a, b), &ov) in sg.graph.links.iter().zip(&sg.graph.overlaps) {
            g.links.push(((a.0 + off, a.1), (b.0 + off, b.1)));
            g.overlaps.push(ov);
        }
    }

    let (mut kmer_joins, mut join_conflicts) = (0, 0);
    let b = blunt_with(&g, |uf, base0| {
        (kmer_joins, join_conflicts) = join_shared(&g, p.k, uf, base0);
        Ok(())
    })?;

    // Which samples hold each piece.
    let n_samples = samples.len();
    let mut holders: Vec<Vec<usize>> = vec![Vec::new(); b.pieces.len()];
    for (s, walk) in b.walks.iter().enumerate() {
        for &(piece, _, _) in walk {
            let h = &mut holders[piece];
            if h.last() != Some(&owner[s]) && !h.contains(&owner[s]) {
                h.push(owner[s]);
            }
        }
    }
    let piece_len = |i: usize| b.pieces[i].end - b.pieces[i].start;
    let sign = |f: bool| if f { '+' } else { '-' };

    // Links: consecutive pieces inside segments, and sample links; with the samples holding each.
    type PieceEnd = (usize, bool);
    let canon = |x: PieceEnd, y: PieceEnd| {
        let rc = ((y.0, !y.1), (x.0, !x.1));
        if (x, y) <= rc {
            (x, y)
        } else {
            rc
        }
    };
    let mut link_samples: BTreeMap<(PieceEnd, PieceEnd), Vec<usize>> = BTreeMap::new();
    let mut add_link = |key: (PieceEnd, PieceEnd), sample: usize| {
        let v = link_samples.entry(key).or_default();
        if !v.contains(&sample) {
            v.push(sample);
        }
    };
    for (s, w) in b.walks.iter().enumerate() {
        for pair in w.windows(2) {
            add_link(
                canon((pair[0].0, pair[0].1), (pair[1].0, pair[1].1)),
                owner[s],
            );
        }
    }
    let mut sample_link_keys = Vec::with_capacity(b.links.len());
    for &((a, c), _) in &b.links {
        let from = *b.after(&g, a, 0)?.last().context("empty segment")?;
        let to = *b.after(&g, c, 0)?.first().context("empty segment")?;
        let key = canon(from, to);
        add_link(key, owner[a.0]);
        sample_link_keys.push(key);
    }

    // Pangenome paths: every sample path, over pieces.
    let mut pan_paths: Vec<PanPath> = Vec::new();
    for (si, sg) in samples.iter().enumerate() {
        for (name, labels, circular) in &sg.paths {
            let local = parse_walk(&sg.graph, labels)
                .with_context(|| format!("{}: path {name}", sg.path.display()))?;
            let walk: Vec<Oriented> = local.iter().map(|&(s, f)| (s + offset[si], f)).collect();
            let steps = piece_walk(&g, &b, &walk, *circular)?;
            pan_paths.push((name.clone(), steps, *circular, si));
        }
    }

    // Verify: every segment, link and path comes back exactly.
    let spell = |steps: &[(usize, bool)]| -> Vec<u8> {
        let mut out = Vec::new();
        for &(id, f) in steps {
            let pc = b.pieces[id];
            let seq = &g.seqs[pc.owner][pc.start..pc.end];
            if f {
                out.extend_from_slice(seq);
            } else {
                out.extend(crate::evidence::revcomp(seq));
            }
        }
        out
    };
    for (s, w) in b.walks.iter().enumerate() {
        let steps: Vec<(usize, bool)> = w.iter().map(|&(id, f, _)| (id, f)).collect();
        if spell(&steps) != g.seqs[s] {
            bail!(
                "internal: segment {} does not come back from its pieces",
                g.names[s]
            );
        }
    }
    for (i, key) in sample_link_keys.iter().enumerate() {
        if !link_samples.contains_key(key) {
            bail!("internal: sample link {i} has no pangenome link");
        }
    }
    for (name, steps, _, si) in &pan_paths {
        let (_, labels, _) = samples[*si]
            .paths
            .iter()
            .find(|(n, _, _)| n == name)
            .context("path vanished")?;
        let local = parse_walk(&samples[*si].graph, labels)?;
        let mut want = Vec::new();
        for &o in &local {
            want.extend(samples[*si].graph.oriented_seq(o));
        }
        if spell(steps) != want {
            bail!("internal: path {name} does not come back from the pangenome");
        }
        for pair in steps.windows(2) {
            if !link_samples.contains_key(&canon(pair[0], pair[1])) {
                bail!("internal: path {name} steps over a missing link");
            }
        }
    }

    // GFA: pieces named 1..N (integers, as ODGI and vg expect), ns:i = samples holding them.
    let mut gfa = String::new();
    writeln!(
        gfa,
        "H\tVN:Z:1.0\tov:Z:{FORMAT}\tsn:i:{n_samples}\tkm:i:{}",
        p.k
    )?;
    for (id, pc) in b.pieces.iter().enumerate() {
        writeln!(
            gfa,
            "S\t{}\t{}\tLN:i:{}\tns:i:{}",
            id + 1,
            std::str::from_utf8(&g.seqs[pc.owner][pc.start..pc.end])?,
            pc.end - pc.start,
            holders[id].len()
        )?;
    }
    for (&(x, y), held) in &link_samples {
        writeln!(
            gfa,
            "L\t{}\t{}\t{}\t{}\t0M\tns:i:{}",
            x.0 + 1,
            sign(x.1),
            y.0 + 1,
            sign(y.1),
            held.len()
        )?;
    }
    for (name, steps, circular, si) in &pan_paths {
        let walk: Vec<String> = steps
            .iter()
            .map(|&(id, f)| format!("{}{}", id + 1, sign(f)))
            .collect();
        let cigar = vec!["0M"; steps.len().saturating_sub(1)].join(",");
        write!(
            gfa,
            "P\t{name}\t{}\t{}",
            walk.join(","),
            if cigar.is_empty() { "*" } else { &cigar }
        )?;
        if *circular {
            gfa.push_str("\tTP:Z:circular");
        }
        if let Some(tags) = samples[*si].path_tags.get(name) {
            for tag in tags { write!(gfa, "\t{tag}")?; }
        }
        gfa.push('\n');
    }
    fs::create_dir_all(out.parent().unwrap_or(Path::new(".")))?;
    fs::write(out, gfa)?;

    // Sample segments as piece walks and sample links as piece links (rebuilds every input).
    let mut tsv = String::from("sample\tsegment\tlength\tpieces\n");
    for (s, w) in b.walks.iter().enumerate() {
        let steps: Vec<String> = w
            .iter()
            .map(|&(id, f, _)| format!("{}{}", id + 1, sign(f)))
            .collect();
        writeln!(
            tsv,
            "{}\t{}\t{}\t{}",
            samples[owner[s]].sample,
            g.names[s],
            g.seqs[s].len(),
            steps.join(",")
        )?;
    }
    let seg_path = out.with_extension("segments.tsv");
    fs::write(&seg_path, tsv)?;
    let mut ltsv = String::from("sample\tfrom\tto\tpangenome_link\n");
    for (i, &((a, c), _)) in b.links.iter().enumerate() {
        let (x, y) = sample_link_keys[i];
        writeln!(
            ltsv,
            "{}\t{}\t{}\t{}{} {}{}",
            samples[owner[a.0]].sample,
            g.label(a),
            g.label(c),
            x.0 + 1,
            sign(x.1),
            y.0 + 1,
            sign(y.1)
        )?;
    }
    fs::write(out.with_extension("links.tsv"), ltsv)?;

    let mut by = vec![(0usize, 0usize); n_samples + 1];
    for (id, h) in holders.iter().enumerate() {
        by[h.len()].0 += 1;
        by[h.len()].1 += piece_len(id);
    }
    let core: Vec<usize> = (0..b.pieces.len())
        .filter(|&i| holders[i].len() == n_samples)
        .collect();
    let report = PanReport {
        format: FORMAT,
        k: p.k,
        samples: samples
            .iter()
            .enumerate()
            .map(|(si, sg)| SampleSummary {
                sample: sg.sample.clone(),
                gfa: sg.path.display().to_string(),
                segments: sg.graph.names.len(),
                links: sg.graph.links.len(),
                paths: sg.paths.iter().map(|x| x.0.clone()).collect(),
                bases: (0..sg.graph.seqs.len())
                    .map(|s| g.seqs[s + offset[si]].len())
                    .sum(),
            })
            .collect(),
        sample_bases: g.seqs.iter().map(Vec::len).sum(),
        pieces: b.pieces.len(),
        pan_bases: (0..b.pieces.len()).map(piece_len).sum(),
        links: link_samples.len(),
        core_pieces: core.len(),
        core_bases: core.iter().map(|&i| piece_len(i)).sum(),
        pieces_by_samples: by,
        kmer_joins,
        join_conflicts,
        verified_segments: b.walks.len(),
        verified_links: sample_link_keys.len(),
        verified_paths: pan_paths.len(),
        elapsed_seconds: t0.elapsed().as_secs_f64(),
        peak_rss_mb: peak_rss_mb(),
    };
    fs::write(out_json, serde_json::to_string_pretty(&report)? + "\n")?;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(n: usize, seed: u64) -> String {
        let mut s = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        (0..n)
            .map(|_| {
                s ^= s << 13;
                s ^= s >> 7;
                s ^= s << 17;
                ['A', 'C', 'G', 'T'][(s % 4) as usize]
            })
            .collect()
    }

    fn rc(s: &str) -> String {
        String::from_utf8(crate::evidence::revcomp(s.as_bytes())).unwrap()
    }

    /// A repeat R between four unique arms; `flip` stores sample's u1 reverse complemented.
    fn sample_gfa(dir: &Path, name: &str, a: &str, b: &str, r: &str, snp: bool) -> PathBuf {
        let mut b2 = b.to_string();
        if snp {
            let mid = b2.len() / 2;
            let c = if &b2[mid..=mid] == "A" { "C" } else { "A" };
            b2.replace_range(mid..=mid, c);
        }
        let text = format!(
            "H\tVN:Z:1.0\tsn:Z:{name}\n\
             S\t{name}.u0\t{a}\n\
             S\t{name}.u1\t{}\n\
             S\t{name}.u2\t{r}\n\
             L\t{name}.u0\t+\t{name}.u2\t+\t0M\n\
             L\t{name}.u2\t+\t{name}.u1\t-\t0M\n\
             L\t{name}.u1\t-\t{name}.u2\t+\t0M\n\
             L\t{name}.u2\t+\t{name}.u0\t+\t0M\n\
             P\t{name}#1#mt\t{name}.u0+,{name}.u2+,{name}.u1-,{name}.u2+\t0M,0M,0M\tTP:Z:circular\n\
             P\t{name}#2#mt\t{name}.u2+,{name}.u1-,{name}.u2+,{name}.u0+\t0M,0M,0M\n",
            rc(&b2)
        );
        let p = dir.join(format!("{name}.ov.gfa"));
        fs::write(&p, text).unwrap();
        p
    }

    #[test]
    fn samples_merge_and_every_segment_link_and_path_comes_back() {
        let dir = std::env::temp_dir().join(format!("ovasm-pan-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let (a, b, r) = (seq(3000, 1), seq(2500, 2), seq(800, 3));
        let s1 = sample_gfa(&dir, "S1", &a, &b, &r, false);
        let s2 = sample_gfa(&dir, "S2", &a, &b, &r, true);
        let out = dir.join("pan.gfa");
        let rep = run(&[s1, s2], &PanParams { k: 31 }, &out, &dir.join("pan.json")).unwrap();
        assert_eq!(rep.verified_paths, 4);
        assert_eq!(rep.verified_segments, 6);
        // shared sequence stored once: a SNP in b adds only the bases around it
        assert!(
            rep.pan_bases < rep.sample_bases / 2 + 100,
            "{} of {}",
            rep.pan_bases,
            rep.sample_bases
        );
        assert!(rep.core_bases >= a.len() + r.len() + b.len() - 2 * 31);
        let text = fs::read_to_string(&out).unwrap();
        assert_eq!(text.lines().filter(|l| l.starts_with("P\t")).count(), 4);
        assert!(text
            .lines()
            .any(|l| l.starts_with("P\tS1#1#mt\t") && l.ends_with("TP:Z:circular")));
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn unresolved_path_provenance_survives_pan_merge() {
        let dir = std::env::temp_dir().join(format!("ovasm-pan-tags-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let input = sample_gfa(&dir, "tags", &seq(300, 1), &seq(300, 2), &seq(100, 3), false);
        let text = fs::read_to_string(&input).unwrap().replace("TP:Z:circular", "TP:Z:circular\tst:Z:unresolved_configuration\trf:Z:source");
        fs::write(&input, text).unwrap();
        let out = dir.join("pan.gfa");
        run(&[input], &PanParams { k: 31 }, &out, &dir.join("pan.json")).unwrap();
        let output = fs::read_to_string(out).unwrap();
        assert!(output.lines().any(|l| l.starts_with("P\ttags#1#mt\t") && l.contains("st:Z:unresolved_configuration") && l.contains("rf:Z:source")));
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn an_unresolved_walk_carries_the_offsets_where_its_own_graph_branches() {
        let dir = std::env::temp_dir().join(format!("ovasm-pan-uj-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let input = sample_gfa(&dir, "S", &seq(300, 1), &seq(300, 2), &seq(100, 3), false);
        let text = fs::read_to_string(&input).unwrap().replace("TP:Z:circular", "TP:Z:circular\tst:Z:unresolved_configuration");
        fs::write(&input, text).unwrap();
        let out = dir.join("pan.gfa");
        run(&[input], &PanParams { k: 31 }, &out, &dir.join("pan.json")).unwrap();
        let tags = |text: &str, name: &str| -> Vec<String> {
            let line = text.lines().find(|l| l.starts_with(&format!("P\t{name}\t"))).unwrap();
            line.split('\t').skip(4).map(str::to_string).collect()
        };
        let merged = fs::read_to_string(&out).unwrap();
        // u0 (300) -> u2 (100) -> u1 (300) -> u2 -> back to u0: the repeat u2 branches at both ends
        assert!(tags(&merged, "S#1#mt").contains(&"uj:Z:0,300,400,700".to_string()), "{:?}", tags(&merged, "S#1#mt"));
        // a walk the reads decided carries none
        assert!(!tags(&merged, "S#2#mt").iter().any(|t| t.starts_with("uj:")));
        // a merged graph is no sample's own: its walks keep the offsets they came with
        let other = sample_gfa(&dir, "T", &seq(300, 1), &seq(300, 2), &seq(100, 3), true);
        let again = dir.join("again.gfa");
        run(&[out, other], &PanParams { k: 31 }, &again, &dir.join("again.json")).unwrap();
        let merged = fs::read_to_string(&again).unwrap();
        assert!(tags(&merged, "S#1#mt").contains(&"uj:Z:0,300,400,700".to_string()));
        assert!(!tags(&merged, "T#1#mt").iter().any(|t| t.starts_with("uj:")));
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn kmers_are_canonical_and_strand_tagged() {
        let s = seq(200, 9);
        let f = kmers(s.as_bytes(), 31);
        let r = kmers(rc(&s).as_bytes(), 31);
        assert_eq!(f.len(), 170);
        // the k-mer at i forward is the one at len - k - i on the reverse strand
        for &(i, km, fw) in &f {
            let (_, km2, fw2) = r[200 - 31 - i];
            assert_eq!(km, km2);
            assert_ne!(fw, fw2);
        }
    }

    #[test]
    fn rejects_even_k_and_overlapping_links() {
        let dir = std::env::temp_dir().join(format!("ovasm-pan-bad-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let (a, b, r) = (seq(500, 4), seq(500, 5), seq(300, 6));
        let s1 = sample_gfa(&dir, "S1", &a, &b, &r, false);
        assert!(run(
            std::slice::from_ref(&s1),
            &PanParams { k: 32 },
            &dir.join("x.gfa"),
            &dir.join("x.json")
        )
        .is_err());
        let text = fs::read_to_string(&s1)
            .unwrap()
            .replacen("\t0M\n", "\t5M\n", 1);
        fs::write(&s1, text).unwrap();
        assert!(run(
            &[s1],
            &PanParams { k: 31 },
            &dir.join("y.gfa"),
            &dir.join("y.json")
        )
        .is_err());
        fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn sample_table_tsv_and_txt_rename_default_sample_graphs_and_resolve_run_dirs() {
        let dir = std::env::temp_dir().join(format!("ovasm-pan-tsv-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        let d1 = dir.join("acc1").join("mitochondrion");
        let d2 = dir.join("acc2").join("mitochondrion");
        fs::create_dir_all(&d1).unwrap();
        fs::create_dir_all(&d2).unwrap();
        let (a, b, r) = (seq(1200, 11), seq(1100, 12), seq(300, 13));
        // Both runs used the default `--sample sample`, which would collide without column 1.
        let g1_tmp = sample_gfa(&d1, "sample", &a, &b, &r, false);
        let g2_tmp = sample_gfa(&d2, "sample", &a, &b, &r, true);
        let g1 = d1.join("organelle.gfa");
        let g2 = d2.join("organelle.gfa");
        let t1 = fs::read_to_string(&g1_tmp)
            .unwrap()
            .replace("TP:Z:circular", "TP:Z:circular\tst:Z:unresolved_configuration");
        fs::write(&g1, t1).unwrap();
        fs::rename(&g2_tmp, &g2).unwrap();
        fs::remove_file(&g1_tmp).unwrap();

        // 1) Two-column TSV with header and comment, pointing directly to GFA files (relative to TSV).
        let tsv = dir.join("samples.tsv");
        fs::write(
            &tsv,
            "# batch merge table\nsample\tpath\nACC_001\tacc1/mitochondrion/organelle.gfa\nACC_002\tacc2/mitochondrion/organelle.gfa\n",
        )
        .unwrap();
        let parsed = parse_sample_list(&tsv).unwrap();
        assert_eq!(parsed.len(), 2);
        assert_eq!(parsed[0].sample.as_deref(), Some("ACC_001"));
        assert_eq!(parsed[1].sample.as_deref(), Some("ACC_002"));

        let out_gfa = dir.join("pan.mt.gfa");
        let out_json = dir.join("pan.mt.json");
        let rep = run_named(&parsed, None, &PanParams { k: 31 }, &out_gfa, &out_json).unwrap();
        assert_eq!(rep.samples.len(), 2);
        assert_eq!(rep.samples[0].sample, "ACC_001");
        assert_eq!(rep.samples[1].sample, "ACC_002");
        assert_eq!(rep.verified_paths, 4);
        let merged = fs::read_to_string(&out_gfa).unwrap();
        assert!(merged.lines().any(|l| l.starts_with("P\tACC_001#1#mt\t") && l.contains("uj:Z:")));
        assert!(merged.lines().any(|l| l.starts_with("P\tACC_002#1#mt\t")));
        let segs = fs::read_to_string(out_gfa.with_extension("segments.tsv")).unwrap();
        assert!(segs.contains("ACC_001\tACC_001.u0\t"));
        assert!(segs.contains("ACC_002\tACC_002.u0\t"));

        // 2) Space-separated TXT pointing to parent run directories (`acc1`, `acc2`) + organelle = mitochondrion.
        let txt = dir.join("samples.txt");
        fs::write(&txt, "sample dir\nACC_001 acc1\nACC_002 acc2\n").unwrap();
        let parsed_txt = parse_sample_list(&txt).unwrap();
        let out_gfa2 = dir.join("pan2.mt.gfa");
        let rep2 = run_named(
            &parsed_txt,
            Some("mitochondrion"),
            &PanParams { k: 31 },
            &out_gfa2,
            &dir.join("pan2.mt.json"),
        )
        .unwrap();
        assert_eq!(rep2.pan_bases, rep.pan_bases);
        assert_eq!(rep2.verified_paths, 4);
        fs::remove_dir_all(&dir).unwrap();
    }
}
