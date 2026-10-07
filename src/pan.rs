//! Stage 3b: graph-to-graph pangenome (`ovasm pan`).
//!
//! Every sample's OV-GFA (its segments, links and configuration paths, isomers included) goes
//! into one graph, and bases the samples share become one base: a canonical k-mer (odd k, so
//! none is its own reverse complement) found in several places joins those places base by
//! base. The pieces are the non-branching runs of the joined base graph, cut at every sample
//! segment end: the base-identity blunting of `ovasm unify`, across samples. Each sample
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

/// One sample graph as read: its graph, P lines (name, step labels, circular) and name.
struct SampleGraph {
    sample: String,
    path: PathBuf,
    graph: Graph,
    paths: Vec<(String, String, bool)>,
}

fn read_sample(path: &Path) -> Result<SampleGraph> {
    let text =
        fs::read_to_string(path).with_context(|| format!("cannot read {}", path.display()))?;
    let graph = parse_gfa(path)?;
    let mut sample = None;
    let mut paths = Vec::new();
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
    if !valid_pansn_field(&sample) {
        bail!(
            "{}: sample {sample:?} must use only letters, digits, '_', '.' and '-'",
            path.display()
        );
    }
    if graph.overlaps.iter().any(|&o| o != 0) {
        bail!(
            "{}: links must be blunt (0M); run `ovasm unify` first",
            path.display()
        );
    }
    Ok(SampleGraph {
        sample,
        path: path.to_path_buf(),
        graph,
        paths,
    })
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
    let t0 = Instant::now();
    if p.k % 2 == 0 || !(15..=63).contains(&p.k) {
        bail!("k must be odd and in 15..=63 (got {})", p.k);
    }
    let samples = inputs
        .iter()
        .map(|f| read_sample(f))
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
    for (name, steps, circular, _) in &pan_paths {
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
}
