//! Stage 3: one graph format for every assembly backend (`ovasm unify`).
//!
//! OV-GFA v1 is GFA1 with blunt (0M) links, segments renamed `<sample>.<name>[.<piece>]`, one
//! depth tag `dp:f` whatever the backend wrote, read support `ev:i` on links (from `ovasm
//! evidence`), and the representative linearization plus near-tied alternatives as PanSN P lines
//! `<sample>#<hap>#<molecule>`. The pangenome module's parser and its overview and conversion
//! stages require exactly this: zero overlaps, P paths and stored sequences.
//!
//! Blunting identifies bases: the bases a link's overlap declares the same become one base,
//! and the pieces are the non-branching runs of the resulting base graph (see `blunt`). Every
//! link and every path of the input survives; no base is written twice.

use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use anyhow::{bail, Context, Result};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use serde_json::Value;

use crate::evidence::{flip, parse_gfa, Graph, Oriented};

pub const FORMAT: &str = "organelleverse.graph.v1";

pub struct UnifyParams<'a> {
    pub sample: &'a str,
    pub backend: &'a str,
    pub organelle: &'a str,
    /// Molecule field of the PanSN path names, e.g. `mt` or `pt`.
    pub molecule: &'a str,
    /// Alternatives written as extra haplotypes when the linearization is not decisive.
    pub max_alternatives: usize,
}

/// A node of the blunt graph: `owner[start..end]` read forward. Every base of every segment
/// lies in exactly one piece, and bases an overlap declares the same are one base.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Piece {
    pub(crate) owner: usize,
    pub(crate) start: usize,
    pub(crate) end: usize,
}

pub(crate) struct Blunt {
    /// Canonical links (a, b) with their overlap, one per adjacency.
    pub(crate) links: Vec<((Oriented, Oriented), usize)>,
    pub(crate) index: FxHashMap<(Oriented, Oriented), usize>,
    pub(crate) pieces: Vec<Piece>,
    /// Each segment, forward, as pieces: (piece, forward?, start offset in the segment).
    pub(crate) walks: Vec<Vec<(usize, bool, usize)>>,
    /// Overlap bases whose two copies differ (an `M` may hold mismatches); the owner's is kept.
    mismatches: usize,
}

pub(crate) fn canonical(a: Oriented, b: Oriented) -> (Oriented, Oriented) {
    let rc = (flip(b), flip(a));
    if (a, b) <= rc {
        (a, b)
    } else {
        rc
    }
}

/// Union-find over bases with a parity bit: a base and the base it is joined to are the same
/// base read on the same strand (0) or on opposite strands (1).
pub(crate) struct StrandUf {
    parent: Vec<u32>,
    parity: Vec<u8>,
}

impl StrandUf {
    pub(crate) fn new(n: usize) -> Self {
        Self {
            parent: (0..n as u32).collect(),
            parity: vec![0; n],
        }
    }

    /// Root of `x` and the strand of `x` relative to it.
    pub(crate) fn find(&mut self, x: usize) -> (usize, u8) {
        let mut path = Vec::new();
        let mut cur = x;
        while self.parent[cur] as usize != cur {
            path.push(cur);
            cur = self.parent[cur] as usize;
        }
        let root = cur;
        // strands to the root, from the node nearest the root outwards
        let mut to_root = 0u8;
        for &node in path.iter().rev() {
            to_root ^= self.parity[node];
            self.parity[node] = to_root;
            self.parent[node] = root as u32;
        }
        (root, if x == root { 0 } else { self.parity[x] })
    }

    /// Declare `x` and `y` the same base, `y` on strand `p` relative to `x`. Returns false when
    /// they already are the same base on the other strand (a base equal to its own complement).
    pub(crate) fn union(&mut self, x: usize, y: usize, p: u8) -> bool {
        let (rx, px) = self.find(x);
        let (ry, py) = self.find(y);
        if rx == ry {
            return px ^ py == p;
        }
        self.parent[ry] = rx as u32;
        self.parity[ry] = px ^ py ^ p;
        true
    }
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        other => other,
    }
}

/// Blunt a graph with overlaps by identifying bases, not by trimming segments.
///
/// A link's overlap says the last `ov` bases of one segment and the first `ov` of the next are
/// the same bases: they are joined (union-find, with the strand), and the joined bases form a
/// graph whose edges are the adjacencies inside segments and the 0M links. Its non-branching
/// runs, cut also at every segment end, are the pieces; each segment is a walk over whole
/// pieces and every path over segments is one over pieces, skipping each entry's overlap.
/// Nothing is dropped or duplicated whatever the overlaps: a short segment whose two overlaps
/// cover it more than once just shares its bases with both neighbours (trimming one side of
/// each link could not express that).
fn blunt(g: &Graph) -> Result<Blunt> {
    blunt_with(g, |_, _| Ok(()))
}

/// `blunt` with more identities: after the links' overlaps, `extra` may join further bases
/// (union-find over all bases of all segments, `base0[s]` the first base of segment `s`).
/// `ovasm pan` joins the bases samples share this way.
pub(crate) fn blunt_with<F>(g: &Graph, extra: F) -> Result<Blunt>
where
    F: FnOnce(&mut StrandUf, &[usize]) -> Result<()>,
{
    let mut links = Vec::new();
    let mut index: FxHashMap<(Oriented, Oriented), usize> = FxHashMap::default();
    for (&(a, b), &ov) in g.links.iter().zip(&g.overlaps) {
        let key = canonical(a, b);
        match index.get(&key) {
            Some(&i) => {
                let (_, prev): ((Oriented, Oriented), usize) = links[i];
                if prev != ov {
                    bail!(
                        "link {} -> {} is listed twice with overlaps {prev} and {ov}",
                        g.label(a),
                        g.label(b)
                    );
                }
            }
            None => {
                if ov > g.seqs[a.0].len() || ov > g.seqs[b.0].len() {
                    bail!(
                        "link {} -> {} overlaps {ov} bp, more than a segment's length",
                        g.label(a),
                        g.label(b)
                    );
                }
                index.insert(key, links.len());
                links.push((key, ov));
            }
        }
    }
    let n = g.seqs.len();
    let mut base0 = Vec::with_capacity(n + 1);
    base0.push(0usize);
    for s in &g.seqs {
        base0.push(base0[base0.len() - 1] + s.len());
    }
    // forward offset of the base at oriented offset `i` of `o`
    let fwd = |o: Oriented, i: usize| {
        if o.1 {
            i
        } else {
            g.seqs[o.0].len() - 1 - i
        }
    };
    let mut uf = StrandUf::new(base0[n]);
    let mut mismatches = 0;
    for &((a, c), ov) in &links {
        let la = g.seqs[a.0].len();
        for t in 0..ov {
            let (fa, fc) = (fwd(a, la - ov + t), fwd(c, t));
            let read = |o: Oriented, f: usize| {
                let b = g.seqs[o.0][f];
                if o.1 {
                    b.to_ascii_uppercase()
                } else {
                    complement(b).to_ascii_uppercase()
                }
            };
            if read(a, fa) != read(c, fc) {
                mismatches += 1;
            }
            let p = u8::from(a.1 != c.1);
            if !uf.union(base0[a.0] + fa, base0[c.0] + fc, p) {
                bail!(
                    "link {} -> {} makes a base of {} its own reverse complement (an odd \
                     palindrome in a hairpin overlap)",
                    g.label(a),
                    g.label(c),
                    g.names[a.0]
                );
            }
        }
    }
    extra(&mut uf, &base0)?;
    // oriented base-class nodes: 2 * root + strand
    let nodes: Vec<Vec<u64>> = (0..n)
        .map(|s| {
            (0..g.seqs[s].len())
                .map(|i| {
                    let (r, p) = uf.find(base0[s] + i);
                    2 * r as u64 + p as u64
                })
                .collect()
        })
        .collect();
    let oriented_nodes = |o: Oriented| -> Vec<u64> {
        if o.1 {
            nodes[o.0].clone()
        } else {
            nodes[o.0].iter().rev().map(|&v| v ^ 1).collect()
        }
    };
    fn add_edge(succ: &mut FxHashMap<u64, Vec<u64>>, u: u64, v: u64) {
        for (x, y) in [(u, v), (v ^ 1, u ^ 1)] {
            let e = succ.entry(x).or_default();
            if !e.contains(&y) {
                e.push(y);
            }
        }
    }
    let mut succ: FxHashMap<u64, Vec<u64>> = FxHashMap::default();
    for v in &nodes {
        for w in v.windows(2) {
            add_edge(&mut succ, w[0], w[1]);
        }
    }
    for &((a, c), ov) in &links {
        if ov == 0 {
            if let (Some(&u), Some(&v)) = (oriented_nodes(a).last(), oriented_nodes(c).first()) {
                add_edge(&mut succ, u, v);
            }
        }
    }
    // a run must stop after every segment's last base and before its first
    let mut ends: FxHashSet<u64> = FxHashSet::default();
    for v in &nodes {
        if let (Some(&first), Some(&last)) = (v.first(), v.last()) {
            ends.insert(last);
            ends.insert(first ^ 1);
        }
    }
    let outdeg = |u: u64| succ.get(&u).map_or(0, Vec::len);
    let cut_after = |u: u64| -> bool {
        if ends.contains(&u) || outdeg(u) != 1 {
            return true;
        }
        let v = succ[&u][0];
        // a self-loop, or a turn onto its own reverse strand (a palindrome's centre: the
        // palindrome is one piece and a hairpin link, not both strands in one piece)
        v == u || v == u ^ 1 || ends.contains(&(v ^ 1)) || outdeg(v ^ 1) != 1
    };

    // pieces: runs between cuts; a run and its reverse are one piece
    let mut pieces: Vec<Piece> = Vec::new();
    let mut by_first: FxHashMap<u64, (usize, bool)> = FxHashMap::default();
    let mut walks: Vec<Vec<(usize, bool, usize)>> = Vec::with_capacity(n);
    for (s, v) in nodes.iter().enumerate() {
        let mut walk = Vec::new();
        let mut start = 0;
        for i in 0..v.len() {
            if i + 1 < v.len() && !cut_after(v[i]) {
                continue;
            }
            let (first, last) = (v[start], v[i]);
            let (id, forward) = match by_first.get(&first) {
                Some(&hit) => hit,
                None => {
                    let id = pieces.len();
                    pieces.push(Piece {
                        owner: s,
                        start,
                        end: i + 1,
                    });
                    by_first.insert(first, (id, true));
                    by_first.entry(last ^ 1).or_insert((id, false));
                    (id, true)
                }
            };
            let p = pieces[id];
            if p.end - p.start != i + 1 - start {
                bail!(
                    "internal: segment {} splits into a run of {} bp where its piece has {} bp",
                    g.names[s],
                    i + 1 - start,
                    p.end - p.start
                );
            }
            walk.push((id, forward, start));
            start = i + 1;
        }
        walks.push(walk);
    }
    Ok(Blunt {
        links,
        index,
        pieces,
        walks,
        mismatches,
    })
}

impl Blunt {
    pub(crate) fn piece_len(&self, id: usize) -> usize {
        self.pieces[id].end - self.pieces[id].start
    }

    /// Pieces of `o` in its orientation: (piece, forward?, start, end) in oriented offsets.
    fn oriented(&self, g: &Graph, o: Oriented) -> Vec<(usize, bool, usize, usize)> {
        let len = g.seqs[o.0].len();
        let mut v: Vec<(usize, bool, usize, usize)> = self.walks[o.0]
            .iter()
            .map(|&(p, f, st)| (p, f, st, st + self.piece_len(p)))
            .collect();
        if !o.1 {
            v.reverse();
            for x in &mut v {
                *x = (x.0, !x.1, len - x.3, len - x.2);
            }
        }
        v
    }

    fn overlap(&self, g: &Graph, prev: Oriented, next: Oriented) -> Result<usize> {
        let &i = self
            .index
            .get(&canonical(prev, next))
            .with_context(|| format!("no link {} -> {}", g.label(prev), g.label(next)))?;
        Ok(self.links[i].1)
    }

    /// Pieces of `o` from oriented offset `entry` on (the bases before it are the previous
    /// segment's, through the overlap).
    pub(crate) fn after(&self, g: &Graph, o: Oriented, entry: usize) -> Result<Vec<(usize, bool)>> {
        let v = self.oriented(g, o);
        if entry < g.seqs[o.0].len() && !v.iter().any(|x| x.2 == entry) {
            bail!("internal: no piece of {} starts at {entry}", g.label(o));
        }
        Ok(v.into_iter()
            .filter(|x| x.2 >= entry)
            .map(|x| (x.0, x.1))
            .collect())
    }
}

/// Piece walk of a segment walk: every segment from its entry overlap on (for a circular walk,
/// the first one too, after the closing link).
pub(crate) fn piece_walk(
    g: &Graph,
    b: &Blunt,
    walk: &[Oriented],
    circular: bool,
) -> Result<Vec<(usize, bool)>> {
    let n = walk.len();
    let mut out = Vec::new();
    for k in 0..n {
        let entry = if k > 0 {
            b.overlap(g, walk[k - 1], walk[k])?
        } else if circular {
            b.overlap(g, walk[n - 1], walk[0])?
        } else {
            0
        };
        out.extend(b.after(g, walk[k], entry)?);
    }
    Ok(out)
}

/// Mean depth of a segment from whichever tag the backend wrote, and the rule used.
pub(crate) fn depth(tags: &FxHashMap<String, String>, len: usize) -> Option<(f64, &'static str)> {
    let num = |t: &str| tags.get(t).and_then(|v| v.parse::<f64>().ok());
    let per_base = |t: &str| num(t).map(|v| v / len.max(1) as f64);
    num("dp")
        .map(|v| (v, "dp"))
        .or_else(|| num("DP").map(|v| (v, "DP")))
        .or_else(|| per_base("KC").map(|v| (v, "KC/LN")))
        .or_else(|| per_base("RC").map(|v| (v, "RC/LN")))
}

/// Tags of every S line, by segment name (`XX:T:value` -> XX -> value).
pub(crate) fn segment_tags(text: &str) -> FxHashMap<String, FxHashMap<String, String>> {
    let mut out = FxHashMap::default();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.first() != Some(&"S") || f.len() < 3 {
            continue;
        }
        let tags = f[3..]
            .iter()
            .filter_map(|t| {
                let mut p = t.splitn(3, ':');
                let key = p.next()?;
                p.next()?; // type
                Some((key.to_string(), p.next()?.to_string()))
            })
            .collect();
        out.insert(f[1].to_string(), tags);
    }
    out
}

pub(crate) fn parse_walk(g: &Graph, path: &str) -> Result<Vec<Oriented>> {
    let index: FxHashMap<&str, usize> = g
        .names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    path.split_whitespace()
        .map(|l| {
            let (name, sign) = l.split_at(l.len().saturating_sub(1));
            let s = *index
                .get(name)
                .with_context(|| format!("path names unknown segment {name:?}"))?;
            match sign {
                "+" => Ok((s, true)),
                "-" => Ok((s, false)),
                _ => bail!("path step {l:?} must end in + or -"),
            }
        })
        .collect()
}

pub(crate) fn valid_pansn_field(s: &str) -> bool {
    !s.is_empty()
        && s.bytes()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, b'_' | b'.' | b'-'))
}

#[derive(Serialize)]
pub struct UnifyReport {
    pub format: &'static str,
    pub source: String,
    pub sample: String,
    pub backend: String,
    pub segments_in: usize,
    pub links_in: usize,
    pub pieces_out: usize,
    pub links_out: usize,
    /// Segments made of several pieces.
    pub segments_cut: usize,
    /// Bases that were a second copy of another through an overlap.
    pub overlap_bp_removed: usize,
    /// Pieces in more than one segment (sequence the overlaps shared).
    pub shared_pieces: usize,
    /// Overlap bases whose two copies differed; the owner segment's base was kept.
    pub overlap_mismatches: usize,
    /// Depth rule per source tag (e.g. `KC/LN`), and segments without any depth tag.
    pub depth_sources: Vec<String>,
    pub segments_without_depth: usize,
    pub links_with_evidence: usize,
    pub paths: Vec<String>,
}

pub fn run(
    gfa: &Path,
    evidence: Option<&Path>,
    linearization: Option<&Path>,
    p: &UnifyParams,
    out: &Path,
    out_json: &Path,
) -> Result<UnifyReport> {
    for (what, v) in [("sample", p.sample), ("molecule", p.molecule)] {
        if !valid_pansn_field(v) {
            bail!("{what} {v:?} must use only letters, digits, '_', '.' and '-'");
        }
    }
    let text = fs::read_to_string(gfa).with_context(|| format!("cannot read {}", gfa.display()))?;
    let g = parse_gfa(gfa)?;
    let b = blunt(&g)?;
    let tags = segment_tags(&text);
    let read_json = |path: &Path| -> Result<Value> {
        Ok(serde_json::from_str(
            &fs::read_to_string(path).with_context(|| format!("cannot read {}", path.display()))?,
        )?)
    };

    // Read support per canonical link.
    let mut support: FxHashMap<usize, u64> = FxHashMap::default();
    if let Some(ev) = evidence {
        let ev = read_json(ev)?;
        for l in ev["links"].as_array().context("evidence has no links")? {
            let walk = parse_walk(
                &g,
                &format!(
                    "{} {}",
                    l["from"].as_str().context("link without from")?,
                    l["to"].as_str().context("link without to")?
                ),
            )?;
            let &i = b
                .index
                .get(&canonical(walk[0], walk[1]))
                .with_context(|| format!("evidence link {:?} is not in the GFA", l))?;
            let e = support.entry(i).or_default();
            *e = (*e).max(l["reads"].as_u64().unwrap_or(0));
        }
    }

    // names: `<sample>.<owner>` or, when the owner segment owns several, `<sample>.<owner>.<k>`
    let mut owned = vec![0usize; g.seqs.len()];
    let mut rank = Vec::with_capacity(b.pieces.len());
    for pc in &b.pieces {
        owned[pc.owner] += 1;
        rank.push(owned[pc.owner]);
    }
    let name = |piece: usize| -> String {
        let o = b.pieces[piece].owner;
        if owned[o] == 1 {
            format!("{}.{}", p.sample, g.names[o])
        } else {
            format!("{}.{}.{}", p.sample, g.names[o], rank[piece])
        }
    };
    let sign = |f: bool| if f { '+' } else { '-' };

    let mut gfa_out = String::new();
    let mut depth_sources = BTreeSet::new();
    let dps: Vec<Option<(f64, &str)>> = g
        .seqs
        .iter()
        .enumerate()
        .map(|(s, seq)| tags.get(&g.names[s]).and_then(|t| depth(t, seq.len())))
        .collect();
    let without_depth = dps.iter().filter(|d| d.is_none()).count();
    for (_, rule) in dps.iter().flatten() {
        depth_sources.insert(*rule);
    }
    writeln!(
        gfa_out,
        "H\tVN:Z:1.0\tov:Z:{FORMAT}\tsn:Z:{}\tbk:Z:{}\tor:Z:{}",
        p.sample, p.backend, p.organelle
    )?;
    for (id, pc) in b.pieces.iter().enumerate() {
        let piece = &g.seqs[pc.owner][pc.start..pc.end];
        write!(
            gfa_out,
            "S\t{}\t{}\tLN:i:{}",
            name(id),
            std::str::from_utf8(piece)?,
            piece.len()
        )?;
        if let Some((d, _)) = dps[pc.owner] {
            write!(gfa_out, "\tdp:f:{d:.3}")?;
        }
        gfa_out.push('\n');
    }
    // Links between pieces: consecutive pieces of a segment, and 0M links. An overlapping link
    // is already the step inside the entered segment from the shared bases onwards.
    type PieceEnd = (usize, bool);
    let canon = |x: PieceEnd, y: PieceEnd| {
        let rc = ((y.0, !y.1), (x.0, !x.1));
        if (x, y) <= rc {
            (x, y)
        } else {
            rc
        }
    };
    let mut piece_links: BTreeSet<(PieceEnd, PieceEnd)> = BTreeSet::new();
    for w in &b.walks {
        for pair in w.windows(2) {
            piece_links.insert(canon((pair[0].0, pair[0].1), (pair[1].0, pair[1].1)));
        }
    }
    let mut ev_of: FxHashMap<(PieceEnd, PieceEnd), u64> = FxHashMap::default();
    for (i, &((a, c), ov)) in b.links.iter().enumerate() {
        let (Some(&from), Some(&to)) = (b.after(&g, a, 0)?.last(), b.after(&g, c, ov)?.first())
        else {
            continue; // `c` lies wholly inside the overlap
        };
        let key = canon(from, to);
        if ov == 0 {
            piece_links.insert(key);
        }
        if let Some(&n) = support.get(&i) {
            let e = ev_of.entry(key).or_default();
            *e = (*e).max(n);
        }
    }
    let mut links_out = 0;
    for &(x, y) in &piece_links {
        write!(
            gfa_out,
            "L\t{}\t{}\t{}\t{}\t0M",
            name(x.0),
            sign(x.1),
            name(y.0),
            sign(y.1)
        )?;
        if let Some(n) = ev_of.get(&(x, y)) {
            write!(gfa_out, "\tev:i:{n}")?;
        }
        gfa_out.push('\n');
        links_out += 1;
    }

    // Paths: every molecule of the unringing, then near-tied alternatives when not decisive.
    let mut paths = Vec::new();
    if let Some(lin) = linearization {
        let lin = read_json(lin)?;
        let near = |best: &Value, a: &Value| {
            best["score"].as_f64().unwrap_or(f64::NEG_INFINITY)
                - a["score"].as_f64().unwrap_or(f64::NEG_INFINITY)
                < std::f64::consts::LN_2
        };
        // (haplotype, molecule number, molecule) in output order
        let mut chosen: Vec<(usize, usize, &Value)> = Vec::new();
        match lin["components"].as_array() {
            // unringing: each component's best cover is one set of molecules; its near-tied
            // alternative covers (same molecule count) are further haplotypes of them
            Some(comps) => {
                let mut offset = 0;
                for c in comps.iter().filter(|c| c["best"].is_object()) {
                    let mut covers = vec![&c["best"]];
                    if c["decisive"] == Value::Bool(false) {
                        covers.extend(
                            c["alternatives"]
                                .as_array()
                                .into_iter()
                                .flatten()
                                .filter(|a| near(&c["best"], a))
                                .take(p.max_alternatives),
                        );
                    }
                    for (hap, cover) in covers.iter().enumerate() {
                        for (j, m) in cover["molecules"]
                            .as_array()
                            .into_iter()
                            .flatten()
                            .enumerate()
                        {
                            chosen.push((hap, offset + j, m));
                        }
                    }
                    offset += c["best"]["molecules"].as_array().map_or(0, Vec::len);
                }
            }
            // a report from before unringing: one path and its alternatives
            None => {
                if lin["best"].is_object() {
                    chosen.push((0, 0, &lin["best"]));
                    if lin["decisive"] == Value::Bool(false) {
                        let alts = lin["alternatives"]
                            .as_array()
                            .into_iter()
                            .flatten()
                            .filter(|a| near(&lin["best"], a))
                            .take(p.max_alternatives);
                        chosen.extend(alts.enumerate().map(|(h, a)| (h + 1, 0, a)));
                    }
                }
            }
        }
        let several = chosen.iter().any(|&(_, m, _)| m > 0);
        for &(hap, mol, c) in &chosen {
            let walk = parse_walk(
                &g,
                c["path"].as_str().context("linearization without a path")?,
            )?;
            let circular = c["circular"].as_bool().unwrap_or(false);
            let steps = piece_walk(&g, &b, &walk, circular)?;
            // PanSN sample#haplotype#contig; several molecules are numbered contigs
            let contig = if several {
                format!("{}.{}", p.molecule, mol + 1)
            } else {
                p.molecule.to_string()
            };
            let path_name = format!("{}#{}#{}", p.sample, hap + 1, contig);
            let names: Vec<String> = steps
                .iter()
                .map(|&(id, f)| format!("{}{}", name(id), sign(f)))
                .collect();
            let cigar = vec!["0M"; steps.len().saturating_sub(1)].join(",");
            write!(
                gfa_out,
                "P\t{path_name}\t{}\t{}",
                names.join(","),
                if cigar.is_empty() { "*" } else { &cigar }
            )?;
            if circular {
                gfa_out.push_str("\tTP:Z:circular");
            }
            gfa_out.push('\n');
            paths.push(path_name);
        }
    }

    fs::create_dir_all(out.parent().unwrap_or(Path::new(".")))?;
    fs::write(out, &gfa_out)?;
    let report = UnifyReport {
        format: FORMAT,
        source: gfa.display().to_string(),
        sample: p.sample.to_string(),
        backend: p.backend.to_string(),
        segments_in: g.seqs.len(),
        links_in: b.links.len(),
        pieces_out: b.pieces.len(),
        links_out,
        segments_cut: b.walks.iter().filter(|w| w.len() > 1).count(),
        overlap_bp_removed: g.seqs.iter().map(Vec::len).sum::<usize>()
            - (0..b.pieces.len()).map(|i| b.piece_len(i)).sum::<usize>(),
        shared_pieces: {
            let mut uses = vec![0usize; b.pieces.len()];
            for w in &b.walks {
                let mut seen: Vec<usize> = w.iter().map(|x| x.0).collect();
                seen.sort_unstable();
                seen.dedup();
                for p in seen {
                    uses[p] += 1;
                }
            }
            uses.iter().filter(|&&u| u > 1).count()
        },
        overlap_mismatches: b.mismatches,
        depth_sources: depth_sources.into_iter().map(String::from).collect(),
        segments_without_depth: without_depth,
        links_with_evidence: support.len(),
        paths,
    };
    fs::create_dir_all(out_json.parent().unwrap_or(Path::new(".")))?;
    fs::write(out_json, serde_json::to_string_pretty(&report)?)?;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::evidence::revcomp;

    fn seq(n: usize, seed: u64) -> Vec<u8> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(s >> 62) as usize]
            })
            .collect()
    }

    /// Spell a piece walk from the blunt pieces.
    fn spell(g: &Graph, b: &Blunt, steps: &[(usize, bool)]) -> Vec<u8> {
        steps
            .iter()
            .flat_map(|&(id, f)| {
                let pc = b.pieces[id];
                let piece = g.seqs[pc.owner][pc.start..pc.end].to_vec();
                if f {
                    piece
                } else {
                    revcomp(&piece)
                }
            })
            .collect()
    }

    fn total(b: &Blunt) -> usize {
        (0..b.pieces.len()).map(|i| b.piece_len(i)).sum()
    }

    /// Every walk spells what the overlapping segments spell, in both directions.
    fn assert_walks(g: &Graph, b: &Blunt, walks: &[Vec<Oriented>]) {
        for walk in walks {
            let back: Vec<Oriented> = walk.iter().rev().map(|&o| flip(o)).collect();
            for w in [walk.clone(), back] {
                let steps = piece_walk(g, b, &w, false).unwrap();
                assert_eq!(spell(g, b, &steps), spell_overlapping(g, &w), "{w:?}");
            }
        }
    }

    /// Spell a segment walk directly with its overlaps.
    fn spell_overlapping(g: &Graph, walk: &[Oriented]) -> Vec<u8> {
        let o = |x: Oriented| {
            if x.1 {
                g.seqs[x.0].clone()
            } else {
                revcomp(&g.seqs[x.0])
            }
        };
        let mut out = o(walk[0]);
        for w in walk.windows(2) {
            let key = canonical(w[0], w[1]);
            let ov = g
                .links
                .iter()
                .zip(&g.overlaps)
                .find(|(&(a, c), _)| canonical(a, c) == key)
                .map(|(_, &v)| v)
                .unwrap();
            out.extend_from_slice(&o(w[1])[ov..]);
        }
        out
    }

    /// A -> X (60M), B -> X (40M), X -> C (30M), plus X- -> D (25M) entering D reversed.
    fn overlapping() -> Graph {
        let x = seq(500, 1);
        let a = [seq(300, 2), x[..60].to_vec()].concat();
        let bb = [seq(300, 3), x[..40].to_vec()].concat();
        let c = [x[x.len() - 30..].to_vec(), seq(300, 4)].concat();
        // X- continues into D-: the end of X- (= revcomp of X's start) overlaps D- 's start.
        let x_rc = revcomp(&x);
        let d_minus = [x_rc[x_rc.len() - 25..].to_vec(), seq(300, 5)].concat();
        Graph {
            names: ["A", "B", "X", "C", "D"].map(String::from).to_vec(),
            seqs: vec![a, bb, x, c, revcomp(&d_minus)],
            links: vec![
                ((0, true), (2, true)),
                ((1, true), (2, true)),
                ((2, true), (3, true)),
                ((2, false), (4, false)),
            ],
            overlaps: vec![60, 40, 30, 25],
        }
    }

    #[test]
    fn blunt_walks_spell_the_same_sequence_without_duplicating_bases() {
        let g = overlapping();
        let b = blunt(&g).unwrap();
        for walk in [
            vec![(0, true), (2, true), (3, true)],
            vec![(1, true), (2, true), (3, true)],
            vec![(3, false), (2, false), (0, false)],
            vec![(0, true), (2, true)],
            vec![(2, false), (4, false)],
        ] {
            let steps = piece_walk(&g, &b, &walk, false).unwrap();
            assert_eq!(
                spell(&g, &b, &steps),
                spell_overlapping(&g, &walk),
                "{walk:?}"
            );
        }
        let total_in: usize = g.seqs.iter().map(Vec::len).sum();
        assert_eq!(
            total(&b),
            total_in - (60 + 40 + 30 + 25),
            "each overlapping base is written once"
        );
        assert_eq!(b.mismatches, 0);
    }

    #[test]
    fn a_segment_shorter_than_its_two_overlaps_shares_its_bases_with_both_neighbours() {
        // S (10 bp) overlaps N by 7 at its start and M by 6 at its end: 13 > 10, so no side of
        // each link could be trimmed alone (a homopolymer stutter pair in ONT graphs).
        let sq = seq(10, 9);
        let n = [seq(300, 10), sq[..7].to_vec()].concat();
        let m = [sq[4..].to_vec(), seq(300, 11)].concat();
        let g = Graph {
            names: ["N", "S", "M"].map(String::from).to_vec(),
            seqs: vec![n, sq, m],
            links: vec![((0, true), (1, true)), ((1, true), (2, true))],
            overlaps: vec![7, 6],
        };
        let b = blunt(&g).unwrap();
        assert_walks(
            &g,
            &b,
            &[
                vec![(0, true), (1, true), (2, true)],
                vec![(1, true), (2, true)],
            ],
        );
        assert_eq!(total(&b), 307 + 10 + 306 - 7 - 6);
    }

    #[test]
    fn a_hairpin_overlap_joins_a_palindrome_to_itself() {
        // X ends in a 50 bp palindrome and continues into its own reverse complement.
        let h = seq(25, 12);
        let x = [seq(30, 13), h.clone(), revcomp(&h)].concat();
        let g = Graph {
            names: vec!["X".into()],
            seqs: vec![x],
            links: vec![((0, true), (0, false))],
            overlaps: vec![50],
        };
        let b = blunt(&g).unwrap();
        assert_walks(&g, &b, &[vec![(0, true), (0, false)]]);
        assert_eq!(
            total(&b),
            80 - 25,
            "the palindrome's halves are one sequence"
        );
    }

    #[test]
    fn a_base_that_must_be_its_own_complement_is_rejected() {
        // An odd palindrome joined to itself pairs its centre base with its own complement.
        let h = seq(25, 14);
        let x = [seq(29, 15), h.clone(), b"A".to_vec(), revcomp(&h)].concat();
        let g = Graph {
            names: vec!["X".into()],
            seqs: vec![x],
            links: vec![((0, true), (0, false))],
            overlaps: vec![51],
        };
        let err = blunt(&g).err().map(|e| e.to_string()).unwrap_or_default();
        assert!(err.contains("its own reverse complement"), "{err}");
    }

    #[test]
    fn a_circular_path_closes_through_its_trimmed_link() {
        // One segment closing on itself with a 50 bp overlap: the circle is len - 50.
        let core = seq(400, 7);
        let s = [core.clone(), core[..50].to_vec()].concat();
        let g = Graph {
            names: vec!["S".into()],
            seqs: vec![s],
            links: vec![((0, true), (0, true))],
            overlaps: vec![50],
        };
        let b = blunt(&g).unwrap();
        let steps = piece_walk(&g, &b, &[(0, true)], true).unwrap();
        let spelled = spell(&g, &b, &steps);
        assert_eq!(spelled.len(), 400);
        let doubled = [core.clone(), core.clone()].concat();
        assert!(
            doubled.windows(400).any(|w| w == spelled.as_slice()),
            "a rotation of the circle"
        );
    }

    #[test]
    fn duplicate_link_lines_must_agree_on_their_overlap() {
        let mut g = overlapping();
        g.links.push(((3, false), (2, false))); // same adjacency as X+ -> C+
        g.overlaps.push(30);
        assert!(blunt(&g).is_ok());
        g.overlaps[4] = 31;
        assert!(blunt(&g).is_err());
    }

    #[test]
    fn depth_follows_the_backend_tag() {
        let t = |pairs: &[(&str, &str)]| {
            pairs
                .iter()
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect::<FxHashMap<_, _>>()
        };
        assert_eq!(depth(&t(&[("dp", "44")]), 100), Some((44.0, "dp")));
        assert_eq!(
            depth(&t(&[("KC", "1000"), ("LN", "100")]), 100),
            Some((10.0, "KC/LN"))
        );
        assert_eq!(depth(&t(&[("RC", "500")]), 100), Some((5.0, "RC/LN")));
        assert_eq!(depth(&t(&[("SC", "1.0")]), 100), None);
    }

    #[test]
    fn pansn_fields_are_restricted() {
        assert!(valid_pansn_field("NIP_hifi.v1-2"));
        assert!(!valid_pansn_field("a#b"));
        assert!(!valid_pansn_field(""));
    }
}
