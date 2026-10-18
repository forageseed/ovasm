//! Optional, graph-bound consumers of U1 evidence. Missing evidence is not evidence of absence.
use crate::evidence::{flip, parse_gfa, Graph, Oriented};
use crate::unzip::graph_fingerprint;
use anyhow::{bail, Context, Result};
use serde::Serialize;
use serde_json::Value;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

type Triple = (Oriented, Oriented, Oriented);

#[derive(Default)]
pub(crate) struct Constraints {
    blocked: BTreeSet<Triple>,
    pub active_nodes: usize,
    pub phases: BTreeMap<Triple, u64>,
    pub resolved: BTreeSet<usize>,
}

impl Constraints {
    pub fn allows(&self, a: Oriented, s: Oriented, b: Oriented) -> bool {
        !self.blocked.contains(&(a, s, b))
    }
    pub fn allows_walk(&self, walk: &[Oriented], circular: bool) -> bool {
        if !walk.windows(3).all(|w| self.allows(w[0], w[1], w[2])) {
            return false;
        }
        if circular && !walk.is_empty() {
            let n = walk.len();
            self.allows(walk[n - 1], walk[0], walk[1 % n])
                && self.allows(walk[(n + n - 2) % n], walk[n - 1], walk[0])
        } else {
            true
        }
    }
    pub fn load(g: &Graph, path: &Path) -> Result<Self> {
        let v: Value = serde_json::from_str(&fs::read_to_string(path)?)?;
        if v["graph_fingerprint"].as_str() != Some(graph_fingerprint(g).as_str()) {
            bail!("phase report does not match graph sequences, links and overlaps; regenerate unzip evidence");
        }
        let min = v["min_count"]
            .as_u64()
            .filter(|&n| n > 0)
            .context("invalid phase min_count")?;
        let unique = v["unique_records"]
            .as_u64()
            .context("missing unique_records")?;
        let index: BTreeMap<_, _> = g
            .names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.as_str(), i))
            .collect();
        let label = |text: &str| -> Result<Oriented> {
            let (name, sign) = text.split_at(
                text.len()
                    .checked_sub(1)
                    .context("empty oriented segment")?,
            );
            if sign != "+" && sign != "-" {
                bail!("invalid oriented segment {text}");
            }
            Ok((
                *index.get(name).context("unknown phase segment")?,
                sign == "+",
            ))
        };
        let mut pred: BTreeMap<_, BTreeSet<_>> = BTreeMap::new();
        let mut succ: BTreeMap<_, BTreeSet<_>> = BTreeMap::new();
        for &(a, b) in &g.links {
            for (x, y) in [(a, b), (flip(b), flip(a))] {
                succ.entry(x).or_default().insert(y);
                pred.entry(y).or_default().insert(x);
            }
        }
        let expected: BTreeSet<_> = (0..g.names.len())
            .filter(|&s| {
                pred.get(&(s, true)).map_or(0, BTreeSet::len) >= 2
                    || succ.get(&(s, true)).map_or(0, BTreeSet::len) >= 2
            })
            .collect();
        let mut seen = BTreeSet::new();
        let mut out = Self::default();
        for row in v["repeat_nodes"]
            .as_array()
            .context("missing repeat_nodes")?
        {
            let name = row["segment"].as_str().context("missing repeat segment")?;
            let s = (*index.get(name).context("unknown repeat")?, true);
            if !expected.contains(&s.0) || !seen.insert(s.0) {
                bail!("duplicate or nonbranching phase node");
            }
            let entries = pred.get(&s).cloned().unwrap_or_default();
            let exits = succ.get(&s).cloned().unwrap_or_default();
            let mut adj = BTreeMap::new();
            for (field, ends) in [("entries", &entries), ("exits", &exits)] {
                let mut actual = BTreeSet::new();
                for end in row[field].as_array().context("missing phase ends")? {
                    let o = label(end["segment"].as_str().context("invalid end label")?)?;
                    if !actual.insert(o) {
                        bail!("duplicate phase end");
                    }
                    let n = end["adjacency_reads"]
                        .as_u64()
                        .filter(|&n| n <= unique)
                        .context("invalid adjacency count")?;
                    adj.insert((field, o), n);
                }
                if &actual != ends {
                    bail!("phase ends differ from graph");
                }
            }
            let mut pairs = BTreeMap::new();
            for p in row["phases"].as_array().context("missing phases")? {
                let a = label(p["entry"].as_str().context("invalid phase entry")?)?;
                let b = label(p["exit"].as_str().context("invalid phase exit")?)?;
                let n = p["reads"]
                    .as_u64()
                    .filter(|&n| n > 0 && n <= unique)
                    .context("invalid phase count")?;
                if !entries.contains(&a)
                    || !exits.contains(&b)
                    || n > adj[&("entries", a)]
                    || n > adj[&("exits", b)]
                    || p["supported"].as_bool() != Some(n >= min)
                    || pairs.insert((a, b), n).is_some()
                {
                    bail!("inconsistent phase evidence");
                }
                out.phases.insert((a, s, b), n);
            }
            let supported: Vec<_> = pairs
                .iter()
                .filter(|(_, n)| **n >= min)
                .map(|(&p, _)| p)
                .collect();
            let covered_entries: BTreeSet<_> = supported.iter().map(|p| p.0).collect();
            let covered_exits: BTreeSet<_> = supported.iter().map(|p| p.1).collect();
            // Fully sampled ends are required before treating a zero count as a prohibition.
            // Weak nonzero phases are never prohibited or promoted to supported phases.
            let complete =
                !supported.is_empty() && covered_entries == entries && covered_exits == exits;
            let resolved = complete
                && supported.len() == entries.len()
                && supported.len() == exits.len()
                && pairs.len() == supported.len();
            let class = if supported.is_empty() {
                "unresolved"
            } else if supported.len() > covered_entries.len()
                || supported.len() > covered_exits.len()
            {
                "mixed"
            } else if pairs.len() > supported.len() {
                "conflicted"
            } else if complete {
                "resolved"
            } else {
                "partial"
            };
            if row["classification"].as_str() != Some(class) {
                bail!("inconsistent phase classification");
            }
            if complete {
                out.active_nodes += 1;
                for &a in &entries {
                    for &b in &exits {
                        if !pairs.contains_key(&(a, b)) {
                            out.blocked.insert((a, s, b));
                            out.blocked.insert((flip(b), flip(s), flip(a)));
                        }
                    }
                }
            }
            if resolved {
                out.resolved.insert(s.0);
            }
        }
        if seen != expected {
            bail!("incomplete phase report");
        }
        Ok(out)
    }
}

#[derive(Serialize)]
pub struct ExportReport {
    pub graph_changed: bool,
    pub split_segments: Vec<String>,
    pub supported_paths: usize,
    pub physical_molecule_count: Option<u32>,
}

/// P overlays preserve the source graph. Materialization replaces only fully resolved
/// one-to-one nodes with context copies. Mixed, partial, conflicted and unresolved nodes stay.
pub fn export(
    gfa: &Path,
    evidence: &Path,
    overlay: Option<&Path>,
    materialized: Option<&Path>,
) -> Result<ExportReport> {
    let g = parse_gfa(gfa)?;
    let c = Constraints::load(&g, evidence)?;
    let text = fs::read_to_string(gfa)?;
    let v: Value = serde_json::from_str(&fs::read_to_string(evidence)?)?;
    let min = v["min_count"].as_u64().unwrap();
    for p in overlay.into_iter().chain(materialized) {
        if p == gfa || p == evidence || p.exists() {
            bail!("phase exports require distinct new output paths");
        }
    }
    if overlay.is_some() && overlay == materialized {
        bail!("phase export paths must differ");
    }
    let supported: Vec<_> = c.phases.iter().filter(|(_, n)| **n >= min).collect();
    let mut paths = String::new();
    let old_names: BTreeSet<_> = text
        .lines()
        .filter_map(|l| l.strip_prefix("P\t").and_then(|s| s.split('\t').next()))
        .collect();
    for (i, (&tr, &n)) in supported.iter().enumerate() {
        let name = format!("ovasm.phase.{}", i + 1);
        if old_names.contains(name.as_str()) {
            bail!("phase path name already exists");
        }
        let overlap = |a, b| -> usize {
            g.links
                .iter()
                .zip(&g.overlaps)
                .find_map(|(&(x, y), &ov)| {
                    ((a, b) == (x, y) || (a, b) == (flip(y), flip(x))).then_some(ov)
                })
                .expect("validated phase link")
        };
        let overlaps = [overlap(tr.0, tr.1), overlap(tr.1, tr.2)];
        paths.push_str(&format!(
            "P\t{name}\t{},{},{}\t{}M,{}M\tRC:i:{n}\n",
            g.label(tr.0),
            g.label(tr.1),
            g.label(tr.2),
            overlaps[0],
            overlaps[1]
        ));
    }
    // Self-linked repeats need multi-visit phasing; triples cannot certify their separation.
    let split: BTreeSet<_> = c
        .resolved
        .iter()
        .copied()
        .filter(|s| !g.links.iter().any(|(a, b)| a.0 == *s && b.0 == *s))
        .collect();
    let mut copies: BTreeMap<usize, Vec<(String, Option<(Oriented, Oriented, u64)>)>> =
        BTreeMap::new();
    for s in 0..g.names.len() {
        let mut list = Vec::new();
        if split.contains(&s) {
            for (&(a, _, b), &n) in c.phases.iter().filter(|((_, m, _), _)| m.0 == s) {
                let name = format!("{}.phase{}", g.names[s], list.len() + 1);
                if g.names.contains(&name) {
                    bail!("phase segment name collision");
                }
                list.push((name, Some((a, b, n))));
            }
        } else {
            list.push((g.names[s].clone(), None));
        }
        copies.insert(s, list);
    }
    if let Some(out) = materialized {
        if split.is_empty() {
            fs::write(out, &text)?;
        } else {
            // Existing P/W records can contain copies whose context is now different. Refuse
            // to silently drop or rewrite them; callers can use the unchanged overlay instead.
            if text
                .lines()
                .any(|l| l.starts_with("P\t") || l.starts_with("W\t") || l.starts_with("C\t"))
            {
                bail!("materialization needs an S/L graph without existing P/W/C records");
            }
            let tags: BTreeMap<_, _> = text
                .lines()
                .filter(|l| l.starts_with("S\t"))
                .map(|l| {
                    let f: Vec<_> = l.split('\t').collect();
                    (
                        f[1],
                        f.iter().skip(3).map(|s| s.to_string()).collect::<Vec<_>>(),
                    )
                })
                .collect();
            let mut buf = String::new();
            for l in text
                .lines()
                .filter(|l| !l.starts_with("S\t") && !l.starts_with("L\t"))
            {
                buf.push_str(l);
                buf.push('\n');
            }
            for (s, list) in &copies {
                let total: u64 = list.iter().filter_map(|(_, p)| p.map(|p| p.2)).sum();
                for (name, phase) in list {
                    let mut ts = tags[g.names[*s].as_str()].clone();
                    if let Some((_, _, n)) = phase {
                        ts = ts
                            .into_iter()
                            .map(|tag| {
                                for prefix in ["dp:f:", "dm:f:"] {
                                    if let Some(d) =
                                        tag.strip_prefix(prefix).and_then(|s| s.parse::<f64>().ok())
                                    {
                                        return format!(
                                            "{prefix}{:.6}",
                                            d * *n as f64 / total as f64
                                        );
                                    }
                                }
                                tag
                            })
                            .filter(|t| !t.starts_with("KC:i:"))
                            .collect();
                        ts.push(format!("OR:Z:{}", g.names[*s]));
                        ts.push(format!("RC:i:{n}"));
                    }
                    buf.push_str(&format!("S\t{name}\t{}", std::str::from_utf8(&g.seqs[*s])?));
                    for t in ts {
                        buf.push('\t');
                        buf.push_str(&t);
                    }
                    buf.push('\n');
                }
            }
            for (&(a, b), &ov) in g.links.iter().zip(&g.overlaps) {
                for (an, ap) in &copies[&a.0] {
                    for (bn, bp) in &copies[&b.0] {
                        if ap.is_some_and(|(x, y, _)| if a.1 { y != b } else { flip(x) != b })
                            || bp.is_some_and(|(x, y, _)| if b.1 { x != a } else { flip(y) != a })
                        {
                            continue;
                        }
                        buf.push_str(&format!(
                            "L\t{an}\t{}\t{bn}\t{}\t{ov}M\n",
                            if a.1 { '+' } else { '-' },
                            if b.1 { '+' } else { '-' }
                        ));
                    }
                }
            }
            fs::write(out, buf)?;
        }
    }
    if let Some(out) = overlay {
        fs::write(
            out,
            format!(
                "{}{}{}",
                text,
                if text.ends_with('\n') { "" } else { "\n" },
                paths
            ),
        )?;
    }
    Ok(ExportReport {
        graph_changed: materialized.is_some() && !split.is_empty(),
        split_segments: split.iter().map(|&s| g.names[s].clone()).collect(),
        supported_paths: supported.len(),
        physical_molecule_count: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;
    use std::sync::atomic::{AtomicUsize, Ordering};
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    fn fixture(counts: [u64; 4], class: &str) -> (std::path::PathBuf, std::path::PathBuf) {
        let dir = std::env::temp_dir().join(format!(
            "ovasm-phase-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&dir).unwrap();
        let gfa = dir.join("graph.gfa");
        fs::write(&gfa,"H\tVN:Z:1.0\nS\tA\tAACG\tdp:f:100\nS\tC\tCAGT\tdp:f:100\nS\tR\tGCTA\tdp:f:200\nS\tB\tTTGC\tdp:f:100\nS\tD\tATGC\tdp:f:100\nL\tA\t+\tR\t+\t0M\nL\tC\t+\tR\t+\t0M\nL\tR\t+\tB\t+\t0M\nL\tR\t+\tD\t+\t0M\n").unwrap();
        let g = parse_gfa(&gfa).unwrap();
        let phases: Vec<_> = [("A+", "B+"), ("A+", "D+"), ("C+", "B+"), ("C+", "D+")]
            .into_iter()
            .zip(counts)
            .filter(|(_, n)| *n > 0)
            .map(|((a, b), n)| json!({"entry":a,"exit":b,"reads":n,"supported":n>=80}))
            .collect();
        let report = json!({"graph_fingerprint":graph_fingerprint(&g),"min_count":80,"unique_records":400,
            "repeat_nodes":[{"segment":"R","classification":class,
            "entries":[{"segment":"A+","adjacency_reads":200},{"segment":"C+","adjacency_reads":200}],
            "exits":[{"segment":"B+","adjacency_reads":200},{"segment":"D+","adjacency_reads":200}],"phases":phases}]});
        let ev = dir.join("phases.json");
        fs::write(&ev, report.to_string()).unwrap();
        (gfa, ev)
    }
    #[test]
    fn unresolved_and_incompletely_sampled_ends_do_not_prohibit_paths() {
        for (counts, class) in [([37, 18, 41, 29], "unresolved"), ([80, 0, 0, 0], "partial")] {
            let (p, e) = fixture(counts, class);
            let g = parse_gfa(&p).unwrap();
            let c = Constraints::load(&g, &e).unwrap();
            assert!(c.blocked.is_empty());
            assert!(c.resolved.is_empty());
            let out = p.with_file_name("materialized.gfa");
            assert!(!export(&p, &e, None, Some(&out)).unwrap().graph_changed);
            assert_eq!(fs::read(p).unwrap(), fs::read(out).unwrap());
        }
    }
    #[test]
    fn zero_pair_is_forbidden_on_both_strands_and_across_circular_origin() {
        let (p, e) = fixture([80, 0, 0, 80], "resolved");
        let g = parse_gfa(&p).unwrap();
        let c = Constraints::load(&g, &e).unwrap();
        let (a, r, b, d) = ((0, true), (2, true), (3, true), (4, true));
        assert!(c.allows(a, r, b));
        assert!(!c.allows(a, r, d));
        assert!(!c.allows(flip(d), flip(r), flip(a)));
        assert!(c.allows_walk(&[r, d, a], false));
        assert!(!c.allows_walk(&[r, d, a], true));
    }
    #[test]
    fn weak_counterexample_is_retained_and_mixed_phases_are_not_split() {
        for (counts, class) in [([80, 1, 0, 80], "conflicted"), ([80, 80, 80, 80], "mixed")] {
            let (p, e) = fixture(counts, class);
            let g = parse_gfa(&p).unwrap();
            let c = Constraints::load(&g, &e).unwrap();
            assert!(c.allows((0, true), (2, true), (4, true)));
            assert!(c.resolved.is_empty());
        }
    }
    #[test]
    fn resolved_materialization_preserves_spelling_and_only_certified_connections() {
        let (p, e) = fixture([80, 0, 0, 80], "resolved");
        let overlay = p.with_file_name("overlay.gfa");
        let out = p.with_file_name("split.gfa");
        let result = export(&p, &e, Some(&overlay), Some(&out)).unwrap();
        assert!(result.graph_changed);
        assert_eq!(result.split_segments, ["R"]);
        assert_eq!(result.supported_paths, 2);
        assert!(fs::read_to_string(overlay)
            .unwrap()
            .starts_with(&fs::read_to_string(&p).unwrap()));
        let g = parse_gfa(&out).unwrap();
        assert_eq!(g.names.len(), 6);
        assert_eq!(g.links.len(), 4);
        let walks = |start: &str| {
            let a = g.names.iter().position(|n| n == start).unwrap();
            let mid = g.links.iter().find(|(x, _)| *x == (a, true)).unwrap().1;
            let b = g.links.iter().find(|(x, _)| *x == mid).unwrap().1;
            (
                crate::linearize::spell(&g, &[(a, true), mid, b], false).unwrap(),
                g.names[b.0].clone(),
            )
        };
        assert_eq!(walks("A"), (b"AACGGCTATTGC".to_vec(), "B".into()));
        assert_eq!(walks("C"), (b"CAGTGCTAATGC".to_vec(), "D".into()));
    }
    #[test]
    fn stale_or_fabricated_evidence_is_rejected_before_writing() {
        let (p, e) = fixture([80, 0, 0, 80], "resolved");
        let mut g = parse_gfa(&p).unwrap();
        g.seqs[0][0] = b'T';
        assert!(Constraints::load(&g, &e).is_err());
        let mut v: Value = serde_json::from_str(&fs::read_to_string(&e).unwrap()).unwrap();
        v["repeat_nodes"][0]["phases"][0]["supported"] = json!(false);
        fs::write(&e, v.to_string()).unwrap();
        assert!(Constraints::load(&parse_gfa(&p).unwrap(), &e).is_err());
    }

    #[test]
    fn primary_enforces_a_pairing_when_the_root_is_the_repeat() {
        let (p, e) = fixture([80, 0, 0, 80], "resolved");
        let out = p.with_file_name("primary.fa");
        let json = p.with_file_name("primary.json");
        let report = crate::primary::run_phased(&p, &out, &json, 10, 10000, 0.9, Some(&e)).unwrap();
        let g = parse_gfa(&p).unwrap();
        let c = Constraints::load(&g, &e).unwrap();
        assert!(report
            .records
            .iter()
            .any(|r| r.path.split_whitespace().count() == 3));
        for record in &report.records {
            let walk: Vec<_> = record
                .path
                .split_whitespace()
                .map(|s| {
                    let (n, sign) = s.split_at(s.len() - 1);
                    (g.names.iter().position(|v| v == n).unwrap(), sign == "+")
                })
                .collect();
            assert!(c.allows_walk(&walk, record.circular), "{}", record.path);
        }
    }
}
