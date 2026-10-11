//! Primary sequence representation: a path rooted in the largest depth-weighted SCC,
//! plus segments that path leaves unplaced. Unlike an anchor cover this does not require
//! joining every unique segment, and it preserves unplaced repeats as sequence records.
//! A primary representation does not establish a unique set of physical molecules.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

use anyhow::{bail, Context, Result};
use serde::Serialize;

use crate::evidence::{flip, parse_gfa, Graph, Oriented};
use crate::linearize::spell;

#[derive(Clone, Debug)]
struct Candidate {
    walk: Vec<Oriented>,
    circular: bool,
}

#[derive(Serialize)]
pub struct PrimaryRecord {
    pub path: String,
    pub circular: bool,
    pub length: usize,
    pub kind: &'static str,
    /// Assembly representation, never an assertion of chromosome completeness.
    pub complete: bool,
}

#[derive(Debug, Serialize)]
pub struct PrimaryComponent {
    pub root: String,
    pub depth: f64,
    pub copy_budget: BTreeMap<String, usize>,
    pub candidates: usize,
    pub expansions: usize,
    pub selected_path: String,
    pub unplaced_segments: Vec<String>,
}

#[derive(Serialize)]
pub struct PrimaryReport {
    pub representation: &'static str,
    pub decisive: bool,
    pub total_length: usize,
    pub unplaced_segments: usize,
    pub components: Vec<PrimaryComponent>,
    pub records: Vec<PrimaryRecord>,
}

fn slot(o: Oriented) -> usize {
    2 * o.0 + usize::from(o.1)
}
fn oriented(i: usize) -> Oriented {
    (i / 2, i % 2 == 1)
}

/// Kosaraju with explicit stacks: a large acyclic graph does not recurse on the call stack.
fn scc(adj: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut seen = vec![false; adj.len()];
    let mut order = Vec::new();
    for root in 0..adj.len() {
        if seen[root] {
            continue;
        }
        seen[root] = true;
        let mut stack = vec![(root, 0)];
        while let Some((v, next)) = stack.last_mut() {
            if *next == adj[*v].len() {
                order.push(*v);
                stack.pop();
                continue;
            }
            let child = adj[*v][*next];
            *next += 1;
            if !seen[child] {
                seen[child] = true;
                stack.push((child, 0));
            }
        }
    }
    let mut reverse = vec![Vec::new(); adj.len()];
    for (v, next) in adj.iter().enumerate() {
        for &u in next {
            reverse[u].push(v);
        }
    }
    seen.fill(false);
    let mut blocks = Vec::new();
    for root in order.into_iter().rev() {
        if seen[root] {
            continue;
        }
        let mut block = Vec::new();
        let mut stack = vec![root];
        seen[root] = true;
        while let Some(v) = stack.pop() {
            block.push(v);
            for &u in &reverse[v] {
                if !seen[u] {
                    seen[u] = true;
                    stack.push(u);
                }
            }
        }
        block.sort_unstable();
        blocks.push(block);
    }
    blocks
}

fn extend(
    initial: Vec<Oriented>,
    budget: &[usize],
    adj: &[Vec<usize>],
    expansions: &mut usize,
    max_steps: usize,
) -> Result<Vec<Candidate>> {
    let mut counts = vec![0; budget.len()];
    for o in &initial {
        counts[o.0] += 1;
    }
    let mut stack = vec![(initial, counts)];
    let mut terminal = Vec::new();
    while let Some((walk, counts)) = stack.pop() {
        *expansions += 1;
        if *expansions > max_steps {
            bail!("primary path search exceeded {max_steps} expansions");
        }
        let next = &adj[slot(*walk.last().unwrap())];
        let circular = next.contains(&slot(walk[0]));
        let children: Vec<Oriented> = next
            .iter()
            .map(|&i| oriented(i))
            .filter(|&o| o != walk[0] && counts[o.0] < budget[o.0])
            .collect();
        if circular || children.is_empty() {
            terminal.push(Candidate {
                walk: walk.clone(),
                circular,
            });
        }
        for child in children {
            let mut extended = walk.clone();
            extended.push(child);
            let mut copied = counts.clone();
            copied[child.0] += 1;
            stack.push((extended, copied));
        }
    }
    Ok(terminal)
}

fn solve(
    g: &Graph,
    depth: &[f64],
    max_copy: usize,
    max_steps: usize,
    circular_fraction: f64,
) -> Result<(Vec<(Candidate, &'static str)>, Vec<PrimaryComponent>)> {
    let n = g.names.len();
    let mut adj = vec![Vec::new(); 2 * n];
    for &(a, b) in &g.links {
        adj[slot(a)].push(slot(b));
        adj[slot(flip(b))].push(slot(flip(a)));
    }
    for neighbors in &mut adj {
        neighbors.sort_unstable();
        neighbors.dedup();
    }
    let blocks = scc(&adj);
    let mut seen = vec![false; n];
    let mut records = Vec::new();
    let mut reports = Vec::new();
    for start in 0..n {
        if seen[start] {
            continue;
        }
        let mut stack = vec![start];
        let mut comp = BTreeSet::new();
        seen[start] = true;
        while let Some(v) = stack.pop() {
            comp.insert(v);
            for i in [2 * v, 2 * v + 1] {
                for &u in &adj[i] {
                    if !seen[u / 2] {
                        seen[u / 2] = true;
                        stack.push(u / 2);
                    }
                }
            }
        }
        let weight = |v: usize| g.seqs[v].len() as f64 * depth[v];
        let block_weight = |block: &&Vec<usize>| {
            block
                .iter()
                .map(|i| i / 2)
                .collect::<BTreeSet<_>>()
                .iter()
                .map(|&v| weight(v))
                .sum::<f64>()
        };
        let block = blocks
            .iter()
            .filter(|b| comp.contains(&(b[0] / 2)))
            .max_by(|a, b| {
                block_weight(a)
                    .total_cmp(&block_weight(b))
                    .then_with(|| b.cmp(a))
            })
            .unwrap();
        let root = *block
            .iter()
            .max_by(|&&a, &&b| {
                weight(a / 2)
                    .total_cmp(&weight(b / 2))
                    .then_with(|| b.cmp(&a))
            })
            .unwrap();
        let mut by_depth: Vec<usize> = comp.iter().copied().collect();
        by_depth.sort_by(|&a, &b| depth[a].total_cmp(&depth[b]).then(a.cmp(&b)));
        let total: usize = comp.iter().map(|&v| g.seqs[v].len()).sum();
        let mut accumulated = 0;
        let median = depth[*by_depth
            .iter()
            .find(|&&v| {
                accumulated += g.seqs[v].len();
                2 * accumulated >= total
            })
            .unwrap()];
        let mut budget = vec![0; n];
        for &v in &comp {
            budget[v] = ((depth[v] / median).round() as usize).clamp(1, max_copy);
        }
        let mut expansions = 0;
        let mut candidates = Vec::new();
        for first in extend(
            vec![oriented(root)],
            &budget,
            &adj,
            &mut expansions,
            max_steps,
        )? {
            if first.circular {
                candidates.push(first);
            } else {
                let reverse = first.walk.into_iter().rev().map(flip).collect();
                candidates.extend(extend(reverse, &budget, &adj, &mut expansions, max_steps)?);
            }
        }
        let mut ranked = Vec::new();
        for candidate in candidates {
            let length = spell(g, &candidate.walk, candidate.circular)?.len();
            let weighted: f64 = candidate.walk.iter().map(|o| weight(o.0)).sum();
            ranked.push((candidate, length, weighted));
        }
        ranked.sort_by(|a, b| {
            b.1.cmp(&a.1)
                .then(b.2.total_cmp(&a.2))
                .then(a.0.walk.cmp(&b.0.walk))
        });
        let mut chosen = 0;
        if let Some(i) = ranked.iter().position(|c| c.0.circular) {
            if ranked[i].1 as f64 >= circular_fraction * ranked[0].1 as f64 {
                chosen = i;
            }
        }
        let chosen = ranked[chosen].0.clone();
        let placed: BTreeSet<usize> = chosen.walk.iter().map(|o| o.0).collect();
        let remaining: Vec<usize> = comp.difference(&placed).copied().collect();
        reports.push(PrimaryComponent {
            root: g.label(oriented(root)),
            depth: median,
            copy_budget: comp
                .iter()
                .map(|&v| (g.names[v].clone(), budget[v]))
                .collect(),
            candidates: ranked.len(),
            expansions,
            selected_path: chosen
                .walk
                .iter()
                .map(|&o| g.label(o))
                .collect::<Vec<_>>()
                .join(" "),
            unplaced_segments: remaining.iter().map(|&v| g.names[v].clone()).collect(),
        });
        records.push((chosen, "primary_path"));
        records.extend(remaining.into_iter().map(|v| {
            (
                Candidate {
                    walk: vec![(v, true)],
                    circular: false,
                },
                "unplaced_segment",
            )
        }));
    }
    Ok((records, reports))
}

pub fn run(
    gfa: &Path,
    out: &Path,
    out_json: &Path,
    max_copy: usize,
    max_steps: usize,
    circular_fraction: f64,
) -> Result<PrimaryReport> {
    if max_copy == 0 || max_steps == 0 || !(0.0..=1.0).contains(&circular_fraction) {
        bail!("invalid primary path search parameters");
    }
    let g = parse_gfa(gfa)?;
    let mut tagged = BTreeMap::new();
    for line in fs::read_to_string(gfa)?.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f[0] != "S" {
            continue;
        }
        let tag = |prefix: &str| {
            f.iter()
                .skip(3)
                .find_map(|t| t.strip_prefix(prefix).and_then(|s| s.parse::<f64>().ok()))
        };
        let d = tag("dm:f:")
            .or_else(|| tag("dp:f:"))
            .or_else(|| tag("KC:i:").map(|c| c / f[2].len() as f64))
            .with_context(|| format!("segment {} has no depth for primary copy inference", f[1]))?;
        if !d.is_finite() || d <= 0.0 {
            bail!("segment {} has invalid primary depth", f[1]);
        }
        tagged.insert(f[1].to_owned(), d);
    }
    let depth: Vec<f64> = g.names.iter().map(|n| tagged[n]).collect();
    let (paths, components) = solve(&g, &depth, max_copy, max_steps, circular_fraction)?;
    let mut records = Vec::new();
    let mut fasta = String::new();
    for (i, (candidate, kind)) in paths.iter().enumerate() {
        let seq = spell(&g, &candidate.walk, candidate.circular)?;
        let path = candidate
            .walk
            .iter()
            .map(|&o| g.label(o))
            .collect::<Vec<_>>()
            .join(" ");
        fasta.push_str(&format!(
            ">primary.{} circular={} length={} kind={} complete=false path={}\n",
            i + 1,
            candidate.circular,
            seq.len(),
            kind,
            path.replace(' ', ",")
        ));
        for chunk in seq.chunks(80) {
            fasta.push_str(std::str::from_utf8(chunk)?);
            fasta.push('\n');
        }
        records.push(PrimaryRecord {
            path,
            circular: candidate.circular,
            length: seq.len(),
            kind,
            complete: false,
        });
    }
    let report = PrimaryReport {
        representation: "primary_path_and_unplaced_segments",
        decisive: false,
        total_length: records.iter().map(|r| r.length).sum(),
        unplaced_segments: records
            .iter()
            .filter(|r| r.kind == "unplaced_segment")
            .count(),
        records,
        components,
    };
    fs::create_dir_all(out.parent().unwrap_or(Path::new(".")))?;
    fs::write(out, fasta)?;
    fs::write(out_json, serde_json::to_string_pretty(&report)?)?;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(lengths: &[usize], links: &[(Oriented, Oriented)]) -> Graph {
        Graph {
            names: (0..lengths.len()).map(|i| format!("u{i}")).collect(),
            seqs: lengths.iter().map(|&n| vec![b'A'; n]).collect(),
            links: links.to_vec(),
            overlaps: vec![0; links.len()],
        }
    }

    #[test]
    fn eight_independent_circles_are_eight_primary_paths() {
        let g = graph(
            &[100; 8],
            &(0..8).map(|i| ((i, true), (i, true))).collect::<Vec<_>>(),
        );
        let (records, _) = solve(&g, &[10.0; 8], 10, 1000, 0.9).unwrap();
        assert_eq!(records.len(), 8);
        assert!(records
            .iter()
            .all(|r| r.0.circular && r.1 == "primary_path"));
    }

    #[test]
    fn root_comes_from_the_largest_strongly_connected_block() {
        let g = graph(
            &[60, 60, 100],
            &[
                ((0, true), (1, true)),
                ((1, true), (0, true)),
                ((1, true), (2, true)),
            ],
        );
        let (_, report) = solve(&g, &[10.0; 3], 10, 1000, 0.9).unwrap();
        assert!(!report[0].root.starts_with("u2"));
    }

    #[test]
    fn truncated_primary_search_does_not_fabricate_a_cover() {
        let g = graph(&[100, 100], &[((0, true), (1, true))]);
        assert!(solve(&g, &[10.0; 2], 10, 1, 0.9)
            .unwrap_err()
            .to_string()
            .contains("exceeded"));
    }
}
