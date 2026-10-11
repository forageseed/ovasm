//! Regression tests and development diagnostics for cover selection across molecule counts.
//!
//! These tests separate candidate discovery from representative selection. They compare the
//! representative result with the existing score on all circular covers of a small
//! graph. This is an audit of a conditional score, not a physical-chromosome inference model.

use super::*;
use serde_json::json;

fn shared_repeat_case(
    name: &str,
    self_reads: u64,
    cross_reads: u64,
) -> (Graph, Value, LinearizeReport) {
    let dir = std::env::temp_dir().join(format!("ovasm-multichrom-{name}-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let gfa = format!(
        "H\tVN:Z:1.0\nS\tA\t{}\nS\tB\t{}\nS\tR\t{}\n\
         L\tA\t+\tR\t+\t0M\nL\tB\t+\tR\t+\t0M\n\
         L\tR\t+\tA\t+\t0M\nL\tR\t+\tB\t+\t0M\n",
        "ACGT".repeat(2000),
        "AGCT".repeat(3000),
        "ATGC".repeat(250)
    );
    let mut paths = vec![
        json!({"path": "A+ R+ A+", "reads": self_reads}),
        json!({"path": "B+ R+ B+", "reads": self_reads}),
    ];
    if cross_reads > 0 {
        paths.extend([
            json!({"path": "A+ R+ B+", "reads": cross_reads}),
            json!({"path": "B+ R+ A+", "reads": cross_reads}),
        ]);
    }
    let ev = json!({
        "repeat_segments": ["R"], "phased_paths": paths, "links": [],
        "segment_depths": {"A": 40.0, "B": 40.0, "R": 80.0}, "median_depth": 40.0,
    });
    fs::write(dir.join("graph.gfa"), gfa).unwrap();
    fs::write(dir.join("evidence.json"), ev.to_string()).unwrap();
    let g = parse_gfa(&dir.join("graph.gfa")).unwrap();
    let result = run(
        &dir.join("graph.gfa"),
        &dir.join("evidence.json"),
        LinearizeParams {
            min_reads: 2,
            max_steps: 100_000,
            max_molecules: 4,
            max_alternatives: 10,
            optional_below: 0.5,
        },
        &dir.join("molecules.fasta"),
        None,
        &dir.join("linearization.json"),
    )
    .unwrap();
    (g, ev, result)
}

/// Exhaust the two-anchor circular covers using production bridge and search functions.
/// Both possible covers have two repeat copies and two closing/connecting bridges; neither
/// wins merely because it is shorter or carries fewer negative log-probability terms.
fn all_circular_covers(g: &Graph, ev: &Value) -> Vec<CoverReport> {
    let (br, _) = bridges(g, ev, &[false, false, true], 2).unwrap();
    let fractions = exit_fractions(&br);
    let mut seen = BTreeSet::new();
    let mut reports = Vec::new();
    for m in 1..=2 {
        let (covers, truncated) = search(3, &[0, 1], &[false; 3], &br, m, 100_000);
        assert!(!truncated);
        for cover in covers {
            if cover.iter().any(|f| !f.circular) {
                continue;
            }
            let mut key: Vec<_> = cover
                .iter()
                .map(|f| canonical_key(&candidate_walk(&br, f), true))
                .collect();
            key.sort();
            if !seen.insert(key) {
                continue;
            }
            let molecules: Vec<_> = cover
                .iter()
                .map(|f| {
                    let walk = candidate_walk(&br, f);
                    let fs: Vec<_> = f.bridges.iter().map(|&i| fractions[i]).collect();
                    Linearization {
                        path: walk
                            .iter()
                            .map(|&o| g.label(o))
                            .collect::<Vec<_>>()
                            .join(" "),
                        circular: true,
                        length: spell(g, &walk, true).unwrap().len(),
                        score: path_score(&fs),
                        bridge_reads: f.bridges.iter().map(|&i| br[i].reads).collect(),
                        exit_fractions: fs,
                    }
                })
                .collect();
            reports.push(CoverReport {
                score: molecules.iter().map(|m| m.score).sum(),
                molecules,
                weak_bridges: 0,
                topological_bridges: cover
                    .iter()
                    .flat_map(|f| &f.bridges)
                    .filter(|&&i| br[i].topological)
                    .count(),
                copy_deviation: copy_deviation(&[None, None, Some(2)], &br, &cover),
            });
        }
    }
    reports.sort_by(|a, b| cover_order(a, b).then(a.molecules.len().cmp(&b.molecules.len())));
    reports
}

fn record(name: &str, baseline: &LinearizeReport, covers: &[CoverReport]) {
    println!(
        "MULTICHROM_DIAGNOSTIC {}",
        json!({
            "case": name, "baseline_molecules": baseline.molecules,
            "baseline_decisive": baseline.decisive,
            "baseline_score": baseline.components[0].best.as_ref().unwrap().ranked_score(),
            "all_circular_covers": covers,
        })
    );
}

#[test]
fn stronger_two_circle_cover_wins_across_molecule_counts() {
    let (g, ev, baseline) = shared_repeat_case("split-100-10", 100, 10);
    let covers = all_circular_covers(&g, &ev);
    record("split-100-10", &baseline, &covers);
    assert_eq!(
        baseline.molecules, 2,
        "the stronger bridge evidence must beat a smaller cover"
    );
    assert_eq!(covers.len(), 2);
    assert!(baseline.decisive);
    assert!(baseline.margin_ln.unwrap() > 4.6);
    assert_eq!(covers[0].molecules.len(), 2);
    assert!(covers[0].ranked_score() - covers[1].ranked_score() > 4.6);
    assert!(covers.iter().all(|c| c.copy_deviation == 0));
    assert_eq!(
        covers
            .iter()
            .map(|c| c.molecules.iter().map(|m| m.length).sum::<usize>())
            .collect::<Vec<_>>(),
        vec![22000, 22000]
    );
}

#[test]
fn stronger_single_circle_cover_is_preserved_when_all_sizes_are_compared() {
    let (g, ev, baseline) = shared_repeat_case("joined-10-100", 10, 100);
    let covers = all_circular_covers(&g, &ev);
    record("joined-10-100", &baseline, &covers);
    assert_eq!(baseline.molecules, 1);
    assert_eq!(covers[0].molecules.len(), 1);
    assert!(covers[0].ranked_score() - covers[1].ranked_score() > 4.6);
}

#[test]
fn equal_support_across_cover_sizes_is_reported_as_ambiguous() {
    let (g, ev, baseline) = shared_repeat_case("balanced-100-100", 100, 100);
    let covers = all_circular_covers(&g, &ev);
    record("balanced-100-100", &baseline, &covers);
    assert_eq!(baseline.molecules, 1);
    assert!(
        !baseline.decisive,
        "equal bridge support across cover sizes is not decisive"
    );
    assert_eq!(covers.len(), 2);
    assert!((covers[0].ranked_score() - covers[1].ranked_score()).abs() < 1e-12);
    assert_ne!(covers[0].molecules.len(), covers[1].molecules.len());
    assert!(baseline.margin_ln.unwrap().abs() < 1e-12);
}

#[test]
fn a_small_read_score_lead_does_not_imply_confidence() {
    let (g, ev, result) = shared_repeat_case("weak-split-12-10", 12, 10);
    let covers = all_circular_covers(&g, &ev);
    assert_eq!(result.molecules, 2);
    assert!(!result.decisive);
    assert!(result.margin_ln.unwrap() > 0.0);
    assert!(result.margin_ln.unwrap() < DECISIVE_LN);
    assert_eq!(covers[0].molecules.len(), 2);
}

#[test]
fn unphased_topology_selects_a_representative_and_remains_ambiguous() {
    let (g, ev, result) = shared_repeat_case("unphased-0-0", 0, 0);
    let covers = all_circular_covers(&g, &ev);
    assert_eq!(result.molecules, 1);
    assert!(!result.decisive);
    assert_eq!(covers.len(), 2);
    assert!(covers.iter().all(|c| c.topological_bridges == 2));
    assert_eq!(result.components[0].alternatives[0].molecules.len(), 2);
}

#[test]
fn omitted_larger_cover_sizes_cannot_establish_confidence() {
    let (_, _, full) = shared_repeat_case("size-limited-100-10", 100, 10);
    let dir = Path::new(&full.graph).parent().unwrap();
    let limited = run(
        Path::new(&full.graph),
        Path::new(&full.evidence),
        LinearizeParams {
            min_reads: 2,
            max_steps: 100_000,
            max_molecules: 1,
            max_alternatives: 10,
            optional_below: 0.5,
        },
        &dir.join("size-limited.fasta"),
        None,
        &dir.join("size-limited.json"),
    )
    .unwrap();
    assert_eq!(limited.molecules, 1);
    assert!(!limited.search_truncated);
    assert!(limited.cover_size_limited);
    assert_eq!(limited.components[0].searched_max_molecules, 1);
    assert!(
        !limited.decisive,
        "the stronger two-circle cover was not searched"
    );
    assert!(!full.cover_size_limited);
    assert_eq!(full.molecules, 2);
    assert!(full.decisive);
}

#[test]
fn a_truncated_search_with_a_valid_cover_is_not_decisive() {
    let (_, _, full) = shared_repeat_case("truncated-100-10", 100, 10);
    let dir = Path::new(&full.graph).parent().unwrap();
    let limited = run(
        Path::new(&full.graph),
        Path::new(&full.evidence),
        LinearizeParams {
            min_reads: 2,
            max_steps: 8,
            max_molecules: 4,
            max_alternatives: 10,
            optional_below: 0.5,
        },
        &dir.join("limited.fasta"),
        None,
        &dir.join("limited.json"),
    )
    .unwrap();
    assert!(limited.molecules > 0, "there is a representative to report");
    assert!(limited.search_truncated);
    assert!(!limited.decisive);
}

#[test]
fn chromosomes_without_cross_bridges_are_already_separate() {
    let (_, _, baseline) = shared_repeat_case("separate-100-0", 100, 0);
    assert_eq!(baseline.molecules, 2);
    assert_eq!(baseline.components.len(), 2);
    assert!(baseline
        .components
        .iter()
        .all(|c| c.best.as_ref().unwrap().molecules[0].circular));
}

/// Read-only replay of existing graphs; never assemble or modify the original run artifacts.
/// The candidate is deliberately only an audit: extend the current cover order across sizes.
/// A truncated search is reported, and is not evidence that its best candidate is optimal.
#[test]
#[ignore = "requires OVASM_MULTICHROM_MANIFEST pointing to server artifacts"]
fn audit_existing_graphs_across_cover_sizes() {
    let manifest = std::env::var("OVASM_MULTICHROM_MANIFEST").unwrap();
    let cases: Value = serde_json::from_str(&fs::read_to_string(manifest).unwrap()).unwrap();
    for case in cases.as_array().unwrap() {
        let name = case["name"].as_str().unwrap();
        let out = std::path::PathBuf::from(case["out"].as_str().unwrap());
        fs::create_dir_all(&out).unwrap();
        let graph = Path::new(case["graph"].as_str().unwrap());
        let evidence_path = Path::new(case["evidence"].as_str().unwrap());
        let baseline = run(
            graph,
            evidence_path,
            LinearizeParams {
                min_reads: 2,
                max_steps: 5_000_000,
                max_molecules: 4,
                max_alternatives: 10,
                optional_below: 0.5,
            },
            &out.join("baseline.fasta"),
            None,
            &out.join("baseline.json"),
        )
        .unwrap();
        let g = parse_gfa(graph).unwrap();
        let ev: Value = serde_json::from_str(&fs::read_to_string(evidence_path).unwrap()).unwrap();
        let repeats: BTreeSet<_> = ev["repeat_segments"]
            .as_array()
            .unwrap()
            .iter()
            .filter_map(Value::as_str)
            .collect();
        let repeat: Vec<_> = g
            .names
            .iter()
            .map(|n| repeats.contains(n.as_str()))
            .collect();
        let optional_names: BTreeSet<_> = baseline.optional_anchors.iter().collect();
        let optional: Vec<_> = g.names.iter().map(|n| optional_names.contains(n)).collect();
        let (br, _) = bridges(&g, &ev, &repeat, 2).unwrap();
        let fractions = exit_fractions(&br);
        let expected: Vec<Option<usize>> = g
            .names
            .iter()
            .enumerate()
            .map(|(s, n)| {
                let d = ev["segment_depths"][n].as_f64()?;
                let md = ev["median_depth"].as_f64()?;
                (repeat[s] && md > 0.0).then(|| ((d / md).round() as usize).max(1))
            })
            .collect();
        let mut component_reports = Vec::new();
        let mut selected = Vec::new();
        let mut any_truncated = false;
        for comp in components(g.names.len(), &br) {
            if comp.iter().all(|&a| optional[a]) {
                continue;
            }
            let mut by_size = Vec::new();
            let mut all = Vec::new();
            for m in 1..=4 {
                let (covers, truncated) = search(g.names.len(), &comp, &optional, &br, m, 500_000);
                any_truncated |= truncated;
                let mut seen = BTreeSet::new();
                let mut ranked = Vec::new();
                for cover in covers {
                    let mut key: Vec<_> = cover
                        .iter()
                        .map(|f| {
                            (
                                f.circular,
                                canonical_key(&candidate_walk(&br, f), f.circular),
                            )
                        })
                        .collect();
                    key.sort();
                    if !seen.insert(key) {
                        continue;
                    }
                    let molecules: Vec<_> = cover
                        .iter()
                        .map(|f| {
                            let walk = candidate_walk(&br, f);
                            let fs: Vec<_> = f.bridges.iter().map(|&i| fractions[i]).collect();
                            Linearization {
                                path: walk
                                    .iter()
                                    .map(|&o| g.label(o))
                                    .collect::<Vec<_>>()
                                    .join(" "),
                                circular: f.circular,
                                length: spell(&g, &walk, f.circular).unwrap().len(),
                                score: path_score(&fs),
                                bridge_reads: f.bridges.iter().map(|&i| br[i].reads).collect(),
                                exit_fractions: fs,
                            }
                        })
                        .collect();
                    let report = CoverReport {
                        score: molecules.iter().map(|m| m.score).sum(),
                        molecules,
                        weak_bridges: cover
                            .iter()
                            .flat_map(|f| &f.bridges)
                            .filter(|&&i| br[i].weak)
                            .count(),
                        topological_bridges: cover
                            .iter()
                            .flat_map(|f| &f.bridges)
                            .filter(|&&i| br[i].topological)
                            .count(),
                        copy_deviation: copy_deviation(&expected, &br, &cover),
                    };
                    ranked.push((report, cover));
                }
                ranked.sort_by(|(a, _), (b, _)| cover_order(a, b));
                by_size.push(json!({"m": m, "candidates": ranked.len(), "truncated": truncated, "best": ranked.first().map(|(r, _)| r)}));
                all.extend(ranked.into_iter().take(1));
            }
            all.sort_by(|(a, _), (b, _)| {
                cover_order(a, b).then(a.molecules.len().cmp(&b.molecules.len()))
            });
            if let Some((_, cover)) = all.first() {
                let mut final_cover = Vec::new();
                for f in cover {
                    if let Some((circle, tail, _)) = cut_tail(&g, &br, f) {
                        final_cover.push(circle);
                        if candidate_walk(&br, &tail)
                            .iter()
                            .any(|o| !repeat[o.0] && !optional[o.0])
                        {
                            final_cover.push(tail);
                        }
                    } else {
                        final_cover.push(f.clone());
                    }
                }
                for f in final_cover {
                    let walk = candidate_walk(&br, &f);
                    let fs: Vec<_> = f.bridges.iter().map(|&i| fractions[i]).collect();
                    let m = Linearization {
                        path: walk
                            .iter()
                            .map(|&o| g.label(o))
                            .collect::<Vec<_>>()
                            .join(" "),
                        circular: f.circular,
                        length: spell(&g, &walk, f.circular).unwrap().len(),
                        score: path_score(&fs),
                        bridge_reads: f.bridges.iter().map(|&i| br[i].reads).collect(),
                        exit_fractions: fs,
                    };
                    selected.push((walk, m));
                }
            }
            component_reports.push(json!({
                "anchors": comp.iter().map(|&a| &g.names[a]).collect::<Vec<_>>(),
                "by_size": by_size, "naive_all_size_best": all.first().map(|(r, _)| r),
            }));
        }
        // Preserve production's unbridged standalone molecules. They do not participate in
        // cover-size selection and must not disappear from the experimental comparison.
        let bridged_names: BTreeSet<_> = br
            .iter()
            .flat_map(|b| [&g.names[b.walk[0].0], &g.names[b.walk.last().unwrap().0]])
            .collect();
        for comp in &baseline.components {
            if comp.anchors.iter().all(|n| !bridged_names.contains(n)) {
                if let Some(best) = &comp.best {
                    for m in &best.molecules {
                        let walk = m
                            .path
                            .split_whitespace()
                            .map(|s| {
                                (
                                    g.names.iter().position(|n| n == &s[..s.len() - 1]).unwrap(),
                                    s.ends_with('+'),
                                )
                            })
                            .collect();
                        selected.push((walk, m.clone()));
                    }
                }
            }
        }
        fs::write(
            out.join("all-size.fasta"),
            fasta_records(&g, &selected, "audit", false, "").unwrap(),
        )
        .unwrap();
        let report = json!({
            "name": name, "scope": "seen graph diagnostic; naive extension of existing order, not a product candidate",
            "budget_per_size": 500_000, "any_search_truncated": any_truncated,
            "baseline_molecules": baseline.molecules, "baseline_decisive": baseline.decisive,
            "baseline_search_truncated": baseline.search_truncated,
            "all_size_molecules": selected.len(), "components": component_reports,
        });
        fs::write(
            out.join("audit.json"),
            serde_json::to_string_pretty(&report).unwrap(),
        )
        .unwrap();
        println!(
            "REAL_MULTICHROM_DIAGNOSTIC {name}: baseline={} all_sizes={} truncated={}",
            baseline.molecules,
            selected.len(),
            any_truncated
        );
    }
}
