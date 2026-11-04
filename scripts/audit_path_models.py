#!/usr/bin/env python3
"""Audit existing linearize candidates using single-record graph walks.

Research diagnostic only: never rewrites a graph/FASTA or certifies topology.
Counts are record-walk occurrences, NOT independent reads or molecule abundance.
"""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def flip(label):
    if not label or label[-1] not in '+-':
        raise ValueError(f'invalid oriented segment: {label!r}')
    return label[:-1] + ('-' if label[-1] == '+' else '+')


def reverse(walk):
    return tuple(flip(x) for x in reversed(walk))


def canonical(walk):
    walk = tuple(walk)
    return min(walk, reverse(walk))


def graph(path):
    nodes, edges = {}, set()
    for line in Path(path).read_text().splitlines():
        f = line.split('\t')
        if f[0] == 'S':
            if len(f) < 3 or f[1] in nodes or f[2] == '*':
                raise ValueError('duplicate segment or missing sequence')
            nodes[f[1]] = f[2]
        elif f[0] == 'L':
            if len(f) < 6:
                raise ValueError('malformed GFA link')
            a, b = f[1] + f[2], f[3] + f[4]
            flip(a); flip(b)
            edges.update(((a, b), (flip(b), flip(a))))
    if not nodes:
        raise ValueError('empty graph')
    predecessors, successors = {}, {}
    for a, b in edges:
        if a[:-1] not in nodes or b[:-1] not in nodes:
            raise ValueError('link references unknown segment')
        predecessors.setdefault(b, set()).add(a)
        successors.setdefault(a, set()).add(b)
    repeats = {name for name in nodes
               if len(predecessors.get(name + '+', ())) > 1
               or len(successors.get(name + '+', ())) > 1}
    return nodes, edges, repeats


def validate_walk(walk, nodes, edges, circular=False):
    if not walk:
        raise ValueError('empty path')
    for x in walk:
        flip(x)
        if x[:-1] not in nodes:
            raise ValueError(f'unknown segment {x}')
    links = list(zip(walk, walk[1:]))
    if circular:
        links.append((walk[-1], walk[0]))
    if any(e not in edges for e in links):
        raise ValueError(f'path uses absent link: {walk}')


def contains(path, circular, observed):
    """Exact node-order linkage; includes both strands and genuine circular wrapping."""
    if not path or not observed:
        return False
    for query in (tuple(observed), reverse(observed)):
        if not circular and len(query) > len(path):
            continue
        stop = len(path) if circular else len(path) - len(query) + 1
        for i in range(stop):
            if all(path[(i + j) % len(path)] == x for j, x in enumerate(query)):
                return True
    return False


def observations(rows, nodes, edges):
    result = {}
    for row in rows:
        walk = tuple(row['path'].split())
        validate_walk(walk, nodes, edges)
        count = row['reads']
        if isinstance(count, bool) or not isinstance(count, int) or count < 0:
            raise ValueError('invalid observation count')
        key = canonical(walk)
        if key in result:
            raise ValueError('duplicate/reverse-complement observation; refusing to double-count')
        if count:
            result[key] = count
    return result


def audit(gfa, evidence, linearization, max_comparisons=2_000_000):
    nodes, edges, repeats = graph(gfa)
    ev = json.loads(Path(evidence).read_text())
    lin = json.loads(Path(linearization).read_text())
    # These legacy reports have no sequence fingerprint. Require their recorded graph
    # to still exist and match byte-for-byte, rather than trusting labels alone.
    root = Path(__file__).resolve().parents[1]
    bindings = {}
    for label, document in [('evidence', ev), ('linearization', lin)]:
        source = Path(document['graph'])
        if not source.is_absolute():
            source = root / source
        if not source.is_file() or sha(source) != sha(gfa):
            raise ValueError(f'{label} recorded graph missing or differs: {source}')
        bindings[label] = str(source.resolve())
    ev_edges = set()
    for row in ev['links']:
        a, b = row['from'], row['to']
        ev_edges.update(((a, b), (flip(b), flip(a))))
    if ev_edges != edges:
        raise ValueError('evidence links do not match the supplied graph')
    local = observations(ev.get('traversals', []), nodes, edges)
    walks = observations(ev.get('observed_walks', []), nodes, edges)
    joint = {w: n for w, n in walks.items() if len({x[:-1] for x in w} & repeats) >= 2}
    comparisons = 0
    components = []
    for number, comp in enumerate(lin['components'], 1):
        candidates = list(comp.get('alternatives', []))
        if comp.get('best'):
            candidates.insert(0, comp['best'])
        listed = []
        for ci, cover in enumerate(candidates):
            molecules = []
            for m in cover['molecules']:
                w = tuple(m['path'].split())
                circular = m['circular']
                if not isinstance(circular, bool):
                    raise ValueError('circular must be a boolean')
                validate_walk(w, nodes, edges, circular)
                molecules.append((w, circular))
            counts, types = {}, {}
            for label, obs in [('local', local), ('walk', walks), ('joint', joint)]:
                matched = []
                for w, n in obs.items():
                    comparisons += sum(len(p) * len(w) for p, _ in molecules)
                    if comparisons > max_comparisons:
                        raise ValueError('audit comparison budget exceeded; no partial ranking emitted')
                    if any(contains(p, c, w) for p, c in molecules):
                        matched.append((w, n))
                counts[label] = sum(n for _, n in matched)
                types[label] = len(matched)
            weak = cover.get('weak_bridges', 0)
            topo = cover.get('topological_bridges', 0)
            copydev = cover.get('copy_deviation', 0)
            # Explicit exploratory ordering, not a probability or structural truth.
            # Neither circularity nor number/length of molecules supplies a bonus.
            key = (counts['joint'], counts['walk'], counts['local'], -(weak + topo), -copydev)
            listed.append(dict(candidate=ci, original_selected=ci == 0,
                               molecules=[dict(path=' '.join(p), circular=c) for p, c in molecules],
                               matched_record_walk_occurrences=counts, matched_walk_types=types,
                               weak_bridges=weak, topological_bridges=topo,
                               copy_deviation=copydev, exploratory_key=key))
        order = sorted(listed, key=lambda c: tuple(-x for x in c['exploratory_key']))
        top = [c['candidate'] for c in order if c['exploratory_key'] == order[0]['exploratory_key']] if order else []
        top_signal = bool(order and any(order[0]['matched_record_walk_occurrences'].values()))
        components.append(dict(component=number, candidates_reported=comp.get('candidates', 0),
                               candidates_audited=len(listed), candidate_list_limited=comp.get('candidates', 0) > len(listed),
                               search_truncated=bool(comp.get('search_truncated', False)),
                               cover_size_limited=bool(comp.get('cover_size_limited', False)),
                               tails_cut=comp.get('tails_cut', []),
                               top_tied_candidates=top, has_matched_sequence_evidence=top_signal,
                               original_excluded_from_top=bool(top_signal and top and 0 not in top),
                               ordered_candidates=order))
    return dict(schema='ovasm.path-model-audit.v1',
                scope='Exploratory re-ranking of reported candidates only; not a new assembly or physical topology call.',
                graph=str(Path(gfa).resolve()), graph_sha256=sha(gfa),
                evidence_sha256=sha(evidence), linearization_sha256=sha(linearization),
                script=str(Path(__file__).resolve()), script_sha256=sha(__file__), graph_bindings=bindings,
                evidence_params=ev.get('params', {}),
                count_scope='Counts are record-walk occurrences. Overlapping chains can share records; legacy evidence does not certify globally unique IDs, absence of chimeras, or organelle origin.',
                score_scope='Lexicographic joint-walk, whole-walk, local-traversal support; then fewer weak/topological bridges and copy deviations. No topology bonus; no calibrated likelihood or abundance.',
                missing_walk_evidence='observed_walks' not in ev,
                local_walk_types=len(local), observed_walk_types=len(walks), joint_walk_types=len(joint),
                topology_branch_nodes=sorted(repeats), comparison_work=comparisons,
                physical_topology_identified=False, components=components)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('gfa', 'evidence', 'linearization', 'out'):
        p.add_argument('--' + name, required=True, type=Path)
    p.add_argument('--max-comparisons', type=int, default=2_000_000)
    a = p.parse_args()
    paths = [a.gfa.resolve(), a.evidence.resolve(), a.linearization.resolve()]
    if a.out.resolve() in paths:
        p.error('output must not overwrite an input')
    if a.max_comparisons < 1:
        p.error('--max-comparisons must be positive')
    try:
        result = audit(a.gfa, a.evidence, a.linearization, a.max_comparisons)
    except (ValueError, KeyError, OSError, TypeError) as e:
        p.exit(2, f'audit failed: {e}\n')
    a.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ('local_walk_types', 'observed_walk_types', 'joint_walk_types')}))
    print(json.dumps([dict(component=c['component'], candidates=c['candidates_audited'],
                           top_ties=len(c['top_tied_candidates']), changed=c['original_excluded_from_top'])
                      for c in result['components']]))


if __name__ == '__main__':
    main()
