"""Replay reachability through graph-bound, externally adjudicated frontier edges."""

from argparse import ArgumentParser
from collections import defaultdict
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile

MANIFEST = "ovasm.wp4.frontier-closure.v1"
GRAPH = "ovasm.wp4.frontier-graph.v1"
SCOPE = "ovasm.wp4.frontier-scope.v1"
EDGE = "ovasm.wp4.frontier-edge.v1"
REPORT = "ovasm.wp4.frontier-closure-report.v1"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def canonical_sha(value):
    return sha(json.dumps(value, sort_keys=True, separators=(",", ":")).encode())


def object_value(value, label):
    if not isinstance(value, dict):
        raise ValueError(f"{label}: expected an object")
    return value


def json_object(data, label):
    def unique_pairs(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"{label}: duplicate JSON key {key}")
            result[key] = value
        return result

    return object_value(json.loads(data, object_pairs_hook=unique_pairs), label)


def name(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label}: expected a nonempty identifier")
    return value


def names(values, label):
    if not isinstance(values, list):
        raise ValueError(f"{label}: expected a list of identifiers")
    return sorted({name(value, label) for value in values})


def lock(ref, directory, label, *, read_data=True):
    ref = object_value(ref, label)
    path = (directory / name(ref.get("path"), label)).resolve()
    digest = ref.get("sha256")
    if not isinstance(digest, str) or not re.fullmatch(r"[a-f0-9]{64}", digest):
        raise ValueError(f"{label}: expected a SHA-256 digest")
    if read_data:
        data = path.read_bytes()
        actual = sha(data)
    else:
        data = None
        with path.open("rb") as stream:
            actual = hashlib.file_digest(stream, "sha256").hexdigest()
    if actual != digest:
        raise ValueError(f"{label}: SHA-256 mismatch")
    return path, data, {"path": str(path), "sha256": digest}


def reference_identity(ref, directory, label):
    ref = object_value(ref, label)
    return ((directory / name(ref.get("path"), label)).resolve(),
            name(ref.get("sha256"), label))


def provenance(refs, directory, label, cache):
    if not isinstance(refs, list) or not refs:
        raise ValueError(f"{label}: provenance references are required")
    receipts = {}
    for ref in refs:
        key = reference_identity(ref, directory, label)
        if key not in cache:
            _, _, cache[key] = lock(ref, directory, label, read_data=False)
        info = cache[key]
        receipts[info["path"], info["sha256"]] = info
    return [receipts[key] for key in sorted(receipts)]


@dataclass(frozen=True)
class Edge:
    id: str
    source: str
    destination: str


def read_graph(raw):
    graph = json_object(raw, "graph")
    if graph.get("schema") != GRAPH:
        raise ValueError("unknown graph schema")
    nodes = set(names(graph.get("nodes"), "graph nodes"))
    if not nodes or not isinstance(graph.get("edges"), list):
        raise ValueError("graph requires nodes and an edge list")
    edges = {}
    for row in graph["edges"]:
        row = object_value(row, "graph edge")
        edge = Edge(name(row.get("id"), "edge id"), name(row.get("from"), "edge origin"),
                    name(row.get("to"), "edge destination"))
        if edge.source not in nodes or edge.destination not in nodes:
            raise ValueError(f"edge {edge.id}: endpoint not in graph")
        if edge.id in edges:
            raise ValueError(f"duplicate graph edge id {edge.id}")
        edges[edge.id] = edge
    return nodes, edges


def bind(receipt, expected, label):
    if any(receipt.get(key) != value for key, value in expected.items()):
        raise ValueError(f"{label}: scope, graph or fixed eligibility binding mismatch")


def read_evidence(ref, directory, binding, edges, provenance_cache):
    path, raw, info = lock(ref, directory, "edge evidence")
    record = json_object(raw, "edge evidence")
    bind(record, {**binding, "schema": EDGE}, "edge evidence")
    edge_id = name(record.get("edge_id"), "evidence edge id")
    if edge_id not in edges:
        raise ValueError(f"evidence edge {edge_id}: edge not in graph")

    decision = record.get("decision")
    basis = record.get("basis")
    if decision not in ("admitted", "excluded", "unresolved"):
        raise ValueError(f"edge {edge_id}: invalid evidence decision")
    basis_options = {
        "admitted": ("target_exact_bridge",),
        "excluded": ("competing_context", "calibrated_absence"),
        "unresolved": ("opportunity_only", "no_opportunity", "unmeasured", "ambiguous"),
    }
    if basis not in basis_options[decision]:
        raise ValueError(f"edge {edge_id}: invalid basis for {decision}")

    opportunities = names(record.get("opportunity_fragment_ids"), "opportunity fragment ids")
    support = names(record.get("support_fragment_ids"), "support fragment ids")
    if not set(support) <= set(opportunities):
        raise ValueError(f"edge {edge_id}: support is outside the opportunity set")
    if basis in ("no_opportunity", "unmeasured") and opportunities:
        raise ValueError(f"edge {edge_id}: unexpected opportunity IDs for {basis}")
    if basis == "opportunity_only" and not opportunities:
        raise ValueError(f"edge {edge_id}: opportunity_only requires identifiable fragments")
    if decision in ("admitted", "excluded"):
        if record.get("full_edge") is not True or not opportunities:
            raise ValueError(f"edge {edge_id}: decisive certificate needs a full edge and opportunities")
        if basis != "calibrated_absence" and not support:
            raise ValueError(f"edge {edge_id}: exact bridge certificate needs supporting fragments")

    detail = {
        "certificate": info, "decision": decision, "basis": basis,
        "_bridge_sequence": record.get("bridge_sequence"),
        "explanation": name(record.get("explanation"), "evidence explanation"),
        "opportunity_fragment_ids": opportunities, "support_fragment_ids": support,
        "opportunity_count": len(opportunities), "support_count": len(support),
        "provenance": provenance(record.get("provenance"), path.parent, f"edge {edge_id}", provenance_cache),
    }
    return edge_id, detail


def inspect_edge(records):
    if not records:
        return "unresolved", "no_edge_evidence"
    outcomes = {row["decision"] for row in records}
    if len(outcomes) > 1:
        return "unresolved", "conflicting_edge_evidence"
    decision, = outcomes
    return decision, "held_by_certificate" if decision == "unresolved" else "certified"


def closure(seeds, edges, ledger, scope_complete, round_budget=None, edge_budget=None):
    reachable = set(seeds)
    outgoing = defaultdict(list)
    for edge in edges.values():
        outgoing[edge.source].append(edge.id)

    evaluated = {}
    rounds = []
    pending = {e for seed in seeds for e in outgoing[seed]}
    while pending:
        stopped_by = []
        if round_budget is not None and len(rounds) == round_budget:
            stopped_by.append("round_budget")
        # An entire wave is the minimum evaluation unit. Partial waves would let edge
        # traversal order influence the result under the same budget.
        if edge_budget is not None and len(evaluated) + len(pending) > edge_budget:
            stopped_by.append("edge_budget")
        if stopped_by:
            break
        current = sorted(pending)
        newly_reachable = set()
        for edge_id in current:
            state, reason = inspect_edge(ledger.get(edge_id, []))
            evaluated[edge_id] = {"round": len(rounds) + 1, "state": state, "reason": reason}
            if state == "admitted":
                newly_reachable.add(edges[edge_id].destination)
        newly_reachable -= reachable
        reachable.update(newly_reachable)
        rounds.append({"round": len(rounds) + 1, "evaluated_edges": current,
                       "new_nodes": sorted(newly_reachable)})
        pending = {e for node in newly_reachable for e in outgoing[node] if e not in evaluated}
    else:
        stopped_by = []

    held = sorted(edge for edge, item in evaluated.items() if item["state"] == "unresolved")
    if stopped_by:
        state = "truncated"
    elif held or not scope_complete:
        state = "unresolved"
    else:
        state = "converged"

    decisions = {
        "seed_nodes": sorted(seeds), "reachable_nodes": sorted(reachable),
        "admitted_edges": sorted(edge for edge, item in evaluated.items() if item["state"] == "admitted"),
        "excluded_edges": sorted(edge for edge, item in evaluated.items() if item["state"] == "excluded"),
        "unresolved_edges": held, "unevaluated_frontier_edges": sorted(pending),
        "scope_complete": scope_complete, "frontier_exhausted": not pending,
        "state": state, "budget_reasons": stopped_by, "rounds": rounds,
    }
    return decisions, evaluated


def bridge_record(receipt, cache, edge_id):
    label = f"edge {edge_id}: bridge sequence"
    if receipt["_bridge_sequence"] is None:
        raise ValueError(f"{label}: admitted certificate is missing bridge_sequence")
    spec = object_value(receipt["_bridge_sequence"], label)
    record_id = name(spec.get("record"), f"{label} FASTA record")
    sequence_sha = spec.get("sequence_sha256")
    if not isinstance(sequence_sha, str) or not re.fullmatch(r"[a-f0-9]{64}", sequence_sha):
        raise ValueError(f"{label}: expected a normalized sequence SHA-256 digest")
    directory = Path(receipt["certificate"]["path"]).parent
    identity = reference_identity(spec, directory, label)
    if identity not in cache:
        _, raw, info = lock(spec, directory, label)
        cache[identity] = raw, info
    raw, info = cache[identity]
    token = record_id.encode("utf-8")
    selected = False
    pieces = []
    matches = 0
    for line in raw.splitlines():
        if line.startswith(b">"):
            fields = line[1:].split(None, 1)
            selected = bool(fields) and fields[0] == token
            matches += selected
        elif selected:
            pieces.append(b"".join(line.split()))
    if matches != 1:
        raise ValueError(f"{label}: FASTA record {record_id!r} must occur exactly once")
    sequence = b"".join(pieces).upper()
    if not re.fullmatch(rb"[ACGT]+", sequence):
        raise ValueError(f"{label}: exact bridge must contain only A, C, G and T")
    if sha(sequence) != sequence_sha:
        raise ValueError(f"{label}: normalized sequence SHA-256 mismatch")
    return sequence, {"edge_id": edge_id, "certificate": receipt["certificate"],
                      "fasta": info, "record": record_id}


def seed_export(decisions, records):
    cache = {}
    by_sequence = {}
    for edge_id in decisions["admitted_edges"]:
        sequence = None
        for receipt in records[edge_id]:
            candidate, source = bridge_record(receipt, cache, edge_id)
            if sequence is not None and candidate != sequence:
                raise ValueError(f"edge {edge_id}: conflicting bridge sequences")
            sequence = candidate
            by_sequence.setdefault(sequence, []).append(source)
    fasta = bytearray()
    metadata = []
    for sequence, sources in sorted(by_sequence.items(), key=lambda item: sha(item[0])):
        digest = sha(sequence)
        identifier = f"bridge_{digest}"
        fasta.extend(f">{identifier}\n".encode("ascii"))
        for start in range(0, len(sequence), 80):
            fasta.extend(sequence[start:start + 80] + b"\n")
        metadata.append({"record": identifier, "bp": len(sequence), "sequence_sha256": digest,
                         "admitted_edges": sorted({source["edge_id"] for source in sources}),
                         "source_records": sources})
    return bytes(fasta), {"record_count": len(metadata), "fasta_sha256": sha(fasta),
                         "budget_truncated": decisions["state"] == "truncated",
                         "records": metadata}


def replay(manifest_path, round_budget=None, edge_budget=None, *, export_seeds=False):
    if any(budget is not None and (type(budget) is not int or budget < 0)
           for budget in (round_budget, edge_budget)):
        raise ValueError("budgets must be nonnegative integers")
    raw = manifest_path.read_bytes()
    m = json_object(raw, "manifest")
    if m.get("schema") != MANIFEST:
        raise ValueError("unknown manifest schema")
    scope_id = name(m.get("scope"), "scope")
    base = manifest_path.parent
    _, graph_raw, graph = lock(m.get("graph"), base, "frozen graph")
    nodes, edges = read_graph(graph_raw)
    seeds = names(m.get("seeds"), "original seeds")
    if not seeds or not set(seeds) <= nodes:
        raise ValueError("original seeds must be present in graph")
    _, _, eligibility = lock(m.get("eligibility"), base, "fixed eligibility", read_data=False)
    _, _, read_pool = lock(m.get("read_pool"), base, "frozen read-pool receipt", read_data=False)
    binding = {"scope": scope_id, "graph_sha256": graph["sha256"],
               "eligibility_sha256": eligibility["sha256"], "read_pool_sha256": read_pool["sha256"],
               "seed_set_sha256": canonical_sha(seeds)}
    scope_path, scope_raw, scope_receipt = lock(m.get("scope_receipt"), base, "scope receipt")
    scope = json_object(scope_raw, "scope receipt")
    bind(scope, {**binding, "schema": SCOPE}, "scope receipt")
    if scope.get("coverage") not in ("complete", "partial"):
        raise ValueError("scope receipt requires complete or partial coverage")
    scope_receipt["coverage"] = scope["coverage"]
    scope_receipt["explanation"] = name(scope.get("explanation"), "scope explanation")
    provenance_cache = {}
    scope_receipt["provenance"] = provenance(
        scope.get("provenance"), scope_path.parent, "scope", provenance_cache,
    )

    refs = m.get("edge_evidence")
    if not isinstance(refs, list):
        raise ValueError("edge_evidence must be a list")
    ledger = defaultdict(dict)
    loaded = set()
    for ref in refs:
        identity = reference_identity(ref, base, "edge evidence")
        if identity in loaded:
            continue
        loaded.add(identity)
        edge_id, row = read_evidence(ref, base, binding, edges, provenance_cache)
        digest = row["certificate"]["sha256"]
        previous = ledger[edge_id].get(digest)
        # Equivalent bytes at an alias never add a vote. Prefer the lexically first
        # path to make the provenance report independent of manifest ordering.
        if previous is None or row["certificate"]["path"] < previous["certificate"]["path"]:
            ledger[edge_id][digest] = row
    records = {edge: [by_sha[digest] for digest in sorted(by_sha)] for edge, by_sha in ledger.items()}
    decisions, evaluated = closure(seeds, edges, records, scope["coverage"] == "complete",
                                   round_budget, edge_budget)
    pending = set(decisions["unevaluated_frontier_edges"])
    edge_audit = []
    for edge_id in sorted(edges):
        edge = edges[edge_id]
        row = {"id": edge_id, "from": edge.source, "to": edge.destination,
               "evidence": [{key: value for key, value in record.items() if key != "_bridge_sequence"}
                            for record in records.get(edge_id, [])]}
        if edge_id in evaluated:
            row.update(evaluated[edge_id])
        elif edge_id in pending:
            row.update(state="unevaluated_frontier", reason="evaluation_budget")
        else:
            row.update(state="not_reached", reason="origin_outside_current_admitted_reachability")
        edge_audit.append(row)
    semantic = {
        "scope": scope_id, "eligibility_sha256": eligibility["sha256"],
        "read_pool_sha256": read_pool["sha256"], "scope_coverage": scope["coverage"], "result": decisions,
        "edges": [{"id": e["id"], "from": e["from"], "to": e["to"],
                   "state": e["state"], "reason": e["reason"]} for e in edge_audit],
    }
    report = {
        "schema": REPORT, "manifest_sha256": sha(raw), "graph": graph, "eligibility": eligibility,
        "read_pool": read_pool, "scope": scope_id, "scope_receipt": scope_receipt, "seed_set_sha256": binding["seed_set_sha256"],
        "round_budget": round_budget, "edge_budget": edge_budget, **decisions, "edges": edge_audit,
        "decision_sha256": canonical_sha(semantic),
        "inference_scope": "finite declared graph and supplied eligibility certificates",
        "source_evidence_recomputed": False,
    }

    if export_seeds:
        fasta, report["seed_export"] = seed_export(decisions, records)
        return report, fasta
    return report, None


def verify(manifest_path, round_budget=None, edge_budget=None):
    report, _ = replay(manifest_path, round_budget, edge_budget)
    return report


def publish(outputs):
    if len({path.resolve() for path, _ in outputs}) != len(outputs):
        raise ValueError("report and seed FASTA outputs must be distinct paths")
    for path, _ in outputs:
        if os.path.lexists(path):
            raise FileExistsError(f"output already exists: {path}")
    staged = []
    linked = []
    try:
        for path, content in outputs:
            with tempfile.NamedTemporaryFile(
                mode="wb", prefix=f".{path.name}.", dir=path.parent, delete=False,
            ) as stream:
                staged.append((path, Path(stream.name)))
                stream.write(content)
        for path, temporary in staged:
            os.link(temporary, path)
            linked.append((path, temporary))
    except BaseException:
        for path, temporary in reversed(linked):
            try:
                if os.path.samestat(path.lstat(), temporary.lstat()):
                    path.unlink()
            except FileNotFoundError:
                pass
        raise
    finally:
        for _, temporary in staged:
            temporary.unlink(missing_ok=True)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--out", required=True, type=Path, help="new report path")
    parser.add_argument("--out-seeds", type=Path, help="new FASTA of verified, admitted bridge sequences")
    parser.add_argument("--round-budget", type=int, help="maximum number of whole waves")
    parser.add_argument("--edge-budget", type=int, help="maximum edges evaluated in whole waves")
    args = parser.parse_args()
    try:
        result, fasta = replay(args.manifest.resolve(), args.round_budget, args.edge_budget,
                               export_seeds=args.out_seeds is not None)
        data = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode("utf-8")
        outputs = [(args.out, data)]
        if args.out_seeds is not None:
            # Publish the report last. Its seed digest attests the preceding FASTA.
            outputs.insert(0, (args.out_seeds, fasta))
        publish(outputs)
    except (OSError, ValueError) as error:
        parser.exit(1, f"{parser.prog}: {error}\n")


if __name__ == "__main__":
    main()
