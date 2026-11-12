#!/usr/bin/env python3
"""Count original records long enough to span a declared literal repeat route."""

from argparse import ArgumentParser
from collections import Counter
import json
from pathlib import Path

from discover_graph_frontier import fastx
from prepare_read_admission import LockedFile, digest, open_input, read_json, require
from validate_source_exclusion import Graph, new_output, publish_outputs


PLAN_SCHEMA = "ovasm.repeat-span-opportunity-plan.v1"
REPORT_SCHEMA = "ovasm.repeat-span-opportunity.v1"


def read_routes(entries, graph):
    require(isinstance(entries, list) and entries, "provide at least one repeat route")
    outgoing = graph.outgoing()
    routes = {}
    for entry in entries:
        require(isinstance(entry, dict) and set(entry) == {"id", "path"}, "invalid repeat route fields")
        name, path = entry["id"], entry["path"]
        require(isinstance(name, str) and name.strip() and name not in routes, "invalid repeat route ID")
        require(isinstance(path, list) and len(path) >= 3
                and all(isinstance(node, str) and node in graph.oriented for node in path),
                "a repeat route requires two exterior segments and at least one interior segment")
        overlaps = []
        for left, right in zip(path, path[1:]):
            options = [overlap for to, overlap in outgoing.get(left, ()) if to == right]
            require(len(options) == 1, f"route {name}: missing or ambiguous link {left} -> {right}")
            overlaps.append(options[0])
        # Flanking bases lie outside the interior spelling. Exterior overlaps locate
        # those bases and do not shorten the interior. Graph.parse checks overhangs.
        interior_bp = sum(len(graph.oriented[node]) for node in path[1:-1]) - sum(overlaps[1:-1])
        routes[name] = {"path": path, "overlaps_bp": overlaps,
                        "interior_spelling_bp": interior_bp, "minimum_literal_read_length": interior_bp + 2}
    return routes


def read_lengths(sources):
    lengths = Counter()
    for source in sources:
        with open_input(source) as stream:
            lengths.update(len(seq) for _, seq, _ in fastx(stream, preserve_original=True))
    return lengths


def audit(plan_path):
    root = plan_path.resolve().parent

    def lock(spec):
        return LockedFile.lock(root, spec)

    plan_file = lock({"path": str(plan_path.resolve()), "sha256": digest(plan_path.read_bytes())})
    plan = read_json(plan_file)
    require(isinstance(plan, dict) and set(plan) == {"schema", "scope", "graph", "raw_reads", "routes"}
            and plan["schema"] == PLAN_SCHEMA, "invalid repeat span plan schema or fields")
    require(isinstance(plan["scope"], str) and plan["scope"].strip(), "provide a scope")
    graph_file = lock(plan["graph"])
    with open_input(graph_file) as stream:
        graph = Graph.parse(stream.read())
    routes = read_routes(plan["routes"], graph)
    require(isinstance(plan["raw_reads"], list) and plan["raw_reads"], "provide original FASTX inputs")
    sources = [lock(entry) for entry in plan["raw_reads"]]
    require(len({source.state[:2] for source in sources}) == len(sources), "repeated original input file")
    lengths = read_lengths(sources)
    for file in (plan_file, graph_file, *sources):
        file.check()
    for route in routes.values():
        n = sum(count for length, count in lengths.items() if length >= route["minimum_literal_read_length"])
        route["raw_records_at_least_minimum"] = n
    return {
        "schema": REPORT_SCHEMA, "scope": plan["scope"], "plan": plan_file.binding(),
        "graph": graph_file.binding(), "raw_reads": [source.binding() for source in sources],
        "length_model": "Full literal interior spelling and one base from each exterior overhang.",
        "evidence_scope": "A length prefilter for declared routes. Read alignment and physical pairing are not evaluated.",
        "raw_records": sum(lengths.values()),
        "raw_bases": sum(length * count for length, count in lengths.items()),
        "maximum_read_length": max(lengths, default=None),
        "routes": routes,
    }


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    try:
        out = new_output(args.out)
        report = audit(args.plan)
        publish_outputs([(out, (json.dumps(report, indent=2, sort_keys=True) + "\n").encode())])
    except (OSError, ValueError, TypeError) as error:
        parser.exit(1, f"repeat span audit: {error}\n")
    print(f'{report["raw_records"]} original records, {len(report["routes"])} repeat routes')


if __name__ == "__main__":
    main()
