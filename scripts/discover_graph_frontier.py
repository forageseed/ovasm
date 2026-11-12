#!/usr/bin/env python3
"""Enumerate source-neutral raw sequence differences from a SHA-bound GFA1 graph."""

from argparse import ArgumentParser
from collections import Counter, defaultdict
from contextlib import ExitStack, contextmanager
from hashlib import file_digest, sha256
from pathlib import Path
import gzip
import json
import os
import re
import tempfile


MANIFEST = "ovasm.graph-frontier-discovery.v1"
REPORT = "ovasm.graph-frontier-discovery-report.v1"
DEFAULT_TRAVERSAL_STATES = 1_000_000
MAX_SCAN_K = 31
COMPLEMENT = bytes.maketrans(b"ACGTRYSWKMBDHVN", b"TGCAYRSWMKVHDBN")
AMBIGUOUS = re.compile(rb"[^ACGT]")
MISSING = object()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def json_object(raw, label):
    def unique(items):
        result = {}
        for key, value in items:
            require(key not in result, f"{label}: duplicate key {key!r}")
            result[key] = value
        return result

    result = json.loads(raw, object_pairs_hook=unique)
    require(isinstance(result, dict), f"{label}: expected a JSON object")
    return result


def fingerprint(info):
    return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns


@contextmanager
def locked_file(entry, base, label):
    require(isinstance(entry, dict) and isinstance(entry.get("path"), str) and entry["path"],
            f"{label}: expected a nonempty file path")
    expected = entry.get("sha256")
    require(isinstance(expected, str) and re.fullmatch(r"[0-9a-f]{64}", expected),
            f"{label}: expected lowercase SHA-256")
    path = (base / entry["path"]).resolve(strict=True)
    require(path.is_file(), f"{label}: expected a regular input file: {path}")
    with path.open("rb") as stream:
        info = os.fstat(stream.fileno())
        snapshot = fingerprint(info)
        require(file_digest(stream, "sha256").hexdigest() == expected, f"{label}: SHA-256 mismatch: {path}")
        stream.seek(0)
        yield {"path": str(path), "sha256": expected, "bytes": info.st_size}, stream, snapshot[:2]
        require(snapshot == fingerprint(os.fstat(stream.fileno())) == fingerprint(path.stat()),
                f"{label}: input changed while scanning: {path}")


def revcomp(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def flip(node):
    return node[:-1] + ("-" if node[-1] == "+" else "+")


def parse_graph(stream):
    sequences = {}
    links = []
    for number, line in enumerate(stream, 1):
        line = line.rstrip(b"\r\n")
        if not line or line.startswith(b"#"):
            continue
        columns = line.split(b"\t")
        tag = columns[0]
        if tag in (b"H", b"P", b"W"):
            continue
        if tag == b"S":
            require(len(columns) >= 3, f"GFA line {number}: incomplete segment")
            name = columns[1].decode("utf-8")
            sequence = columns[2].upper()
            require(name and not any(c.isspace() for c in name) and name + "+" not in sequences,
                    f"GFA line {number}: invalid or duplicate segment")
            require(re.fullmatch(rb"[ACGTRYSWKMBDHVN]+", sequence),
                    f"GFA line {number}: expected a nonempty inline nucleotide sequence")
            sequences[name + "+"] = sequence
            sequences[name + "-"] = revcomp(sequence)
        elif tag == b"L":
            require(len(columns) >= 6 and columns[2] in (b"+", b"-") and columns[4] in (b"+", b"-"),
                    f"GFA line {number}: invalid oriented link")
            overlap = re.fullmatch(rb"(\d+)M", columns[5])
            require(overlap is not None, f"GFA line {number}: only exact nM overlaps are supported")
            links.append(((columns[1] + columns[2]).decode("utf-8"),
                          (columns[3] + columns[4]).decode("utf-8"), int(overlap[1]), number))
        else:
            raise ValueError(f"GFA line {number}: unsupported graph record {tag!r}")
    require(sequences, "GFA has no segments")
    adjacency = {node: set() for node in sequences}
    for left, right, overlap, number in links:
        require(left in sequences and right in sequences, f"GFA line {number}: unknown segment")
        a, b = sequences[left], sequences[right]
        require(overlap <= min(len(a), len(b)) and (overlap == 0 or a[-overlap:] == b[:overlap]),
                f"GFA line {number}: overlap disagrees with segment sequences")
        adjacency[left].add((right, overlap))
        adjacency[flip(right)].add((flip(left), overlap))
    return sequences, {node: sorted(edges) for node, edges in adjacency.items()}


def graph_words(sequences, adjacency, k, state_budget):
    words = {}
    for node, sequence in sorted(sequences.items()):
        for offset in range(len(sequence) - k + 1):
            word = sequence[offset:offset + k]
            if not AMBIGUOUS.search(word):
                words[word] = None if word in words else (node, offset)
    states = 0
    for node, sequence in sorted(sequences.items()):
        if not adjacency[node]:
            continue
        for length in range(1, min(k - 1, len(sequence)) + 1):
            tail = sequence[-length:]
            if AMBIGUOUS.search(tail):
                continue
            if states == state_budget:
                return words, states, False
            states += 1
            start = (node, len(sequence), tail)
            todo = [start]
            visited = {start}
            while todo:
                current, offset, prefix = todo.pop()
                extension = sequences[current][offset:offset + (k - len(prefix))]
                if AMBIGUOUS.search(extension):
                    continue
                word = prefix + extension
                if len(word) == k:
                    # This occurrence consumes novel bases on another oriented node.
                    # It invalidates a node-only uniqueness claim even if it is an alias.
                    words[word] = None
                    continue
                for target, overlap in reversed(adjacency[current]):
                    state = (target, overlap, word)
                    if state not in visited:
                        if states == state_budget:
                            return words, states, False
                        states += 1
                        visited.add(state)
                        todo.append(state)
    return words, states, True


def graph_position(node, offset, size, sequences):
    length = len(sequences[node])
    forward = [offset, offset + size] if node[-1] == "+" else [length - offset - size, length - offset]
    return {
        "segment": node[:-1],
        "orientation": node[-1],
        "oriented_interval0": [offset, offset + size],
        "forward_segment_interval0": forward,
    }


def anchor(raw_start, position, sequences, k):
    node, offset = position
    return dict(graph_position(node, offset, k, sequences), raw_interval0=[raw_start, raw_start + k])


def graph_topology(sequences, adjacency, k):
    ports = []
    links = []
    for node, sequence in sorted(sequences.items()):
        edges = adjacency[node]
        if len(edges) != 1:
            size = min(k, len(sequence))
            ports.append(dict(graph_position(node, len(sequence) - size, size, sequences),
                              type="terminal" if not edges else "branch", outdegree=len(edges),
                              sequence=sequence[-size:].decode("ascii")))
        for target, overlap in edges:
            links.append(dict(source=node, target=target, overlap_bp=overlap,
                              left_suffix=sequence[-k:].decode("ascii"),
                              right_after_overlap=sequences[target][overlap:overlap + k].decode("ascii")))
    return dict(ports=ports, oriented_links=links)


def observable_runs(sequence, k):
    start = 0
    for match in AMBIGUOUS.finditer(sequence):
        stop = match.start() - k + 1
        if start < stop:
            yield start, stop
        start = match.end()
    stop = len(sequence) - k + 1
    if start < stop:
        yield start, stop


def events(sequence, words, sequences, adjacency, k):
    sequence = sequence.upper()
    word_count = max(0, len(sequence) - k + 1)
    observable = represented = unique = 0
    found = []
    for run_start, run_stop in observable_runs(sequence, k):
        observable += run_stop - run_start
        present = []
        for offset in range(run_start, run_stop):
            hit = words.get(sequence[offset:offset + k], MISSING)
            if hit is not MISSING:
                present.append((offset, hit))
        represented += len(present)
        unique += sum(hit is not None for _, hit in present)
        if not present:
            continue
        bounds = [(run_start - 1, None), *present, (run_stop, None)]
        for (start, left), (stop, right) in zip(bounds, bounds[1:]):
            if stop <= start + 1 or left is None and right is None:
                continue
            left_port = left is not None and left[1] + k == len(sequences[left[0]]) and len(adjacency[left[0]]) != 1
            right_port = right is not None and right[1] == 0 and len(adjacency[flip(right[0])]) != 1
            two = left is not None and right is not None
            noncollinear = two and (left[0] != right[0] or right[1] - left[1] != stop - start)
            found.append(dict(
                absent_word_start_interval0=[start + 1, stop],
                affected_raw_interval0=[start + 1, stop + k - 1],
                left_anchor=anchor(start, left, sequences, k) if left is not None else None,
                right_anchor=anchor(stop, right, sequences, k) if right is not None else None,
                same_read_two_unique_anchors=two,
                two_unique_noncollinear_oriented_anchors=noncollinear,
                terminal_or_branch_port=left_port or right_port,
            ))
    return found, represented, unique, observable, word_count - observable


def fastx(stream, *, preserve_original=False):
    def without_line_ending(line):
        if preserve_original:
            return line.removesuffix(b"\n").removesuffix(b"\r")
        return line.rstrip(b"\r\n")

    marker = stream.peek(1)[:1]
    if marker == b"@":
        for header in stream:
            sequence, plus, quality = (without_line_ending(stream.readline()) for _ in range(3))
            identifier = without_line_ending(header[1:])
            require(header.startswith(b"@") and identifier.strip() and sequence and plus.startswith(b"+")
                    and len(quality) == len(sequence), "invalid four-line FASTQ record")
            yield identifier, sequence, quality
    elif marker == b">":
        identifier = None
        pieces = []
        for line in stream:
            if line.startswith(b">"):
                if identifier is not None:
                    require(pieces, "empty FASTA sequence")
                    yield identifier, b"".join(pieces), None
                identifier = without_line_ending(line[1:])
                require(identifier.strip(), "empty FASTA identifier")
                pieces = []
            else:
                require(identifier is not None, "FASTA sequence before first identifier")
                # Native FASTA records remove CR/LF and keep other bytes.
                piece = line.translate(None, b"\r\n") if preserve_original else line.strip()
                if piece:
                    pieces.append(piece)
        if identifier is not None:
            require(pieces, "empty FASTA sequence")
            yield identifier, b"".join(pieces), None
    else:
        require(not marker, "input must be FASTA or four-line FASTQ")


def raw_records(source):
    if source.peek(2).startswith(b"\x1f\x8b"):
        with gzip.GzipFile(fileobj=source, mode="rb") as stream:
            yield from fastx(stream)
    else:
        yield from fastx(source)


EVENT_FLAGS = ("same_read_two_unique_anchors", "two_unique_noncollinear_oriented_anchors", "terminal_or_branch_port")


def scan_input(stream, output, input_index, words, sequences, adjacency, k, include_no_graph):
    counts = Counter(dict(records_scanned=0, no_graph_words=0, shorter_than_scan_k=0,
                          any_graph_words=0, any_unique_node_word=0, candidate_records=0,
                          records_with_graph_words_no_candidate=0, output_records=0, candidate_events=0,
                          no_observable_words=0, records_with_unobservable_words=0,
                          unobservable_word_positions=0, observable_word_positions=0,
                          no_graph_words_with_complete_observation=0))
    for flag in EVENT_FLAGS:
        counts[flag + "_records"] = 0
    for ordinal, (identifier, sequence, quality) in enumerate(raw_records(stream)):
        counts["records_scanned"] += 1
        counts["shorter_than_scan_k"] += len(sequence) < k
        gaps, represented, unique, observable, unobservable = events(sequence, words, sequences, adjacency, k)
        counts["observable_word_positions"] += observable
        counts["unobservable_word_positions"] += unobservable
        counts["records_with_unobservable_words"] += unobservable > 0
        counts["no_observable_words"] += observable == 0
        counts["any_graph_words"] += represented > 0
        counts["any_unique_node_word"] += unique > 0
        counts["no_graph_words"] += represented == 0 and observable > 0
        counts["no_graph_words_with_complete_observation"] += represented == 0 and observable > 0 and unobservable == 0
        counts["records_with_graph_words_no_candidate"] += represented > 0 and not gaps
        if not gaps and (represented or not include_no_graph):
            continue
        counts["output_records"] += 1
        counts["candidate_records"] += bool(gaps)
        counts["candidate_events"] += len(gaps)
        for flag in EVENT_FLAGS:
            counts[flag + "_records"] += any(event[flag] for event in gaps)
        row = dict(
            input_index0=input_index, record_index0=ordinal, read_id_hex=identifier.hex(), read_length=len(sequence),
            sequence_sha256=sha256(sequence).hexdigest(),
            quality_sha256=sha256(quality).hexdigest() if quality is not None else None,
            status="graph_anchored_candidate" if gaps else "no_graph_words" if observable else "no_observable_words",
            represented_word_positions=represented, unique_node_word_positions=unique, events=gaps,
            observable_word_positions=observable, unobservable_word_positions=unobservable,
        )
        output.write(json.dumps(row, separators=(",", ":"), sort_keys=True).encode() + b"\n")
    return dict(counts)


def evaluate(manifest_path, evidence, state_budget, include_no_graph):
    raw = manifest_path.read_bytes()
    manifest = json_object(raw, "manifest")
    require(manifest.get("schema") == MANIFEST, "unexpected discovery manifest schema")
    scope = manifest.get("scope")
    require(isinstance(scope, str) and scope.strip(), "manifest scope is required")
    require(("recruit" in manifest) != ("scan_k" in manifest), "specify exactly one of recruit or scan_k")
    base = manifest_path.parent
    if "recruit" in manifest:
        with locked_file(manifest["recruit"], base, "recruit report") as (recruit, stream, _):
            data = json_object(stream.read(), "recruit report")
            require(data.get("tool") == "ovasm find recruit", "recruit report tool does not match OVASM recruitment")
            k = data.get("k")
    else:
        recruit = None
        k = manifest["scan_k"]
    require(type(k) is int and 1 <= k <= MAX_SCAN_K, f"scan k must be an integer in [1, {MAX_SCAN_K}]")
    require(type(state_budget) is int and state_budget >= 0, "traversal state budget must be nonnegative")
    with locked_file(manifest.get("graph"), base, "graph") as (graph, stream, _):
        sequences, adjacency = parse_graph(stream)
        words, states, complete = graph_words(sequences, adjacency, k, state_budget)
    sources = manifest.get("raw_reads")
    require(isinstance(sources, list) and sources, "at least one raw_reads file is required")
    reads = []
    identities = set()
    for index, item in enumerate(sources):
        with locked_file(item, base, f"raw_reads[{index}]") as (receipt, stream, identity):
            require(identity not in identities, "duplicate raw input file or hardlink")
            identities.add(identity)
            stats = scan_input(stream, evidence, index, words, sequences, adjacency, k, include_no_graph) if complete else None
            reads.append(dict(receipt, counts=stats))
    return dict(
        schema=REPORT, scope=scope, manifest_sha256=sha256(raw).hexdigest(), graph=graph, recruit=recruit,
        scan_k=k, raw_reads=reads, graph_topology=graph_topology(sequences, adjacency, k),
        graph_index=dict(status="complete" if complete else "truncated",
                         traversal_states=states, traversal_state_budget=state_budget,
                         represented_oriented_words=len(words) if complete else None,
                         unique_oriented_node_words=sum(x is not None for x in words.values()) if complete else None),
        raw_scan_status="complete" if complete else "not_run_incomplete_graph_index",
        no_graph_word_rows_emitted=include_no_graph and complete,
        coordinate_system="all intervals are zero-based and half-open",
        scope_limits="Graph-compatible words and raw differences only. No reference or source classification. "
                     "Graph occurrences spanning edges invalidate node-only unique anchors. FASTX counts are records.",
    )


@contextmanager
def staged(path):
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        tmp = Path(stream.name)
        try:
            yield tmp, stream
        finally:
            tmp.unlink(missing_ok=True)


def publish(report_path, evidence_path, manifest_path, budget, include_no_graph):
    paths = [report_path, evidence_path]
    require(len({p.resolve() for p in paths}) == len(paths), "output paths must differ")
    for path in paths:
        require(not (path.exists() or path.is_symlink()), f"output already exists: {path}")
    with staged(evidence_path) as (evidence_tmp, evidence), staged(report_path) as (report_tmp, report_stream):
        if evidence_path.suffix == ".gz":
            with gzip.GzipFile(filename="", fileobj=evidence, mode="wb", compresslevel=4, mtime=0) as output:
                report = evaluate(manifest_path, output, budget, include_no_graph)
        else:
            report = evaluate(manifest_path, evidence, budget, include_no_graph)
        evidence.flush()
        os.fsync(evidence.fileno())
        with evidence_tmp.open("rb") as source:
            report["evidence"] = dict(path=str(evidence_path.resolve()), bytes=evidence_tmp.stat().st_size,
                                      sha256=file_digest(source, "sha256").hexdigest())
        report_stream.write(json.dumps(report, indent=2, sort_keys=True).encode() + b"\n")
        report_stream.flush()
        os.fsync(report_stream.fileno())
        linked = []
        try:
            for temp, path in ((evidence_tmp, evidence_path), (report_tmp, report_path)):
                os.link(temp, path)
                linked.append(path)
        except BaseException:
            for path in reversed(linked):
                path.unlink()
            raise


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--out", type=Path, required=True, help="new JSON report")
    parser.add_argument("--out-evidence", type=Path, required=True, help="new JSONL or .gz evidence file")
    parser.add_argument("--max-traversal-states", type=int, default=DEFAULT_TRAVERSAL_STATES)
    parser.add_argument("--include-no-graph-words", action="store_true",
                        help="also write each record with no represented or no observable graph words")
    args = parser.parse_args()
    try:
        publish(args.out, args.out_evidence, args.manifest, args.max_traversal_states, args.include_no_graph_words)
    except (OSError, ValueError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
