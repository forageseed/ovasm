"""Validate scoped source-exclusion edges against a frozen exact-spelling GFA."""
from argparse import ArgumentParser
from collections import defaultdict
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import os
import re
import tempfile

MANIFEST = "ovasm.wp2-d.source-exclusion.v1"
CERTIFICATE = "ovasm.wp2-d.certified-window.v1"
DNA = frozenset("ACGT")
COMP = str.maketrans("ACGT", "TGCA")
ROLES = {"competing", "protected", "forbidden"}


class Abstain(ValueError):
    """The input cannot authorize an exclusion."""


def sha(data):
    return hashlib.sha256(data).hexdigest()


def rc(sequence):
    return sequence.translate(COMP)[::-1]


def flip(node):
    return node[:-1]+("-" if node[-1]=="+" else "+")


def edge_key(a,b,overlap):
    return min((a,b,overlap),(flip(b),flip(a),overlap))


def load_locked(record, base, label):
    if not isinstance(record, dict) or not isinstance(record.get("path"), str):
        raise Abstain(f"{label}: missing path")
    expected=record.get("sha256")
    if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-f]{64}",expected):
        raise Abstain(f"{label}: missing SHA-256")
    path=(base/record["path"]).resolve()
    data=path.read_bytes()
    actual=sha(data)
    if actual!=expected:
        raise Abstain(f"{label}: SHA-256 mismatch")
    return path,data


@dataclass
class Graph:
    oriented: dict
    edges: set

    @classmethod
    def parse(cls,data):
        nodes, links={}, []
        for line in data.decode("utf-8").splitlines():
            if not line or line.startswith("#"):
                continue
            f=line.split("\t")
            if f[0] in {"H","P","W"}:
                continue
            if f[0]=="S":
                if len(f)<3 or not f[1] or f[1] in nodes or not f[2] or not set(f[2].upper())<=DNA:
                    raise Abstain("GFA has duplicate, missing, or non-ACGT segment sequence")
                nodes[f[1]]=f[2].upper()
            elif f[0]=="L":
                if len(f)<6 or f[2] not in {"+","-"} or f[4] not in {"+","-"}:
                    raise Abstain("GFA has an invalid oriented link")
                ov=re.fullmatch(r"([0-9]+)M",f[5])
                if not ov:
                    raise Abstain("GFA link does not specify an exact nM overlap")
                links.append((f[1]+f[2],f[3]+f[4],int(ov[1])))
            else:
                raise Abstain(f"GFA record {f[0]} may change the supported spelling model")
        if not nodes:
            raise Abstain("GFA has no segments")
        seqs={name+strand: q if strand=="+" else rc(q)
              for name,q in nodes.items() for strand in ("+","-")}
        physical, by_ports=set(), {}
        for a,b,ov in links:
            if a not in seqs or b not in seqs or ov>=min(len(seqs[a]),len(seqs[b])):
                raise Abstain("GFA has an unknown endpoint or non-progressing overlap")
            if ov and seqs[a][-ov:]!=seqs[b][:ov]:
                raise Abstain("GFA overlap disagrees with its oriented segment sequences")
            e=edge_key(a,b,ov)
            if e[:2] in by_ports and by_ports[e[:2]]!=ov:
                raise Abstain("GFA has distinct parallel overlaps for one physical port pair")
            by_ports[e[:2]]=ov
            physical.add(e)
        return cls(seqs,physical)

    def outgoing(self):
        out=defaultdict(set)
        for a,b,n in self.edges:
            out[a].add((b,n))
            out[flip(b)].add((flip(a),n))
        return {a:sorted(edges) for a,edges in out.items()}


def offsets(sequence,word):
    """Find all starts, including those whose seed crosses a segment boundary."""
    seed=word[:32]
    off=sequence.find(seed)
    while off>=0:
        yield off
        off=sequence.find(seed,off+1)
    for off in range(max(0,len(sequence)-len(seed)+1),len(sequence)):
        if seed.startswith(sequence[off:]):
            yield off


def route_signature(graph,walk,overlaps,start,stop):
    forward=(walk,overlaps,start,stop)
    reverse=(tuple(flip(n) for n in reversed(walk)),tuple(reversed(overlaps)),
             len(graph.oriented[walk[-1]])-stop,len(graph.oriented[walk[0]])-start)
    return min(forward,reverse)


def enumerate_routes(graph,word,state_cap,route_cap):
    adjacency=graph.outgoing()
    representatives={}
    states=starts=0

    def exhausted(reason):
        return {"complete":False,"reason":reason,"states":states,"starts_examined":starts,
                "routes":None}

    reverse=rc(word)
    directions=[("+",word)] if word==reverse else [("+",word),("-",reverse)]
    for strand,query in directions:
        for start_node,seq in sorted(graph.oriented.items()):
            for start in offsets(seq,query):
                starts+=1
                stack=[((start_node,),(),start,0)]
                while stack:
                    walk,ovs,lo,consumed=stack.pop()
                    states+=1
                    if states>state_cap:
                        return exhausted("state_cap")
                    node=walk[-1]
                    text=graph.oriented[node]
                    n=min(len(query)-consumed,len(text)-lo)
                    if text[lo:lo+n]!=query[consumed:consumed+n]:
                        continue
                    if consumed+n==len(query):
                        stop=lo+n
                        key=route_signature(graph,walk,ovs,start,stop)
                        if key not in representatives:
                            if len(representatives)>=route_cap:
                                return exhausted("route_cap")
                            representatives[key]={"walk":walk,"overlaps":ovs,"entry_offset0":start,
                                          "window_strand":strand,
                                                  "exit_offset0":stop,
                                                  "physical_edges":sorted({edge_key(a,b,o)
                                                        for a,b,o in zip(walk,walk[1:],ovs)})}
                        continue
                    assert n>0 and lo+n==len(text)
                    next_edges=adjacency.get(node,())
                    if len(stack)+len(next_edges)+states>state_cap:
                        return exhausted("state_cap")
                    for z,o in reversed(next_edges):
                        stack.append((walk+(z,),ovs+(o,),o,consumed+n))
    return {"complete":True,"reason":None,"states":states,"starts_examined":starts,
            "routes":[representatives[key] for key in sorted(representatives)]}


def keys(route):
    return {tuple(e) for e in route["physical_edges"]}


def require_windows(manifest,base,graph_sha,result):
    rows=manifest.get("windows")
    if manifest.get("protected_scope_state")!="complete_for_declared_scope":
        raise Abstain("protected panel has no completeness attestation for the declared scope")
    if not isinstance(rows,list) or not rows:
        raise Abstain("no certified windows")
    seen=set()
    for row in rows:
        name,role,word=row["id"],row["role"],row["sequence"]
        if (not isinstance(name,str) or not name or name in seen or role not in ROLES
            or not isinstance(word,str) or not word or not set(word)<=DNA):
            raise Abstain("window has invalid or duplicate identity, role, or exact sequence")
        seen.add(name)
        digest=sha(word.encode("ascii"))
        cp,body=load_locked(row.get("certificate"),base,f"window {name} certificate")
        cert=json.loads(body)
        if (cert.get("schema"),cert.get("graph_sha256"),cert.get("window_id"),cert.get("role"),
            cert.get("window_sha256"),cert.get("full_window") is True) != (
                CERTIFICATE,graph_sha,name,role,digest,True):
            raise Abstain(f"{name}: certificate does not bind the full window to this graph")
        proofs=cert.get("provenance")
        if not isinstance(proofs,list) or not proofs:
            raise Abstain(f"{name}: certificate has no independent provenance")
        for ref in proofs:
            load_locked(ref,cp.parent,f"{name} provenance")
        required="derived" if role=="forbidden" else "certified"
        state=cert.get("status")
        result["windows"][name]={"role":role,"bp":len(word),"sha256":digest,
                                  "certificate_sha256":sha(body),"certification_status":state}
        if state!=required:
            raise Abstain(f"{name}: full source/target context is unresolved or uncertified")
        yield name,role,word


def validate(manifest_path,state_cap=200_000,route_cap=4096):
    result={"schema":"ovasm.wp2-d.source-exclusion-audit.v1","state":"abstained","reason":None,
            "safe_edges":[],"windows":{},"caps":{"states_per_window":state_cap,"routes_per_window":route_cap}}
    try:
        if state_cap<1 or route_cap<1:
            raise Abstain("enumeration caps must be positive")
        raw=manifest_path.read_bytes()
        m=json.loads(raw)
        result["manifest_sha256"]=sha(raw)
        scope=m.get("scope")
        if m.get("schema")!=MANIFEST or not isinstance(scope,str) or not scope.strip():
            raise Abstain("unknown manifest schema or missing local scope")
        result["scope"]=scope
        path,data=load_locked(m.get("graph"),manifest_path.parent,"frozen graph")
        graph_sha=sha(data)
        result["graph"]={"path":str(path),"sha256":graph_sha}
        graph=Graph.parse(data)
        windows=list(require_windows(m,manifest_path.parent,graph_sha,result))
        if not {"competing","protected"}<={role for _,role,_ in windows}:
            raise Abstain("need independently certified competing and protected full windows")
        competing_words={min(w,rc(w)) for _,role,w in windows if role=="competing"}
        protected_words={min(w,rc(w)) for _,role,w in windows if role=="protected"}
        if competing_words&protected_words:
            raise Abstain("same full spelling certified for competing and protected source roles")
        for name,role,word in windows:
            result["windows"][name]["original"]=enumerate_routes(graph,word,state_cap,route_cap)
        for name,role,_ in windows:
            orig=result["windows"][name]["original"]
            if not orig["complete"]:
                raise Abstain(f"{name}: exact route enumeration exceeded {orig['reason']}")
            count=len(orig["routes"])
            if role!="forbidden" and count!=1:
                raise Abstain(f"{name}: expected one exact physical placement, found {count}")
            if role=="forbidden" and count>1:
                raise Abstain(f"{name}: forbidden full window has {count} ambiguous exact placements")
        protected, competing=set(),set()
        for name,role,_ in windows:
            if role=="forbidden":
                continue
            route,=result["windows"][name]["original"]["routes"]
            required=keys(route)
            if role=="competing":
                if not required:
                    raise Abstain(f"{name}: competing spelling has no required graph edge")
                competing.update(required)
            else:
                protected.update(required)
        result["protected_edges"]=sorted(protected)
        result["competing_edges"]=sorted(competing)
        result["competing_protected_edge_overlap"]=sorted(competing&protected)
        if competing&protected:
            raise Abstain("competing and protected physical routes share an edge")
        reduced=Graph(graph.oriented,graph.edges-competing)
        for name,role,word in windows:
            entry=result["windows"][name]
            entry["masked"]=enumerate_routes(reduced,word,state_cap,route_cap)
        for name,role,_ in windows:
            entry=result["windows"][name]
            masked=entry["masked"]
            if not masked["complete"]:
                raise Abstain(f"{name}: masked graph enumeration exceeded {masked['reason']}")
            if role=="protected" and masked["routes"]!=entry["original"]["routes"]:
                raise Abstain(f"{name}: protected exact placement changed")
            if role!="protected" and masked["routes"]:
                raise Abstain(f"{name}: competing or forbidden exact spelling survives the proposed mask")
        view_payload=json.dumps([graph_sha,sorted(reduced.edges)],separators=(",",":")).encode()
        result["masked_physical_view_sha256"]=sha(view_payload)
        result["safe_edges"]=sorted(competing)
        result["state"]="eligible_scoped_graph_exclusion"
    except (Abstain, OSError, ValueError, TypeError, KeyError, AttributeError) as error:
        result["reason"]=str(error)
    return result


def materialize_view(result):
    """Spell the audited edge-only view without rewriting unrelated GFA records."""
    graph_record = result["graph"]
    _, data = load_locked(graph_record, Path("."), "frozen graph at publication")
    graph = Graph.parse(data)
    blocked = {tuple(edge) for edge in result["safe_edges"]}
    remaining = graph.edges - blocked
    if not blocked or not blocked <= graph.edges:
        raise Abstain("verified exclusion no longer matches the frozen graph")
    view_payload = json.dumps(
        [graph_record["sha256"], sorted(remaining)], separators=(",", ":")
    ).encode()
    if sha(view_payload) != result["masked_physical_view_sha256"]:
        raise Abstain("physical view identity changed before publication")

    kept = []
    removed = 0
    for line in data.splitlines(keepends=True):
        fields = line.rstrip(b"\r\n").decode("utf-8").split("\t")
        kind = fields[0]
        # P and W annotate paths already chosen on the input topology. A raw view
        # must not publish potentially stale path annotations.
        if kind in {"P", "W"}:
            raise Abstain("materializing a view requires a raw GFA without P or W records")
        if kind == "L":
            a, b = fields[1] + fields[2], fields[3] + fields[4]
            if edge_key(a, b, int(fields[5][:-1])) in blocked:
                removed += 1
                continue
        kept.append(line)
    view = b"".join(kept)
    parsed = Graph.parse(view)
    if parsed.oriented != graph.oriented or parsed.edges != remaining:
        raise Abstain("serialized view differs from the audited graph")
    return view, removed


def new_output(path):
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"output already exists: {path}")
    return path.resolve()


def publish_outputs(outputs):
    """Stage bytes before publication and never replace an existing file."""
    staged = []
    published = []
    try:
        for path, data in outputs:
            with tempfile.NamedTemporaryFile(
                mode="wb", prefix=f".{path.name}.", dir=path.parent, delete=False
            ) as stream:
                temporary = Path(stream.name)
                staged.append((temporary, path))
                stream.write(data)
                stream.flush()
                os.fsync(stream.fileno())
        for temporary, path in staged:
            os.link(temporary, path)
            published.append((temporary, path))
    except BaseException:
        for temporary, path in reversed(published):
            # An output that another process replaced is outside this transaction.
            if not path.is_symlink() and path.is_file() and temporary.samefile(path):
                path.unlink()
        raise
    finally:
        for temporary, _ in staged:
            temporary.unlink(missing_ok=True)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--out", required=True, type=Path, help="New JSON audit path.")
    parser.add_argument(
        "--out-gfa", type=Path,
        help="New raw GFA path for the certified view, omitted when the audit abstains.",
    )
    parser.add_argument("--state-cap", type=int, default=200_000)
    parser.add_argument("--route-cap", type=int, default=4096)
    args = parser.parse_args()
    try:
        audit_path = new_output(args.out)
        view_path = new_output(args.out_gfa) if args.out_gfa else None
        if audit_path == view_path:
            raise ValueError("audit and graph view need different output paths")
        result = validate(args.manifest.resolve(), args.state_cap, args.route_cap)
        outputs = []
        if view_path is not None and result["state"] == "eligible_scoped_graph_exclusion":
            try:
                view, removed = materialize_view(result)
            except (Abstain, OSError, ValueError) as error:
                result.update(state="abstained", reason=str(error), safe_edges=[])
            else:
                result["graph_view"] = {
                    "path": str(view_path), "sha256": sha(view),
                    "removed_link_records": removed,
                }
                outputs.append((view_path, view))
        audit = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode("utf-8")
        outputs.append((audit_path, audit))
        publish_outputs(outputs)
    except (OSError, ValueError) as error:
        parser.exit(1, f"{parser.prog}: {error}\n")


if __name__ == "__main__":
    main()
