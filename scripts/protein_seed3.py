"""Protein seeding for OVASM, v3 (11-05): protein_seed2.py plus a cluster call that needs enough hit reads.

1. Depth filter: only reads that `ovasm find discover` put into depth clusters (organelle-level depth, uniform along the
   read) are looked at, so nuclear insertions of organelle genes (NUMT/NUPT, at nuclear depth) never start a seed.
2. Protein start: up to SAMPLE reads per cluster are translated in six frames (standard code: land-plant organelles) and
   their amino-acid K-mers looked up in OVASM's organelle protein panel (data/organelle_proteins.faa). A read is a hit
   when it has >= MIN_HITS K-spaced hits on one organelle's proteins, more than on the other's.
3. Each cluster is called mitochondrion or plastid by the majority of its hit reads (depth separates the two organelles:
   a cluster is one depth level), and only when at least MIN_HIT_SHARE of its scanned reads are hits: an organelle
   cluster is full of gene-bearing reads, a deep repeat cluster is not (Polytrichum commune: a 6,317x cluster with 24
   mitochondrial hits in 3,000 reads seeded the mitochondrion and its recruitment precision fell to 0.09).
4. Seeds are only the hit stretch of each read (first to last hit) plus FLANK bp on each side, not the whole read.
5. Control: the same reads reversed (not complemented), scored and called the same way; reported per cluster.
Writes OUT/protein_seeds.{mitochondrion,plastid}.fasta and OUT/protein_seeds.json.

usage: python3 -I protein_seed2.py DISCOVER_JSON PROTEINS OUT [SAMPLE] [K] [MIN_HITS] [FLANK]
"""
import json
import sys
import time

disc_json, prot_path, out = sys.argv[1:4]
SAMPLE = int(sys.argv[4]) if len(sys.argv) > 4 else 3000
K = int(sys.argv[5]) if len(sys.argv) > 5 else 9
MIN_HITS = int(sys.argv[6]) if len(sys.argv) > 6 else 3
FLANK = int(sys.argv[7]) if len(sys.argv) > 7 else 2000
MAJORITY = 0.8          # share of a cluster's hit reads that must agree on the organelle
MIN_HIT_SHARE = 0.10    # share of a cluster's scanned reads that must be hits

CODE = {}
bases, aas = "TCAG", "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
i = 0
for a in bases:
    for b in bases:
        for c in bases:
            CODE[a + b + c] = aas[i]
            i += 1
COMP = str.maketrans("ACGTN", "TGCAN")

index = {}           # aa k-mer -> set of organelles


def add(name, seq):
    org = "mito" if name.split("|")[0].startswith("mit") else "pltd"
    for j in range(len(seq) - K + 1):
        index.setdefault(seq[j:j + K], set()).add(org)


name, buf = None, []
for line in open(prot_path):
    line = line.strip()
    if line.startswith(">"):
        if name and buf:
            add(name, "".join(buf))
        name, buf = line[1:], []
    elif line:
        buf.append(line.upper())
if name and buf:
    add(name, "".join(buf))


def hits(s):
    """K-spaced hits per organelle and the nucleotide span (first, last) they cover on s"""
    n = len(s)
    rc = s.translate(COMP)[::-1]
    count = {"mito": 0, "pltd": 0}
    span = {"mito": [n, -1], "pltd": [n, -1]}
    for strand, seq in ((1, s), (-1, rc)):
        for f in range(3):
            p = "".join(CODE.get(seq[j:j + 3], "X") for j in range(f, len(seq) - 2, 3))
            last = {"mito": -K, "pltd": -K}
            for j in range(len(p) - K + 1):
                orgs = index.get(p[j:j + K])
                if not orgs:
                    continue
                a = f + 3 * j
                b = a + 3 * K
                if strand == -1:
                    a, b = n - b, n - a
                for o in orgs:
                    if j - last[o] >= K:
                        count[o] += 1
                        last[o] = j
                    span[o][0] = min(span[o][0], a)
                    span[o][1] = max(span[o][1], b)
    return count, span


def call(count):
    for o, other in (("mito", "pltd"), ("pltd", "mito")):
        if count[o] >= MIN_HITS and count[o] > count[other]:
            return o
    return None


def fastq(path):
    with open(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().strip().upper()
            fh.readline()
            fh.readline()
            yield h[1:].split()[0], s


t0 = time.time()
disc = json.load(open(disc_json))
clusters, seeds = [], {"mito": [], "pltd": []}
for c in disc["clusters"]:
    step = max(1, c["kept_reads"] // SAMPLE)
    called = {"mito": [], "pltd": []}
    ctrl = {"mito": 0, "pltd": 0}
    scanned = 0
    for n, (rid, s) in enumerate(fastq(c["output"])):
        if n % step:
            continue
        scanned += 1
        cnt, span = hits(s)
        o = call(cnt)
        if o:
            a, b = span[o]
            called[o].append((rid, s[max(0, a - FLANK):min(len(s), b + FLANK)]))
        oc = call(hits(s[::-1])[0])
        if oc:
            ctrl[oc] += 1
        if scanned >= SAMPLE:
            break
    nm, npl = len(called["mito"]), len(called["pltd"])
    tot = nm + npl
    verdict = None
    if tot >= 3 and tot >= MIN_HIT_SHARE * scanned:
        if nm / tot >= MAJORITY:
            verdict = "mito"
        elif npl / tot >= MAJORITY:
            verdict = "pltd"
    if verdict:
        seeds[verdict].extend(called[verdict])
    clusters.append({"name": c["name"], "estimated_depth": round(c["estimated_depth"], 1), "reads": c["reads"],
                     "scanned": scanned, "hit_reads": {"mitochondrion": nm, "plastid": npl},
                     "reversed_control_hits": {"mitochondrion": ctrl["mito"], "plastid": ctrl["pltd"]},
                     "called": {"mito": "mitochondrion", "pltd": "plastid", None: None}[verdict]})
files = {"mito": "mitochondrion", "pltd": "plastid"}
for o, reads in seeds.items():
    with open(f"{out}/protein_seeds.{files[o]}.fasta", "w") as fh:
        for rid, s in reads:
            fh.write(f">{rid}\n{s}\n")
summary = {"design": "depth clusters (find discover) -> aa %d-mer hits (>= %d) on up to %d reads per cluster -> cluster "
                     "majority %.1f with hit share >= %.2f -> hit stretch +- %d bp" % (K, MIN_HITS, SAMPLE, MAJORITY, MIN_HIT_SHARE, FLANK),
           "clusters": clusters,
           "seed_reads": {files[o]: len(v) for o, v in seeds.items()},
           "seed_bases": {files[o]: sum(len(s) for _, s in v) for o, v in seeds.items()},
           "seconds": round(time.time() - t0, 1)}
json.dump(summary, open(f"{out}/protein_seeds.json", "w"), indent=1)
print(json.dumps(summary))
