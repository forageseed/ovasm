"""Sensitivity and chance hits of a mitochondrial protein panel at several hit thresholds, scored the way panel.rs
assigns a stretch (k-mer hits against the union of a gene's panel proteins):
* sensitivity: held-out test genomes (one per genus, seed 1, the tests of polish_panel.py); a test gene counts as found at
  threshold T when >= T of its protein's 6-mers are in the gene's index built from panel proteins of other genera;
* chance: random DNA (GC as given), six frames, stop-free stretches >= 80 aa; a stretch is a false hit at T when its best
  gene has >= T hits. Reported per Mb and per 15 kb read.
usage: python3 -I panel_tradeoff.py GBFF_GZ KIND PANEL_FAA TABLE EXCLUDE_GENERA [N_TEST] [MB] [GC]"""
import random
import sys
from collections import defaultdict

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from build_lineage_panels import norm, parse, records  # noqa: E402

gbff, kind, faa, table, excl = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
N_TEST = int(sys.argv[6]) if len(sys.argv) > 6 else 1500
MB = float(sys.argv[7]) if len(sys.argv) > 7 else 3
GC = float(sys.argv[8]) if len(sys.argv) > 8 else 0.3
TS = (5, 6, 7, 8, 10, 12, 15)
K = 6
KINGDOM = {"metazoa": "Metazoa", "fungi": "Fungi"}[kind]

# tests exactly as polish_panel.py draws them
exclude = {l.strip() for l in open(excl) if l.strip()}
genomes, seen = [], set()
for lines in records(gbff):
    acc, organism, lineage, cds = parse(lines)
    if not acc or KINGDOM not in lineage or not organism:
        continue
    genus = organism.split()[0]
    if genus in exclude or organism in seen or organism.split()[-1] in ("sp.",) or " sp. " in organism:
        continue
    prots = {}
    for q in cds:
        g = norm(q.get("gene"), q.get("product"))
        aa = q.get("translation", "").replace(" ", "")
        if g and len(aa) >= 40 and g not in prots:
            prots[g] = aa
    if len(prots) < 5:
        continue
    seen.add(organism)
    genomes.append({"genus": genus, "prots": prots})
rng = random.Random(1)
by_genus = defaultdict(list)
for i, g in enumerate(genomes):
    by_genus[g["genus"]].append(i)
genera = sorted(by_genus)
rng.shuffle(genera)
tests = [rng.choice(by_genus[g]) for g in genera[:N_TEST]]

index = defaultdict(lambda: defaultdict(set))   # gene -> kmer -> genera
any_gene = defaultdict(set)                     # kmer -> genes
name, buf = None, []


def add(n, s):
    f = n.split("|")
    gene, genus = f[1], f[2].split("_")[0]
    for i in range(len(s) - K + 1):
        index[gene][s[i:i + K]].add(genus)
        any_gene[s[i:i + K]].add(gene)


for line in open(faa):
    line = line.strip()
    if line.startswith(">"):
        if name:
            add(name, "".join(buf))
        name, buf = line[1:], []
    else:
        buf.append(line)
if name:
    add(name, "".join(buf))
found = {t: 0 for t in TS}
total = 0
for t in tests:
    g = genomes[t]
    for gene, aa in g["prots"].items():
        if gene not in index:
            continue
        total += 1
        idx = index[gene]
        hits = sum(1 for i in range(len(aa) - K + 1) if (idx.get(aa[i:i + K]) or set()) - {g["genus"]})
        for T in TS:
            found[T] += hits >= T

B = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
code = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(B) for j, b in enumerate(B) for k, c in enumerate(B)}
if table in (2, 4, 5):
    code["TGA"] = "W"
    if table != 4:
        code["ATA"] = "M"
        code["AGA"] = code["AGG"] = "S" if table == 5 else "*"
r2 = random.Random(7)
n = int(MB * 1e6)
seq = "".join(r2.choices("ACGT", weights=[(1 - GC) / 2, GC / 2, GC / 2, (1 - GC) / 2], k=n))
chance = {t: 0 for t in TS}
comp = str.maketrans("ACGT", "TGCA")
for s in (seq, seq.translate(comp)[::-1]):
    for f in range(3):
        prot = "".join(code[s[i:i + 3]] for i in range(f, len(s) - 2, 3))
        for orf in prot.split("*"):
            if len(orf) < 80:
                continue
            hits = defaultdict(int)
            for i in range(len(orf) - K + 1):
                for gene in any_gene.get(orf[i:i + K], ()):
                    hits[gene] += 1
            best = max(hits.values()) if hits else 0
            for T in TS:
                chance[T] += best >= T
print(f"# {faa} table {table}: {total} held-out test genes; random {MB} Mb at GC {GC}")
print("T\tsensitivity\tchance_per_Mb\tchance_per_15kb_read")
for T in TS:
    print(f"{T}\t{found[T] / total:.4f}\t{chance[T] / MB:.2f}\t{chance[T] / MB * 0.015:.4f}")
