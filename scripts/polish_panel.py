"""Polish a mitochondrial protein panel (metazoa or fungi): the fewest proteins that find each core gene in species of
genera the panel does not hold. Leave-genus-out greedy cover, per gene:

* pool: RefSeq mitochondrion genomes of the kingdom with >= 5 core genes (one per species), minus every genus in
  EXCLUDE_GENERA (the benchmark species, so the panel never holds a tested genus);
* tests: one genome per genus (seeded random, up to N_TEST), the panel's test of "a species it has not seen";
* a panel protein finds a test protein of the same gene when they share >= MIN_SHARED distinct amino-acid 6-mers
  (the read-side rule is >= 5 k-mer hits on the best gene, `panel::PanelParams::min_hits`) and the two are of different
  genera;
* greedy: add the candidate finding the most tests not yet found, until a pick adds fewer than MIN_GAIN tests.
Reports per gene: tests, found by the greedy panel, by the previous panel (PREV_FAA, scored the same way), proteins.
Writes OUT/<kind>_mito_proteins.faa (>mito|gene|Genus_species|accession), .tsv and polish.json.

usage: python3 -I polish_panel.py GBFF_GZ KIND OUT EXCLUDE_GENERA_FILE [PREV_FAA] [N_TEST] [MIN_SHARED] [MIN_GAIN] [THREADS]"""
import gzip
import json
import random
import re
import sys
from collections import defaultdict
from multiprocessing import Pool

sys.path.insert(0, __file__.rsplit("/", 1)[0])
gbff, kind, out, excl_file = sys.argv[1:5]
prev = sys.argv[5] if len(sys.argv) > 5 and sys.argv[5] != "-" else None
N_TEST = int(sys.argv[6]) if len(sys.argv) > 6 else 1500
MIN_SHARED = int(sys.argv[7]) if len(sys.argv) > 7 else 5
MIN_GAIN = int(sys.argv[8]) if len(sys.argv) > 8 else 1
THREADS = int(sys.argv[9]) if len(sys.argv) > 9 else 15
MAX_CAND = 6000
K = 6
KINGDOM = {"metazoa": "Metazoa", "fungi": "Fungi"}[kind]

from build_lineage_panels import norm, parse, records  # noqa: E402  (same parser and gene names)

exclude = {l.strip() for l in open(excl_file) if l.strip()}
genomes, seen_species = [], set()
for lines in records(gbff):
    acc, organism, lineage, cds = parse(lines)
    if not acc or KINGDOM not in lineage or not organism:
        continue
    genus = organism.split()[0]
    if genus in exclude or organism in seen_species or organism.split()[-1] in ("sp.",) or " sp. " in organism:
        continue
    prots = {}
    for q in cds:
        g = norm(q.get("gene"), q.get("product"))
        aa = q.get("translation", "").replace(" ", "")
        if g and len(aa) >= 40 and g not in prots:
            prots[g] = aa
    if len(prots) < 5:
        continue
    seen_species.add(organism)
    genomes.append({"acc": acc, "organism": organism, "genus": genus, "lineage": lineage, "prots": prots})
rng = random.Random(1)
by_genus = defaultdict(list)
for i, g in enumerate(genomes):
    by_genus[g["genus"]].append(i)
genera = sorted(by_genus)
rng.shuffle(genera)
tests = [rng.choice(by_genus[g]) for g in genera[:N_TEST]]
cand = list(range(len(genomes)))
if len(cand) > MAX_CAND:  # one per genus first, then random others
    first = [rng.choice(by_genus[g]) for g in sorted(by_genus)]
    rest = [i for i in cand if i not in set(first)]
    rng.shuffle(rest)
    cand = (first + rest)[:MAX_CAND] if len(first) < MAX_CAND else rng.sample(first, MAX_CAND)


def kmers(s):
    return {s[i:i + K] for i in range(len(s) - K + 1)}


def read_faa(path):
    out_, name, buf = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                out_.append((name, "".join(buf)))
            name, buf = line[1:], []
        else:
            buf.append(line)
    if name:
        out_.append((name, "".join(buf)))
    return out_


prev_prots = defaultdict(list)
if prev:
    for name, seq in read_faa(prev):
        f = name.split("|")
        prev_prots[f[1]].append((f[2].split("_")[0], seq))


import os
REDUNDANCY = tuple(int(x) for x in os.environ.get("REDUNDANCY", "1,2,3,5,8").split(","))


def sum_bits(cands):
    b = 0
    for _, x in cands:
        b |= x
    return b


def one_gene(gene):
    tids = [t for t in tests if gene in genomes[t]["prots"]]
    index = defaultdict(list)
    for j, t in enumerate(tids):
        for km in kmers(genomes[t]["prots"][gene]):
            index[km].append(j)
    tgenus = [genomes[t]["genus"] for t in tids]

    def found(seq, genus):
        cnt = defaultdict(int)
        for km in kmers(seq):
            for j in index.get(km, ()):
                cnt[j] += 1
        bits = 0
        for j, c in cnt.items():
            if c >= MIN_SHARED and tgenus[j] != genus:
                bits |= 1 << j
        return bits

    # even tests choose the panel, odd tests only score it (held out from the choice)
    sel = sum(1 << j for j in range(0, len(tids), 2))
    val = sum(1 << j for j in range(1, len(tids), 2))
    pc = lambda b: bin(b).count("1")
    cands = [(i, found(genomes[i]["prots"][gene], genomes[i]["genus"])) for i in cand if gene in genomes[i]["prots"]]

    def greedy(mask, r):
        """Fewest candidates so that every test in mask a candidate can find is found r times (or by all that can)."""
        level = [mask] + [0] * r          # level[k]: tests found at least k times
        picks, used = [], set()
        while True:
            need = mask & ~level[r]
            best, gain = None, 0
            for i, b in cands:
                if i in used:
                    continue
                g = pc(b & need)
                if g > gain or (g == gain and g and i < best):
                    best, gain, bits = i, g, b
            if best is None or gain < MIN_GAIN:
                return picks
            used.add(best)
            picks.append(best)
            for k in range(r, 0, -1):
                level[k] |= level[k - 1] & bits

    bits_of = dict(cands)
    rep = {"select_tests": pc(sel), "validate_tests": pc(val), "redundancy": {}}
    final = {}
    for r in REDUNDANCY:
        p_sel = greedy(sel, r)
        cov = 0
        for i in p_sel:
            cov |= bits_of[i]
        rep["redundancy"][r] = {"proteins": len(p_sel), "found": pc(cov & val)}
        final[r] = greedy(sel | val, r)   # the shipped panel: chosen on all tests
        rep["redundancy"][r]["proteins_all_tests"] = len(final[r])
    pbits = 0
    for genus, seq in prev_prots.get(gene, []):
        pbits |= found(seq, genus)
    rep["previous_found"] = pc(pbits & val) if prev else None
    rep["previous_proteins"] = len(prev_prots.get(gene, [])) if prev else None
    rep["findable"] = pc(val & sum_bits(cands))
    return gene, rep, final


genes = sorted({g for t in tests for g in genomes[t]["prots"]})
with Pool(min(THREADS, len(genes))) as pool:
    results = pool.map(one_gene, genes)
summary = {"kind": kind, "pool_genomes": len(genomes), "candidates": len(cand), "tests": len(tests),
           "excluded_genera": len(exclude), "min_shared_6mers": MIN_SHARED, "min_gain": MIN_GAIN, "genes": {}}
tot_val = sum(rep["validate_tests"] for _, rep, _ in results)
prev_found = sum(rep["previous_found"] or 0 for _, rep, _ in results)
by_r = {r: {"held_out_found": sum(rep["redundancy"][r]["found"] for _, rep, _ in results),
            "proteins": sum(rep["redundancy"][r]["proteins_all_tests"] for _, rep, _ in results)} for r in REDUNDANCY}
for r in REDUNDANCY:
    by_r[r]["held_out_rate"] = round(by_r[r]["held_out_found"] / tot_val, 4)
# the shipped redundancy: the smallest that finds as many held-out tests as the previous panel, else the best
ok = [r for r in REDUNDANCY if by_r[r]["held_out_found"] >= prev_found]
chosen = min(ok) if ok else max(REDUNDANCY, key=lambda r: by_r[r]["held_out_found"])
if os.environ.get("CHOSEN"):
    chosen = int(os.environ["CHOSEN"])  # a choice made from an earlier run of the curve
summary.update(redundancy=by_r, chosen_redundancy=chosen, validate_tests=tot_val, previous_found=prev_found,
               previous_rate=round(prev_found / tot_val, 4),
               findable=sum(rep["findable"] for _, rep, _ in results))
with open(f"{out}/{kind}_mito_proteins.faa", "w") as fa, open(f"{out}/{kind}_mito_proteins.tsv", "w") as tsv:
    tsv.write("gene\taccession\torganism\tlineage\n")
    for gene, rep, finals in results:
        summary["genes"][gene] = rep
        for i in finals[chosen]:
            g = genomes[i]
            sp = re.sub(r"[^A-Za-z0-9]+", "_", g["organism"]).strip("_")
            fa.write(f">mito|{gene}|{sp}|{g['acc']}\n")
            aa = g["prots"][gene]
            for k in range(0, len(aa), 80):
                fa.write(aa[k:k + 80] + "\n")
            tsv.write(f"{gene}\t{g['acc']}\t{g['organism']}\t{'; '.join(g['lineage'])}\n")
summary["previous_proteins"] = sum(r["previous_proteins"] or 0 for r in summary["genes"].values())
json.dump(summary, open(f"{out}/polish.json", "w"), indent=1)
print(json.dumps({k: v for k, v in summary.items() if k != "genes"}))
for g, r in summary["genes"].items():
    print(g, r["validate_tests"], "findable", r["findable"], "| previous", r["previous_found"], "with", r["previous_proteins"],
          "|", " ".join(f"r{k}:{v['found']}/{v['proteins_all_tests']}" for k, v in r["redundancy"].items()))
