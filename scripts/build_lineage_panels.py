"""Mitochondrial protein panels for OVASM's metazoa and fungi protein seeds, from the RefSeq mitochondrion release
(mitochondrion.1.genomic.gbff.gz). Only the core respiratory-chain and ATP-synthase genes (plus rps3) are kept; their
CDS translations come from RefSeq (each record's own genetic code). Representatives: per taxonomic group (the first
DEPTH lineage terms) up to PER_GROUP genomes, preferring genomes with the most core genes.
Headers: >mito|GENE|Genus_species|ACCESSION  (the panel format class|gene|species).

usage: python3 -I build_lineage_panels.py GBFF_GZ OUTDIR [METAZOA_DEPTH] [FUNGI_DEPTH] [PER_GROUP]"""
import gzip
import json
import re
import sys
from collections import Counter, defaultdict

if __name__ == "__main__":
    gbff, out = sys.argv[1:3]
    DEPTH = {"metazoa": int(sys.argv[3]) if len(sys.argv) > 3 else 8,
             "fungi": int(sys.argv[4]) if len(sys.argv) > 4 else 7}
    PER_GROUP = int(sys.argv[5]) if len(sys.argv) > 5 else 1

ROMAN = {"i": "1", "ii": "2", "iii": "3", "1": "1", "2": "2", "3": "3"}


def norm(gene, product):
    g = (gene or "").lower().replace("-", "").replace("_", "").replace(" ", "")
    m = re.fullmatch(r"(?:cox|coi|co|cox)(i{1,3}|[123])", g)
    if m:
        return "cox" + ROMAN[m.group(1)]
    if g in ("cob", "cytb", "cyb"):
        return "cob"
    m = re.fullmatch(r"(?:nad|nd|nadh)([1-6]|4l)", g)
    if m:
        return "nad" + m.group(1)
    m = re.fullmatch(r"(?:atp|atpase)([689])", g)
    if m:
        return "atp" + m.group(1)
    if g in ("rps3", "var1"):
        return "rps3"
    p = (product or "").lower()
    m = re.search(r"cytochrome c oxidase subunit (iii|ii|i|[123])\b", p)
    if m:
        return "cox" + ROMAN[m.group(1)]
    if re.search(r"\bcytochrome b\b", p) and "b5" not in p and "c1" not in p:
        return "cob"
    m = re.search(r"nadh dehydrogenase subunit (4l|[1-6])\b", p)
    if m:
        return "nad" + m.group(1)
    m = re.search(r"atp synthase (?:f0 )?subunit ([689])\b", p)
    if m:
        return "atp" + m.group(1)
    if re.search(r"ribosomal protein s3\b", p):
        return "rps3"
    return None


def records(path):
    rec = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("//"):
                yield rec
                rec = []
            else:
                rec.append(line.rstrip("\n"))


def parse(lines):
    acc, organism, lineage, cds = None, None, [], []
    i, n = 0, len(lines)
    while i < n:
        l = lines[i]
        if l.startswith("VERSION"):
            acc = l.split()[1]
        elif l.startswith("  ORGANISM"):
            organism = l[12:].strip()
            j = i + 1
            tax = []
            while j < n and lines[j].startswith(" " * 12) and not lines[j].startswith("REFERENCE"):
                tax.append(lines[j].strip())
                j += 1
            lineage = [t.strip().rstrip(".") for t in " ".join(tax).split(";") if t.strip()]
        elif l.startswith("     CDS "):
            quals, j = {}, i + 1
            key = None
            while j < n and lines[j].startswith(" " * 21):
                t = lines[j][21:]
                if t.startswith("/"):
                    key, _, v = t[1:].partition("=")
                    quals[key] = v.strip('"')
                elif key:
                    quals[key] = quals[key] + (t.strip('"') if key == "translation" else " " + t.strip('"'))
                j += 1
            cds.append(quals)
            i = j - 1
        i += 1
    return acc, organism, lineage, cds


if __name__ == "__main__":
    groups = {"metazoa": defaultdict(list), "fungi": defaultdict(list)}
    for lines in records(gbff):
        acc, organism, lineage, cds = parse(lines)
        if not acc or not lineage:
            continue
        kind = "metazoa" if "Metazoa" in lineage else "fungi" if "Fungi" in lineage else None
        if kind is None or "mitochondrion" not in " ".join(lines[:40]).lower() and "mitochondri" not in organism.lower():
            pass
        if kind is None:
            continue
        prots = {}
        for q in cds:
            gname = norm(q.get("gene"), q.get("product"))
            aa = q.get("translation", "").replace(" ", "")
            if gname and len(aa) >= 40 and gname not in prots:
                prots[gname] = aa
        if len(prots) < 5:
            continue
        key = tuple(lineage[:DEPTH[kind]])
        groups[kind][key].append((len(prots), acc, organism, lineage, prots))

    summary = {}
    for kind, g in groups.items():
        chosen = []
        for key, members in g.items():
            members.sort(key=lambda x: (-x[0], x[1]))
            chosen.extend(members[:PER_GROUP])
        genes = Counter()
        phyla = Counter()
        with open(f"{out}/{kind}_mito_proteins.faa", "w") as fh:
            for _, acc, organism, lineage, prots in sorted(chosen, key=lambda x: x[1]):
                sp = re.sub(r"[^A-Za-z0-9]+", "_", organism).strip("_")
                phyla[lineage[2] if len(lineage) > 2 else "?"] += 1
                for gname, aa in sorted(prots.items()):
                    genes[gname] += 1
                    fh.write(f">mito|{gname}|{sp}|{acc}\n")
                    for k in range(0, len(aa), 80):
                        fh.write(aa[k:k + 80] + "\n")
        with open(f"{out}/{kind}_mito_proteins.tsv", "w") as fh:
            fh.write("accession\torganism\tgenes\tlineage\n")
            for n_, acc, organism, lineage, prots in sorted(chosen, key=lambda x: x[1]):
                fh.write(f"{acc}\t{organism}\t{n_}\t{'; '.join(lineage)}\n")
        summary[kind] = {"genomes_with_5_core_genes": sum(len(m) for m in g.values()), "groups": len(g),
                         "chosen": len(chosen), "proteins": sum(genes.values()), "genes": dict(sorted(genes.items())),
                         "top_phyla": phyla.most_common(15), "depth": DEPTH[kind], "per_group": PER_GROUP}
    json.dump(summary, open(f"{out}/lineage_panels.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))
