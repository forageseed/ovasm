"""Compare ovasm hifi-v1 atlas results with the published Oatk assemblies (Zenodo 15053552).

For every finished DToL sample: match its Zenodo ToLID by Genus[:3]+species[:4], align the
ovasm molecules against Oatk MITO.fasta / PLTD.fasta with minimap2 asm5, and report how much
of each side the other covers (union of aligned target intervals) and the matched identity.
usage: python3 -I compare_oatk.py ATLAS_DIR ZENODO_DATA MINIMAP2 OUT_TSV
"""
import json
import os
import subprocess
import sys

atlas, zdata, mm, out = sys.argv[1:5]
tolids = sorted(os.listdir(zdata))


def tolid_for(species):
    g, s = species.split("_")[:2]
    key = g[:3].capitalize() + s[:4].capitalize()
    hits = [t for t in tolids if key in t]
    return hits[0] if len(hits) == 1 else None


def fasta_lengths(path):
    n, lens = None, {}
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                n = line[1:].split()[0]
                lens[n] = 0
            elif n:
                lens[n] += len(line.strip())
    return lens


def covered(query, target):
    """Fraction of target bp covered by query alignments, and matched identity."""
    paf = subprocess.run([mm, "-x", "asm5", "-t", "2", "-N", "50", "-p", "0.5", target, query],
                         capture_output=True, text=True).stdout
    tlen = fasta_lengths(target)
    iv = {t: [] for t in tlen}
    match = alen = 0
    for row in paf.splitlines():
        c = row.split("\t")
        if int(c[11]) == 0 and int(c[10]) < 1000:
            continue
        iv[c[5]].append((int(c[7]), int(c[8])))
        match += int(c[9])
        alen += int(c[10])
    cov = 0
    for spans in iv.values():
        end = -1
        for a, b in sorted(spans):
            if b > end:
                cov += b - max(a, end)
                end = b
    total = sum(tlen.values())
    return (cov / total if total else 0.0), (match / alen if alen else 0.0)


rows = []
for sample in sorted(os.listdir(atlas)):
    st = os.path.join(atlas, sample, "status.json")
    if sample.endswith(".tmp") or not os.path.exists(st):
        continue
    status = json.load(open(st))
    if status.get("source") != "dtol" or status.get("exit_code") != 0:
        continue
    tol = tolid_for(status["species"])
    row = {"sample": sample, "tolid": tol or "", "wall_min": round(status["wall_seconds"] / 60, 1),
           "peak_gb": status["peak_rss_gb"]}
    for org, ref in (("mitochondrion", "MITO.fasta"), ("plastid", "PLTD.fasta")):
        ov = os.path.join(atlas, sample, "ovasm", org, "molecules.fasta")
        ok = os.path.join(zdata, tol, ref) if tol else ""
        ovl = sum(fasta_lengths(ov).values()) if os.path.exists(ov) else 0
        okl = sum(fasta_lengths(ok).values()) if ok and os.path.exists(ok) else 0
        row[f"{org}_ovasm_bp"] = ovl
        row[f"{org}_oatk_bp"] = okl
        if ovl and okl:
            c1, i1 = covered(ov, ok)
            c2, _ = covered(ok, ov)
            row[f"{org}_oatk_covered"] = round(c1, 4)
            row[f"{org}_ovasm_covered"] = round(c2, 4)
            row[f"{org}_identity"] = round(i1, 5)
        else:
            row[f"{org}_oatk_covered"] = row[f"{org}_ovasm_covered"] = row[f"{org}_identity"] = ""
    rows.append(row)
    print(sample, tol, row.get("mitochondrion_oatk_covered"), row.get("plastid_oatk_covered"), flush=True)

with open(out, "w") as f:
    keys = list(rows[0])
    f.write("\t".join(keys) + "\n")
    for r in rows:
        f.write("\t".join(str(r[k]) for k in keys) + "\n")
