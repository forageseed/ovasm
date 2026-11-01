"""ovasm vs published Oatk assemblies on the DToL HiFi panel, drawn with OrganelleVerse's journal style.

usage: python plot_oatk_compare.py COMPARE_TSV MANIFEST_TSV OUT_PREFIX
"""
import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from organelleverse.visualization.style import FONT_PT, add_panel_label, publication_style, width_inches  # noqa: E402

cmp_tsv, manifest_tsv, out = sys.argv[1:4]
rows = list(csv.DictReader(open(cmp_tsv), delimiter="\t"))
clade = {r["sample"]: r["clade"] for r in csv.DictReader(open(manifest_tsv), delimiter="\t")}
COLOR = {"mitochondrion": "#2a78d6", "plastid": "#eb6834"}
NAME = {"mitochondrion": "Mitochondrion", "plastid": "Plastid"}
LABEL = {"Viscum_album", "Ballota_nigra", "Lathyrus_linifolius", "Linaria_vulgaris", "Medicago_arabica"}


def species(sample):
    return " ".join(sample.split("_")[:2])


with publication_style():
    fig, axes = plt.subplots(1, 3, figsize=(width_inches("double"), 62 / 25.4),
                             gridspec_kw={"wspace": 0.42, "left": 0.07, "right": 0.985, "bottom": 0.19, "top": 0.9})

    # a: assembled length, ovasm vs Oatk
    ax = axes[0]
    for org in ("plastid", "mitochondrion"):
        pts = [(int(r[f"{org}_oatk_bp"]) / 1e3, int(r[f"{org}_ovasm_bp"]) / 1e3, r["sample"])
               for r in rows if r[f"{org}_oatk_covered"] != ""]
        ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=9, color=COLOR[org], linewidths=0.4,
                   edgecolors="white", label=f"{NAME[org]} (n = {len(pts)})", zorder=3)
        for x, y, s in pts:
            if "_".join(s.split("_")[:2]) in LABEL and y < 0.9 * x:
                left = s.startswith("Medicago")
                ax.annotate(species(s), (x, y), xytext=(-3 if left else 3, -1), textcoords="offset points",
                            ha="right" if left else "left", va="center",
                            fontsize=FONT_PT["min"] + 0.5, style="italic", color="#444444")
    lim = (25, 2000)
    ax.plot(lim, lim, color="#999999", lw=0.6, ls="--", zorder=1)
    ax.set(xscale="log", yscale="log", xlim=lim, ylim=lim,
           xlabel="Oatk assembly (kb)", ylabel="ovasm assembly (kb)")
    ax.legend(loc="upper left", handletextpad=0.2, borderaxespad=0.2)
    add_panel_label(ax, "a")

    # b: how much of each Oatk assembly ovasm recovers, in bins
    ax = axes[1]
    bins = [("≥ 99", 99, 101), ("95–99", 95, 99), ("50–95", 50, 95), ("< 50", -1, 50)]
    width = 0.38
    for k, org in enumerate(("mitochondrion", "plastid")):
        vals = [float(r[f"{org}_oatk_covered"]) * 100 for r in rows if r[f"{org}_oatk_covered"] != ""]
        counts = [sum(lo <= v < hi for v in vals) for _, lo, hi in bins]
        xs = [i + (k - 0.5) * (width + 0.03) for i in range(len(bins))]
        ax.bar(xs, counts, width=width, color=COLOR[org], label=NAME[org], zorder=3)
        for x, c in zip(xs, counts):
            ax.text(x, c + 1.5, str(c), ha="center", va="bottom", fontsize=FONT_PT["tick"])
    ax.set_xticks(range(len(bins)), [b[0] for b in bins])
    ax.set(ylim=(0, 118), xlabel="Oatk bases covered by ovasm (%)", ylabel="Species")
    ax.legend(loc="upper right", handletextpad=0.3, borderaxespad=0.2)
    add_panel_label(ax, "b")

    # c: cost per species
    ax = axes[2]
    groups = {}
    for r in rows:
        c = clade.get(r["sample"]) or "Other"
        groups.setdefault("Bryophytes" if c == "Bryophytes" else "Vascular plants", []).append(r)
    for g, col in (("Vascular plants", "#1baf7a"), ("Bryophytes", "#e87ba4")):
        pts = groups.get(g, [])
        ax.scatter([float(r["wall_min"]) for r in pts], [float(r["peak_gb"]) for r in pts], s=9, color=col,
                   edgecolors="white", linewidths=0.4, label=f"{g} (n = {len(pts)})", zorder=3)
    from matplotlib.ticker import FuncFormatter
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set(xscale="log", xlabel="Wall time, 8 threads (min)", ylabel="Peak memory (GB)", ylim=(0, 8))
    ax.legend(loc="upper left", handletextpad=0.2, borderaxespad=0.2)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    add_panel_label(ax, "c")

    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    for ext in ("svg", "pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=300 if ext == "png" else None)
print("wrote", out)
