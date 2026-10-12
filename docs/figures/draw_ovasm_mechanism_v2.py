"""OVASM mechanism schematic, v2 (redesigned layout). Geometry is illustrative; no measured values.
Content follows ovasm-mechanism-integrated (public OVASM b529d67, OrganelleVerse integration 330ce767).
Layout: numbered serpentine flow, row 1 left-to-right (1-3), row 2 right-to-left (4-7). 183 x 124 mm."""
from pathlib import Path
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle, Circle, PathPatch, Arc
from matplotlib.path import Path as MPath

OUT = Path(__file__).resolve().parent
NAME = "ovasm-mechanism-v2"
INK, GREY, SOFT = "#24343D", "#6E7C85", "#C9D3D8"
TEAL, VIOLET, AMBER = "#137F88", "#74699A", "#CF9B48"
TINT, NUC = "#EEF5F6", "#9AA5AB"
mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 6.5, "svg.fonttype": "none", "pdf.fonttype": 42})
W, H = 183, 124
fig = plt.figure(figsize=(7.204724409448819, 4.881889763779528))  # 183 x 124 mm
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, W); ax.set_ylim(H, 0); ax.set_axis_off()

def T(x, y, s, size=6.2, color=INK, weight="normal", ha="left", va="center", style="normal", z=6):
    return ax.text(x, y, s, fontsize=size, color=color, fontweight=weight, ha=ha, va=va, linespacing=1.25, fontstyle=style, zorder=z)
def seg(pts, color, lw=1.6, dash=False, z=2):
    p = PathPatch(MPath(pts, [MPath.MOVETO] + [MPath.LINETO] * (len(pts) - 1)), fill=False, edgecolor=color, lw=lw,
                  capstyle="round", joinstyle="round", linestyle=(0, (2.5, 1.8)) if dash else "-", zorder=z)
    ax.add_patch(p)
def arrow(pts, color=INK, lw=0.8, dash=False):
    ax.add_patch(FancyArrowPatch(path=MPath(pts, [MPath.MOVETO] + [MPath.LINETO] * (len(pts) - 1)), arrowstyle="-|>",
                 mutation_scale=7, lw=lw, color=color, linestyle=(0, (2.5, 1.8)) if dash else "-", joinstyle="round", zorder=3))
def card(x, y, w, h, num, title, fill="white"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.6", fc=fill, ec=SOFT, lw=0.6, zorder=1))
    ax.add_patch(Circle((x + 4.2, y + 4.4), 2.3, fc=INK, ec="none", zorder=4))
    T(x + 4.2, y + 4.45, str(num), 6.4, "white", "bold", ha="center")
    T(x + 8.0, y + 4.4, title, 7.2, INK, "bold")
def rep(x, y, w=4.6, h=3.0):
    ax.add_patch(Rectangle((x - w / 2, y - h / 2), w, h, fc=AMBER, ec="none", zorder=3))
def node(x, y, col, r=1.5):
    ax.add_patch(Circle((x, y), r, fc=col, ec="white", lw=0.5, zorder=4))
def chip(x, y, w, h, s, fc, tc="white", size=5.6, ec="none", weight="bold"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.0", fc=fc, ec=ec, lw=0.6, zorder=2))
    T(x + w / 2, y + h / 2, s, size, tc, weight, ha="center")

# ------------------------------------------------------------------ row 1
R1, H1 = 3, 50
# 1 reads
card(3, R1, 38, H1, 1, "Sequencing reads")
for k, (x0, ln, col) in enumerate([(7, 28, TEAL), (10, 27, VIOLET), (6, 30, TEAL), (9, 26, VIOLET)]):
    y = R1 + 13 + k * 4.2
    seg([(x0, y), (x0 + ln, y)], col, 1.5); seg([(x0 + 11, y), (x0 + 16, y)], AMBER, 1.5)
T(7, R1 + 34, "HiFi · ONT / CLR · Illumina", 6.2)
T(7, R1 + 39.5, "whole-genome or\norganelle-enriched input", 5.8, GREY)

# 2 read selection and graph
card(46, R1, 56, H1, 2, "Read selection and graph")
lanes = [("HiFi", "native closed-syncmer graph\n(Syncasm, from Oatk) or portable graph", TEAL),
         ("ONT / CLR", "self-correction,\nthen graph construction", VIOLET),
         ("Illumina", "dense k-mer graph", GREY)]
for k, (lab, desc, col) in enumerate(lanes):
    y = R1 + 11 + k * 9.6
    chip(49, y, 14.5, 6.2, lab, col, size=5.6)
    T(65.5, y + 3.1, desc, 5.7, INK)
T(49, R1 + 43.5, "reads recruited by organelle seeds or discovery;\nmultichrom HiFi can use the full WGS input", 5.6, GREY)

# 3 repeat resolution (dominant mechanism)
card(107, R1, 73, H1, 3, "Read-supported repeat resolution", fill=TINT)
# collapsed repeat with two read threads
cx, cy = 125, R1 + 23
L1, L2, Rr1, Rr2 = (111, cy - 7), (111, cy + 7), (139, cy - 7), (139, cy + 7)
for a in (L1, L2):
    seg([a, (cx - 2.3, cy)], SOFT, 0.9)
for b in (Rr1, Rr2):
    seg([(cx + 2.3, cy), b], SOFT, 0.9)
seg([L1, (cx - 2.3, cy - 0.6), (cx + 2.3, cy - 0.6), Rr1], TEAL, 1.3)
seg([L2, (cx - 2.3, cy + 0.6), (cx + 2.3, cy + 0.6), Rr2], VIOLET, 1.3)
rep(cx, cy, 4.8, 3.2)
for p_, col in ((L1, TEAL), (Rr1, TEAL), (L2, VIOLET), (Rr2, VIOLET)):
    node(*p_, col)
T(cx, R1 + 37, "collapsed repeat;\nreads span flank–repeat–flank", 5.6, GREY, ha="center")
arrow([(142.5, cy), (148.5, cy)])
# resolved paths
for y, col in ((cy - 5, TEAL), (cy + 5, VIOLET)):
    seg([(152, y), (176, y)], col, 1.5); node(152, y, col); node(176, y, col); rep(164, y, 4.6, 3.0)
T(164, R1 + 37, "repeat copies unzipped\nwhere reads support a pairing", 5.6, GREY, ha="center")
T(110.5, R1 + 45.5, "ambiguous or partial walks are kept as alternatives in the evidence", 5.6, INK)

# ------------------------------------------------------------------ row 2 (right to left)
R2, H2 = 64, 46
# 4 consensus (under 3, right)
card(143, R2, 37, H2, 4, "Consensus")
for k, (a, b) in enumerate([(147, 172), (150, 176), (148, 174)]):
    y = R2 + 12 + k * 3.6
    seg([(a, y), (b, y)], TEAL, 1.0); seg([(159, y), (164, y)], AMBER, 1.0)
seg([(147, R2 + 25), (176, R2 + 25)], TEAL, 2.0); seg([(159, R2 + 25), (164, R2 + 25)], AMBER, 2.0)
T(147, R2 + 28.3, "consensus", 5.2, TEAL)
T(146.5, R2 + 35, "original observations pooled\nat shared coordinates; each\nread position counted once", 5.6, INK)
T(146.5, R2 + 42.2, "native multichrom strategy", 5.4, GREY, style="italic")

# 5 component assessment
card(66, R2, 72, H2, 5, "Organelle component assessment")
# component: cycle with marker genes and coverage
ccx, ccy = 76, R2 + 20
ax.add_patch(Circle((ccx, ccy), 5.2, fc="none", ec=TEAL, lw=1.6, zorder=3))
for ang in (20, 140, 260):
    ax.add_patch(Arc((ccx, ccy), 10.4, 10.4, theta1=ang, theta2=ang + 28, color=AMBER, lw=2.4, zorder=4))
T(ccx, R2 + 30, "marker genes\n+ coverage", 5.5, GREY, ha="center")
# source assessment: original-read contexts vs nuclear flank / periodic fits
T(87, R2 + 11.5, "optional source assessment of original reads", 5.6, INK)
seg([(88, R2 + 17), (108, R2 + 17)], TEAL, 1.5); seg([(97, R2 + 17), (101, R2 + 17)], AMBER, 1.5)
T(110, R2 + 17, "organelle context", 5.4, TEAL)
seg([(88, R2 + 22), (95, R2 + 22)], NUC, 1.5); seg([(95, R2 + 22), (101, R2 + 22)], AMBER, 1.5); seg([(101, R2 + 22), (108, R2 + 22)], NUC, 1.5)
T(110, R2 + 22, "nuclear-flank fit", 5.4, GREY)
for k in range(4):
    seg([(88 + k * 5.2, R2 + 27), (91.5 + k * 5.2, R2 + 27)], AMBER, 1.5)
T(110, R2 + 27, "periodic fit", 5.4, GREY)
# decision outcomes
chip(70, R2 + 35, 30, 8.4, "", "#ECEFF1", ec=SOFT)
T(85, R2 + 37.4, "source-exclusive evidence,", 5.4, INK, "bold", ha="center")
T(85, R2 + 40.6, "eligible singleton cycle → exclude", 5.4, INK, ha="center")
chip(103, R2 + 35, 32, 8.4, "", "#E1EFF0", ec=TEAL)
T(119, R2 + 37.4, "compatible or unresolved", 5.4, TEAL, "bold", ha="center")
T(119, R2 + 40.6, "→ retain with evidence status", 5.4, TEAL, ha="center")

# 6 export
card(32, R2, 30, H2, 6, "Export")
for k, (lab, col) in enumerate((("mitochondrion", TEAL), ("plastid", VIOLET))):
    y = R2 + 13 + k * 8
    seg([(36, y), (58, y)], col, 1.8); seg([(45, y), (49, y)], AMBER, 1.8)
    T(36, y + 3.2, lab, 5.4, col)
T(35.5, R2 + 32, "overlaps at base\ncoordinates trimmed", 5.6, INK)
T(35.5, R2 + 40.5, "GFA · FASTA ·\nevidence report", 5.6, GREY)

# 7 evaluation (optional)
card(3, R2, 25, H2, 7, "Scoring")
T(15.5, R2 + 9.6, "optional, after assembly", 5.2, GREY, ha="center", style="italic")
chip(5.5, R2 + 13, 20, 6.4, "same-sample\nreference", "white", tc=GREY, ec=GREY, size=5.2, weight="normal")
T(5.5, R2 + 24.5, "relaxed", 5.6, INK, "bold"); T(5.5, R2 + 26.8, "complete, accurate\nper-reference recovery", 5.2, GREY, va="top")
T(5.5, R2 + 35.2, "strict", 5.6, INK, "bold"); T(5.5, R2 + 37.5, "+ specificity,\ncorrespondence,\nlength limits", 5.2, GREY, va="top")

# ------------------------------------------------------------------ flow arrows
arrow([(41.3, R1 + 25), (45.7, R1 + 25)])
arrow([(102.3, R1 + 25), (106.7, R1 + 25)])
arrow([(161.5, R1 + H1 + 0.3), (161.5, R2 - 0.3)])                 # 3 -> 4
arrow([(142.7, R2 + 22), (138.3, R2 + 22)])                        # 4 -> 5
arrow([(65.7, R2 + 22), (62.3, R2 + 22)])                          # 5 -> 6
arrow([(31.7, R2 + 22), (28.3, R2 + 22)], dash=True, color=GREY)   # 6 -> 7 (optional)
# evidence input: original reads feed source assessment
arrow([(22, R1 + H1 + 0.3), (22, R1 + H1 + 5.5), (98, R1 + H1 + 5.5), (98, R2 - 0.3)], color=GREY, dash=True)
T(60, R1 + H1 + 3.9, "original reads (evidence input)", 5.4, GREY, ha="center", va="bottom")

# ------------------------------------------------------------------ key + caveat
yk = 120.5
seg([(3, yk), (9, yk)], TEAL, 1.6); seg([(10, yk), (16, yk)], VIOLET, 1.6); T(17.5, yk, "distinguishable sequence contexts", 5.6, GREY)
rep(57, yk, 4, 2.4); T(60, yk, "shared or repeated sequence", 5.6, GREY)
arrow([(91, yk), (98, yk)]); T(99.5, yk, "processing", 5.6, GREY)
arrow([(114, yk), (121, yk)], color=GREY, dash=True); T(122.5, yk, "evidence input or optional step", 5.6, GREY)
T(180, yk, "Graph cycles and FASTA records do not identify physical molecule number or topology.", 5.4, GREY, ha="right", va="center") if False else None
T(3, yk - 5.2, "Graph cycles and FASTA records do not identify physical molecule number or topology; geometry is schematic.", 5.4, GREY)

fig.savefig(OUT / "ovasm-mechanism-v2.svg")
fig.savefig(OUT / "ovasm-mechanism-v2.pdf")
fig.savefig(OUT / "ovasm-mechanism-v2.png", dpi=600)
print("saved", OUT / NAME)
