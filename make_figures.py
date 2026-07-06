"""
Generate two figures from Allen Brain Atlas striatal gene expression data:
1. Heatmap — expression energy of dopamine pathway genes across striatal subregions
2. Gradient plot — DRD1/DRD2 ratio and DAT expression across the D→V axis
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path

OUT = Path(__file__).parent

# ── Load data ─────────────────────────────────────────────────────────────────
data = {}  # data[gene][struct] = expression_energy
with open(OUT / "atlas_expression.tsv") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        g = row["gene"]
        s = row["structure"]
        v = row["expression_energy"]
        if v and v != "None":
            data.setdefault(g, {})[s] = float(v)

GENES_ORDER = ["Drd1", "Drd2", "Drd3", "Adora2a", "Ppp1r1b",
               "Penk", "Pdyn", "Th", "Slc6a3", "Slc18a2"]
GENES_NICE  = ["Drd1 (D1R)", "Drd2 (D2R)", "Drd3 (D3R)", "Adora2a (A2aR)", "Ppp1r1b (DARPP-32)",
               "Penk (Enkephalin)", "Pdyn (Dynorphin)", "Th (TH)", "Slc6a3 (DAT)", "Slc18a2 (VMAT2)"]

# Dorsal → Ventral ordering (matches Berke lab gradient axis)
STRUCTS = ["CP", "ACB", "OT"]
STRUCT_LABELS = ["Caudoputamen\n(dorsal)", "Nucleus accumbens\n(ventral)", "Olfactory tubercle\n(most ventral)"]

# Build matrix (genes × structures)
matrix = np.zeros((len(GENES_ORDER), len(STRUCTS)))
for i, g in enumerate(GENES_ORDER):
    for j, s in enumerate(STRUCTS):
        matrix[i, j] = data.get(g, {}).get(s, 0)

# ── Figure 1: Heatmap ─────────────────────────────────────────────────────────
# Row-normalize so patterns are visible across genes with very different scales
row_max = matrix.max(axis=1, keepdims=True)
row_max[row_max == 0] = 1
mat_norm = matrix / row_max

fig, ax = plt.subplots(figsize=(5.5, 4.2))
im = ax.imshow(mat_norm, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)

ax.set_xticks(range(len(STRUCTS)))
ax.set_xticklabels(STRUCT_LABELS, fontsize=8)
ax.set_yticks(range(len(GENES_ORDER)))
ax.set_yticklabels(GENES_NICE, fontsize=8)
ax.xaxis.set_ticks_position("top")
ax.xaxis.set_label_position("top")

# Annotate cells with raw values (log-scale for readability)
for i in range(len(GENES_ORDER)):
    for j in range(len(STRUCTS)):
        v = matrix[i, j]
        txt = f"{v:.2f}" if v < 10 else f"{v:.1f}"
        col = "white" if mat_norm[i, j] > 0.65 else "#333"
        ax.text(j, i, txt, ha="center", va="center", fontsize=6.5, color=col)

cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
cb.set_label("Normalized expression energy\n(row max = 1)", fontsize=7)
cb.ax.tick_params(labelsize=7)

ax.set_title("Dopamine pathway gene expression\nacross striatal subregions (Allen Mouse Brain Atlas)",
             fontsize=8.5, fontweight="bold", pad=14)
ax.tick_params(axis="x", length=0)
ax.tick_params(axis="y", length=0)
fig.tight_layout()
fig.savefig(OUT / "heatmap.svg", format="svg")
fig.savefig(OUT / "heatmap_preview.png", format="png", dpi=140)
plt.close(fig)
print("Saved heatmap.svg")

# ── Figure 2: Gradient plot ───────────────────────────────────────────────────
# Show DRD1, DRD2, and DRD1/DRD2 ratio along the D→V axis
# This is the molecular correlate of Berke lab's reward time horizon gradient

fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.2))
x = np.arange(len(STRUCTS))

# Panel A: absolute expression of D1R and D2R
ax = axes[0]
drd1 = [data.get("Drd1", {}).get(s, 0) for s in STRUCTS]
drd2 = [data.get("Drd2", {}).get(s, 0) for s in STRUCTS]
penk = [data.get("Penk", {}).get(s, 0) for s in STRUCTS]
pdyn = [data.get("Pdyn", {}).get(s, 0) for s in STRUCTS]

ax.plot(x, drd1, "o-", color="#e07b39", linewidth=1.8, markersize=6, label="Drd1 (D1R, direct path)")
ax.plot(x, drd2, "s-", color="#5b8ecf", linewidth=1.8, markersize=6, label="Drd2 (D2R, indirect path)")
ax.plot(x, pdyn, "^--", color="#c0392b", linewidth=1.2, markersize=5, alpha=0.7, label="Pdyn (dynorphin)")
ax.plot(x, penk, "v--", color="#2980b9", linewidth=1.2, markersize=5, alpha=0.7, label="Penk (enkephalin)")

ax.set_xticks(x)
ax.set_xticklabels(["Caudoputamen\n(dorsal)", "Accumbens\n(ventral)", "Olf. tubercle\n(most ventral)"],
                   fontsize=7.5)
ax.set_ylabel("Expression energy (ISH)", fontsize=8)
ax.set_title("D1R vs D2R along\ndorsal→ventral axis", fontsize=8.5, fontweight="bold")
ax.legend(fontsize=6.8, frameon=True, framealpha=0.9)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(labelsize=8)

# Add dorsal-ventral arrow annotation
ax.annotate("", xy=(2, ax.get_ylim()[0] * 0.05 + ax.get_ylim()[1] * 0.95),
            xytext=(0, ax.get_ylim()[0] * 0.05 + ax.get_ylim()[1] * 0.95),
            arrowprops=dict(arrowstyle="->", color="#888", lw=1.2))

# Panel B: DRD1/DRD2 ratio (proxy for direct/indirect pathway balance)
ax2 = axes[1]
ratio = [d1 / d2 if d2 > 0 else 0 for d1, d2 in zip(drd1, drd2)]
colors_bar = ["#c0392b" if r == max(ratio) else "#5d8aa8" for r in ratio]
bars = ax2.bar(x, ratio, color=colors_bar, edgecolor="white", linewidth=0.6, width=0.55)
ax2.set_xticks(x)
ax2.set_xticklabels(["Caudoputamen\n(dorsal)", "Accumbens\n(ventral)", "Olf. tubercle\n(most ventral)"],
                    fontsize=7.5)
ax2.set_ylabel("Drd1 / Drd2 expression ratio", fontsize=8)
ax2.set_title("Direct/indirect pathway balance\nacross striatal subregions", fontsize=8.5, fontweight="bold")
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.tick_params(labelsize=8)

# Value labels on bars
for bar, v in zip(bars, ratio):
    ax2.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
             f"{v:.2f}", ha="center", va="bottom", fontsize=8)

fig.tight_layout(pad=1.5)
fig.savefig(OUT / "gradient.svg", format="svg")
fig.savefig(OUT / "gradient_preview.png", format="png", dpi=140)
plt.close(fig)
print("Saved gradient.svg")
