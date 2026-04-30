"""
Rep2 Step 2: 整合 cellbin / bin20 / nucleus 的折线图
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REP2_DIR = "/data/scratch/projects/punim1236/AAA_chen_here/AAA克隆_XCR/Rep2"
OUT_DIR  = os.path.join(REP2_DIR, "Integrated_line_plot")
os.makedirs(OUT_DIR, exist_ok=True)

datasets = {
    "Bin20":              {"line_csv": os.path.join(REP2_DIR, "bin20_BCR_QC",  "B_proportion_by_coords_line.csv"),
                           "n_col": "n_bins_matched", "unit": "bins",
                           "color": "#F8766D", "offset": -0.080},
    "Cell bin":           {"line_csv": os.path.join(REP2_DIR, "cellbin_BCR_QC","B_proportion_by_coords_line.csv"),
                           "n_col": "n_cells_matched", "unit": "cells",
                           "color": "#00BFC4", "offset": -0.060},
    "Nucleus-strict bin": {"line_csv": os.path.join(REP2_DIR, "nucleus_BCR_QC","B_proportion_by_coords_line.csv"),
                           "n_col": "n_cells_matched", "unit": "cells",
                           "color": "#C77CFF", "offset": +0.090},
}

X_MAX, MIN_N = 15, 5

def set_black(fig, ax):
    fig.patch.set_facecolor("black"); ax.set_facecolor("black")
    ax.tick_params(colors="white", labelsize=12)
    for el in (ax.xaxis.label, ax.yaxis.label, ax.title): el.set_color("white")
    for s in ax.spines.values(): s.set_color("white")
    ax.grid(True, axis="y", ls="--", alpha=0.25, color="white")

fig, ax = plt.subplots(figsize=(15, 9))
set_black(fig, ax)

for name, info in datasets.items():
    if not os.path.exists(info["line_csv"]):
        print(f"  Skip {name}: {info['line_csv']} not found"); continue
    line = pd.read_csv(info["line_csv"]).copy()
    if info["n_col"] not in line.columns:
        n_col = next(c for c in line.columns if c.startswith("n_") and "matched" in c)
    else:
        n_col = info["n_col"]
    line["n_matched"] = line[n_col]
    plot_df = line[line["n_coords"] <= X_MAX].copy()
    enough = plot_df["n_matched"] >= MIN_N

    ax.plot(plot_df["n_coords"], plot_df["B_proportion"],
            ls="--", lw=1.2, alpha=0.35, color=info["color"])
    sub_ok = plot_df[enough]
    ax.plot(sub_ok["n_coords"], sub_ok["B_proportion"],
            ls="-", lw=3.0, marker="o", markersize=9,
            markerfacecolor=info["color"],
            markeredgecolor="white", markeredgewidth=0.8,
            color=info["color"], label=name, zorder=4)
    for _, r in plot_df.iterrows():
        if pd.notna(r["B_proportion"]):
            ax.text(r["n_coords"], r["B_proportion"]+info["offset"],
                    f"{r['B_proportion']*100:.0f}%\n{int(r['n_matched']):,} {info['unit']}",
                    ha="center", va="center", fontsize=8.5,
                    fontweight="bold", color=info["color"],
                    zorder=6, linespacing=1.1)

ax.set_xlim(-0.5, X_MAX+0.5); ax.set_ylim(0, 1.18)
ax.set_xticks(range(0, X_MAX+1))
ax.set_yticks(np.arange(0,1.01,0.2))
ax.set_yticklabels([f"{int(v*100)}%" for v in np.arange(0,1.01,0.2)])
ax.set_xlabel("BCR coordinates per bin", fontsize=14)
ax.set_ylabel("B cell proportion", fontsize=14)
ax.set_title("Rep2: More BCR signal → more B cells", fontsize=18, pad=14)
leg = ax.legend(frameon=False, fontsize=13, loc="lower right")
for t in leg.get_texts(): t.set_color("white")
fig.tight_layout()

for ext in ("png","pdf"):
    p = os.path.join(OUT_DIR, f"Rep2_integrated_B_proportion.{ext}")
    fig.savefig(p, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    print("saved:", p)
plt.close(fig)