"""
===========================================================================
BIN20 version BCR QC analysis (IGH / IGK / IGL)
===========================================================================
Analysis for the bin20 (fixed 20x20 square bin) repeated cell-bin version:
  1. Per-chain unique targetSeq frequency + spatial binning + dominant composition
  2. Per-chain unique coords frequency + spatial binning + dominant composition
  3. Dominant-B heatmap for IGH x IGK unique targetSeq
  4. Dominant-B heatmap for IGH x IGK unique coords

RCTD convention:
  - Use first_type directly
  - first_type in B_related is treated as B cell (dominant)

Output directory structure:
  BASE_OUT/
    bin20_IGH/
    bin20_IGK/
    bin20_IGL/
    bin20_sequence_heatmap/
    bin20_coordinates_heatmap/
===========================================================================
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from itertools import product
import anndata as ad


# =========================================================================
# User configuration
# =========================================================================
H5AD_PATH = "/data/scratch/projects/punim1236/spleen_stomics/processed/spleen_out/outs/visualization/C04278C6.bin20_1.0.h5ad"
META_PATH = "/data/scratch/projects/punim1236/AAA_chen_here/AAA克隆_XCR/C6_read1.fq.gz.meta.gz"
RCTD_PATH = "/data/scratch/projects/punim1236/AAA_chen_here/STOmics官方指南/RCTD_Clustering比较/GCsurrounding.csv"
BASE_OUT  = "/data/scratch/projects/punim1236/AAA_chen_here/AAA克隆_XCR/扫描高占比B/BIn20_扫描高占比"

CHAINS = ["IGH", "IGK", "IGL"]
SEQ_COL = "targetSequences"

B_RELATED = [
    "Activated", "Atypical", "B ASC", "B B1", "B follicular", "B MZ", "B1",
    "DZd", "DZp", "Follicular Naive", "GC", "Immature_Transitional",
    "LZ", "LZ_DZ_intermediate", "Memory", "MZ", "PB", "preMem", "preMem intermediate",
    # Common names from GCsurrounding.csv are also included here
    "B GC", "B naive",
]

# Heatmap axis ranges (set wide; will be clipped to actual data max at runtime)
SEQ_AXIS_MAX_IGH   = 8
SEQ_AXIS_MAX_IGK   = 10
COORD_AXIS_MAX_IGH = 6
COORD_AXIS_MAX_IGK = 12

# Gold-border highlight range
TARGET_B_LOW  = 0.80
TARGET_B_HIGH = 0.90

# Subfolders
DIRS = {
    "IGH":      os.path.join(BASE_OUT, "bin20_IGH"),
    "IGK":      os.path.join(BASE_OUT, "bin20_IGK"),
    "IGL":      os.path.join(BASE_OUT, "bin20_IGL"),
    "seq_hm":   os.path.join(BASE_OUT, "bin20_sequence_heatmap"),
    "coord_hm": os.path.join(BASE_OUT, "bin20_coordinates_heatmap"),
}
for d in DIRS.values():
    os.makedirs(d, exist_ok=True)


# =========================================================================
# Step 1: Load bin20 adata and auto-detect bin coordinate convention
# =========================================================================

def load_bin20(h5ad_path):
    print(f"Loading bin20 h5ad: {h5ad_path}")
    adata = ad.read_h5ad(h5ad_path)

    obs = adata.obs
    if {"x", "y"}.issubset(obs.columns):
        bx = obs["x"].to_numpy(dtype=float)
        by = obs["y"].to_numpy(dtype=float)
    else:
        sp = adata.obsm["spatial"]
        bx = sp[:, 0].astype(float)
        by = sp[:, 1].astype(float)

    # Auto-infer bin size: use the minimum positive difference
    # between adjacent unique values along x and y
    ux = np.sort(np.unique(bx))
    uy = np.sort(np.unique(by))
    dx = np.min(np.diff(ux)) if len(ux) > 1 else 20.0
    dy = np.min(np.diff(uy)) if len(uy) > 1 else 20.0
    bin_size = float(round((dx + dy) / 2))
    if bin_size <= 0:
        bin_size = 20.0
    print(f"  Inferred bin_size ~ {bin_size}")

    # Auto-detect: are the coordinates at the lower-left corner or the center of the bin?
    # If lower-left corner, coordinates are mostly integer multiples of bin_size.
    # If center, they are typically bin_size*k + bin_size/2.
    frac_x = (bx / bin_size) - np.floor(bx / bin_size)
    near_integer = np.mean((frac_x < 0.05) | (frac_x > 0.95))
    near_half    = np.mean(np.abs(frac_x - 0.5) < 0.05)
    if near_integer >= near_half:
        convention = "corner"
    else:
        convention = "center"
    print(f"  Bin coordinate convention: {convention} "
          f"(near_integer={near_integer:.2f}, near_half={near_half:.2f})")

    bin_df = pd.DataFrame({"bx": bx, "by": by})
    return bin_df, bin_size, convention


def point_to_bin(px, py, bin_size, convention):
    """Map point coordinates to the reference coordinates of their owning bin
    (consistent with the obs convention)."""
    if convention == "corner":
        # Lower-left corner: floor to the bin_size grid
        rx = np.floor(px / bin_size) * bin_size
        ry = np.floor(py / bin_size) * bin_size
    else:
        # Center
        rx = np.floor(px / bin_size) * bin_size + bin_size / 2.0
        ry = np.floor(py / bin_size) * bin_size + bin_size / 2.0
    return rx, ry


# =========================================================================
# Step 2: Compute unique targetSeq / unique coords per bin for each chain
# =========================================================================

def build_chain_bin_stats(df, chain, bin_df, bin_size, convention):
    """
    Returns a DataFrame with columns:
      bx, by, n_unique_{chain}_targetSeq, n_unique_{chain}_coords
    Only records that fall within valid bin20 bins are kept.
    """
    sub = df[df["topChains"].astype(str).str.strip() == chain].copy()
    sub = sub.dropna(subset=["x", "y"]).copy()

    # Point -> bin
    rx, ry = point_to_bin(
        sub["x"].to_numpy(dtype=float),
        sub["y"].to_numpy(dtype=float),
        bin_size, convention
    )
    sub["bx"] = rx
    sub["by"] = ry

    # Inner merge keeps only points that fall within adata bins
    valid_bins = bin_df[["bx", "by"]].drop_duplicates()
    sub = pd.merge(sub, valid_bins, on=["bx", "by"], how="inner")

    # --- Unique targetSeq per bin ---
    seq = sub.dropna(subset=[SEQ_COL]).copy()
    seq[SEQ_COL] = seq[SEQ_COL].astype(str).str.strip()
    seq = seq[~seq[SEQ_COL].isin(["", "nan", "None", "none", "NA"])]
    seq_unique = seq[["bx", "by", SEQ_COL]].drop_duplicates()
    n_seq = (
        seq_unique.groupby(["bx", "by"]).size()
        .reset_index(name=f"n_unique_{chain}_targetSeq")
    )

    # --- Unique coords per bin ---
    coord_unique = sub[["bx", "by", "x", "y"]].drop_duplicates()
    n_coord = (
        coord_unique.groupby(["bx", "by"]).size()
        .reset_index(name=f"n_unique_{chain}_coords")
    )

    # Merge onto all bins (bins without records for this chain are filled with 0)
    stats = bin_df[["bx", "by"]].drop_duplicates().copy()
    stats = pd.merge(stats, n_seq,   on=["bx", "by"], how="left")
    stats = pd.merge(stats, n_coord, on=["bx", "by"], how="left")
    stats[f"n_unique_{chain}_targetSeq"] = stats[f"n_unique_{chain}_targetSeq"].fillna(0).astype(int)
    stats[f"n_unique_{chain}_coords"]    = stats[f"n_unique_{chain}_coords"].fillna(0).astype(int)
    return stats


# =========================================================================
# Step 3: RCTD (first_type) -> dominant_is_B flag
# =========================================================================

def prepare_rctd(rctd_path, b_related, bin_size, convention):
    print(f"Loading RCTD: {rctd_path}")
    rctd = pd.read_csv(rctd_path)
    # Map RCTD (x, y) to bin reference coordinates
    rx, ry = point_to_bin(
        rctd["x"].to_numpy(dtype=float),
        rctd["y"].to_numpy(dtype=float),
        bin_size, convention
    )
    rctd["bx"] = rx
    rctd["by"] = ry
    rctd["is_B_dominant"] = rctd["first_type"].isin(b_related)
    print(f"  RCTD rows: {len(rctd)},  B-dominant: {rctd['is_B_dominant'].sum()}")
    return rctd[["bx", "by", "first_type", "is_B_dominant"]]


# =========================================================================
# Step 4: Per-chain analysis (frequency + spatial binning + dominant composition)
# =========================================================================

def plot_chain_frequency(stats, chain, value_col, xlabel_kind, outdir):
    freq = (
        stats[value_col].value_counts().sort_index().reset_index()
    )
    freq.columns = [value_col, "frequency"]
    freq.to_csv(os.path.join(outdir, f"{chain}_{xlabel_kind}_frequency.csv"), index=False)

    plot_freq = freq[freq[value_col] <= 30].copy()
    total = plot_freq["frequency"].sum()

    plt.figure(figsize=(12, 10))
    plt.barh(plot_freq[value_col], plot_freq["frequency"], height=0.8)
    plt.ylabel(f"Number of unique {chain} {xlabel_kind} per bin20", fontsize=13)
    plt.xlabel("Number of bins", fontsize=13)
    plt.title(f"Frequency distribution of unique {chain} {xlabel_kind} per bin20", fontsize=15)
    ax = plt.gca()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.yticks(plot_freq[value_col])
    for y, x in zip(plot_freq[value_col], plot_freq["frequency"]):
        pct = x / total * 100 if total > 0 else 0
        plt.text(x, y, f"  {x} ({pct:.2f}%)", ha="left", va="center", fontsize=12)
    plt.tight_layout()
    png = os.path.join(outdir, f"{chain}_{xlabel_kind}_frequency.png")
    plt.savefig(png, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"    -> {png}")


def plot_chain_spatial_bins(stats, chain, value_col, xlabel_kind, outdir):
    df_plot = stats[["bx", "by", value_col]].copy()
    df_plot["bin_cat"] = df_plot[value_col].clip(upper=4)

    color_map = {0: "lightgrey", 1: "#4DBBD5", 2: "#00A087",
                 3: "#E64B35",   4: "#3C5488"}
    label_map = {0: "0", 1: "1", 2: "2", 3: "3", 4: ">=4"}

    plt.figure(figsize=(10, 10))
    for k in [0, 1, 2, 3, 4]:
        sub = df_plot[df_plot["bin_cat"] == k]
        plt.scatter(sub["bx"], sub["by"], s=3, c=color_map[k],
                    label=f"{label_map[k]} (n={len(sub)})",
                    linewidths=0, alpha=0.9)
    plt.gca().set_aspect("equal", "box")
    plt.gca().invert_yaxis()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"Spatial distribution of {chain} {xlabel_kind} count per bin20")
    plt.legend(frameon=False, markerscale=3)
    plt.tight_layout()
    png = os.path.join(outdir, f"{chain}_{xlabel_kind}_spatial_bins.png")
    plt.savefig(png, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"    -> {png}")

def plot_bin_dominant_composition(stats, rctd, chain, value_col,
                                  xlabel_kind, outdir):
    merged = pd.merge(
        stats[["bx", "by", value_col]],
        rctd[["bx", "by", "first_type"]],
        on=["bx", "by"], how="left"
    )
    matched = merged[merged["first_type"].notna()].copy()
    matched["bin_cat"] = matched[value_col].clip(upper=4)
    bin_label_map = {0: "0", 1: "1", 2: "2", 3: "3", 4: ">=4"}
    matched["bin_cat_label"] = matched["bin_cat"].map(bin_label_map)

    # === New: merge all cell types starting with "B" into "B cells" ===
    def _merge_b_types(ct):
        s = str(ct).strip()
        # Starts with "B " (B GC, B naive, B MZ, B ASC, B B1, B follicular...)
        # or equals "B" / "B1" on its own
        if s.startswith("B ") or s == "B" or s == "B1":
            return "B cells"
        return s
    matched["cell_type_merged"] = matched["first_type"].apply(_merge_b_types)
    # ===============================================

    all_comp = []
    for b in ["0", "1", "2", "3", ">=4"]:
        sub = matched[matched["bin_cat_label"] == b]
        if len(sub) == 0:
            continue
        count_df = (sub["cell_type_merged"].value_counts()   # <-- use merged types
                    .rename_axis("cell_type").reset_index(name="count"))
        count_df["fraction"] = count_df["count"] / count_df["count"].sum()
        count_df["bin_label"] = b
        all_comp.append(count_df)

        # Draw individual barplot
        plot_sub = count_df.copy()
        main_df = plot_sub[plot_sub["fraction"] > 0.01].copy()
        others_fraction = plot_sub.loc[plot_sub["fraction"] <= 0.01, "fraction"].sum()
        others_count = plot_sub.loc[plot_sub["fraction"] <= 0.01, "count"].sum()
        if others_fraction > 0:
            main_df = pd.concat([main_df, pd.DataFrame({
                "cell_type": ["Others"], "count": [others_count],
                "fraction": [others_fraction], "bin_label": [b]
            })], ignore_index=True)
        main_df = main_df.sort_values("fraction", ascending=True)

        plt.figure(figsize=(12, 10))
        # Accent color for "B cells", default blue for others
        colors = ["#E64B35" if ct == "B cells" else "#4C72B0"
                  for ct in main_df["cell_type"]]
        plt.barh(main_df["cell_type"], main_df["fraction"],
                 height=0.8, color=colors)
        plt.ylabel("Dominant cell type (first_type, B-types merged)", fontsize=13)
        plt.xlabel("Fraction", fontsize=13)
        plt.title(
            f"Dominant cell-type distribution in bin20 with {b} "
            f"unique {chain} {xlabel_kind}",
            fontsize=14
        )
        ax = plt.gca()
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for y, x, c in zip(main_df["cell_type"], main_df["fraction"], main_df["count"]):
            plt.text(x, y, f"  {c} ({x*100:.2f}%)",
                     ha="left", va="center", fontsize=11)
        plt.tight_layout()
        safe_b = b.replace(">=", "ge")
        png = os.path.join(outdir,
                           f"{chain}_{xlabel_kind}_bin{safe_b}_dominant_composition.png")
        plt.savefig(png, dpi=200, bbox_inches="tight")
        plt.close()

    if all_comp:
        out_csv = os.path.join(
            outdir, f"{chain}_{xlabel_kind}_dominant_composition_all.csv"
        )
        pd.concat(all_comp, ignore_index=True).to_csv(out_csv, index=False)
        print(f"    -> {out_csv}")


def run_chain_analysis(chain, stats, rctd, outdir):
    print(f"\n--- {chain} single-chain analysis ---")
    # targetSeq
    plot_chain_frequency(stats, chain, f"n_unique_{chain}_targetSeq",
                         "targetSeq", outdir)
    plot_chain_spatial_bins(stats, chain, f"n_unique_{chain}_targetSeq",
                            "targetSeq", outdir)
    plot_bin_dominant_composition(stats, rctd, chain,
                                  f"n_unique_{chain}_targetSeq",
                                  "targetSeq", outdir)
    # coords
    plot_chain_frequency(stats, chain, f"n_unique_{chain}_coords",
                         "coords", outdir)
    plot_chain_spatial_bins(stats, chain, f"n_unique_{chain}_coords",
                            "coords", outdir)
    plot_bin_dominant_composition(stats, rctd, chain,
                                  f"n_unique_{chain}_coords",
                                  "coords", outdir)

    out_csv = os.path.join(outdir, f"{chain}_bin20_stats.csv")
    stats.to_csv(out_csv, index=False)
    print(f"    -> {out_csv}")


# =========================================================================
# Step 5: IGH x IGK heatmap (shared across both conventions)
# =========================================================================

def compute_igh_igk_heatmap(stats_igh, stats_igk, rctd,
                            col_igh, col_igk,
                            igh_max, igk_max):
    joint = pd.merge(
        stats_igh[["bx", "by", col_igh]],
        stats_igk[["bx", "by", col_igk]],
        on=["bx", "by"], how="outer"
    )
    joint[col_igh] = joint[col_igh].fillna(0).astype(int)
    joint[col_igk] = joint[col_igk].fillna(0).astype(int)

    joint = pd.merge(
        joint, rctd[["bx", "by", "first_type", "is_B_dominant"]],
        on=["bx", "by"], how="left"
    )

    results = []
    for igh_n, igk_n in product(range(igh_max + 1), range(igk_max + 1)):
        sub = joint[(joint[col_igh] == igh_n) & (joint[col_igk] == igk_n)]
        n_total = len(sub)
        matched = sub[sub["first_type"].notna()]
        n_matched = len(matched)
        if n_matched == 0:
            b_frac = np.nan
            n_B = 0
        else:
            n_B = int(matched["is_B_dominant"].sum())
            b_frac = n_B / n_matched
        results.append({
            "IGH_n": igh_n, "IGK_n": igk_n,
            "n_bins": n_total, "n_matched": n_matched,
            "n_dominant_B": n_B, "B_fraction": b_frac,
        })
    return pd.DataFrame(results)


def plot_igh_igk_heatmap(sweep_df, kind, outdir, overall_b_frac):
    """kind: 'targetSeq' or 'coords'"""
    pivot_b = sweep_df.pivot(index="IGH_n", columns="IGK_n",
                             values="B_fraction").sort_index()
    pivot_matched = sweep_df.pivot(index="IGH_n", columns="IGK_n",
                                   values="n_matched").sort_index()
    pivot_domB = sweep_df.pivot(index="IGH_n", columns="IGK_n",
                                values="n_dominant_B").sort_index()

    fig, ax = plt.subplots(figsize=(16, 9))
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "bcell", ["#D9D9D9", "#4DBBD5", "#00A087", "#E64B35"], N=256
    )
    norm = mcolors.TwoSlopeNorm(vmin=0.3, vcenter=0.8, vmax=1.0)
    im = ax.imshow(pivot_b.values, cmap=cmap, norm=norm,
                   aspect="auto", origin="lower")

    for i in range(pivot_b.shape[0]):
        for j in range(pivot_b.shape[1]):
            b_val = pivot_b.values[i, j]
            dom_b = pivot_domB.values[i, j]
            matched = pivot_matched.values[i, j]
            if np.isnan(b_val) or matched == 0:
                ax.text(j, i, "-", ha="center", va="center",
                        fontsize=8, color="gray")
                continue
            color = "white" if b_val > 0.75 else "black"
            weight = "bold" if TARGET_B_LOW <= b_val <= TARGET_B_HIGH else "normal"
            ax.text(j, i, f"{b_val:.1%}\n{int(dom_b)}/{int(matched)}",
                    ha="center", va="center", fontsize=7,
                    color=color, fontweight=weight)

    ax.set_xticks(range(len(pivot_b.columns)))
    ax.set_xticklabels(pivot_b.columns.astype(int))
    ax.set_yticks(range(len(pivot_b.index)))
    ax.set_yticklabels(pivot_b.index.astype(int))
    ax.set_xlabel(f"Unique IGK {kind} per bin20", fontsize=12)
    ax.set_ylabel(f"Unique IGH {kind} per bin20", fontsize=12)
    ax.set_title(
        f"Dominant-B proportion by (IGH, IGK) {kind}  -  bin20\n"
        f"(overall dominant-B among matched bins = {overall_b_frac:.1%})",
        fontsize=13
    )
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Fraction of bins with dominant = B cell", fontsize=11)

    for i in range(pivot_b.shape[0]):
        for j in range(pivot_b.shape[1]):
            b_val = pivot_b.values[i, j]
            if not np.isnan(b_val) and TARGET_B_LOW <= b_val <= TARGET_B_HIGH:
                rect = plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                     fill=False, edgecolor="gold",
                                     linewidth=2.5)
                ax.add_patch(rect)

    plt.tight_layout()
    out_png = os.path.join(outdir, f"IGH_IGK_{kind}_heatmap_dominant.png")
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"  -> {out_png}")


# =========================================================================
# Main pipeline
# =========================================================================

def run():
    print(f"Output root directory: {BASE_OUT}")
    for k, v in DIRS.items():
        print(f"  {k}: {v}")
    print()

    # 1. bin20 adata
    bin_df, bin_size, convention = load_bin20(H5AD_PATH)
    print(f"  Valid bin20 bin count: {len(bin_df)}\n")

    # 2. meta
    print(f"Loading meta.gz: {META_PATH}")
    df = pd.read_csv(META_PATH, sep=",", compression="gzip")
    print(f"  records: {len(df)}\n")

    # 3. RCTD
    rctd = prepare_rctd(RCTD_PATH, B_RELATED, bin_size, convention)
    # Keep only RCTD records that fall within adata bins
    rctd = pd.merge(rctd, bin_df[["bx", "by"]].drop_duplicates(),
                    on=["bx", "by"], how="inner")
    print(f"  RCTD rows after alignment with bin20: {len(rctd)}\n")

    # 4. Per-chain analysis for the three chains
    stats_by_chain = {}
    for chain in CHAINS:
        print(f"=== {chain} ===")
        stats = build_chain_bin_stats(df, chain, bin_df, bin_size, convention)
        stats_by_chain[chain] = stats
        run_chain_analysis(chain, stats, rctd, DIRS[chain])

    # 5. IGH x IGK heatmaps
    stats_igh = stats_by_chain["IGH"]
    stats_igk = stats_by_chain["IGK"]

    # Overall dominant-B (across all matched bins)
    all_bins = bin_df[["bx", "by"]].drop_duplicates()
    merged_all = pd.merge(all_bins, rctd, on=["bx", "by"], how="left")
    matched_all = merged_all[merged_all["first_type"].notna()]
    overall_b = (matched_all["is_B_dominant"].sum() / len(matched_all)
                 if len(matched_all) else np.nan)
    print(f"\nOverall dominant-B (among bin20 bins matched to RCTD): "
          f"{overall_b:.1%}  ({int(matched_all['is_B_dominant'].sum())}"
          f"/{len(matched_all)})")

    print(f"\n=== IGH x IGK heatmap: unique targetSeq ===")
    sweep_seq = compute_igh_igk_heatmap(
        stats_igh, stats_igk, rctd,
        "n_unique_IGH_targetSeq", "n_unique_IGK_targetSeq",
        SEQ_AXIS_MAX_IGH, SEQ_AXIS_MAX_IGK
    )
    sweep_seq.to_csv(os.path.join(DIRS["seq_hm"],
                                  "IGH_IGK_targetSeq_sweep.csv"),
                     index=False)
    plot_igh_igk_heatmap(sweep_seq, "targetSeq", DIRS["seq_hm"], overall_b)

    print(f"\n=== IGH x IGK heatmap: unique coords ===")
    sweep_coord = compute_igh_igk_heatmap(
        stats_igh, stats_igk, rctd,
        "n_unique_IGH_coords", "n_unique_IGK_coords",
        COORD_AXIS_MAX_IGH, COORD_AXIS_MAX_IGK
    )
    sweep_coord.to_csv(os.path.join(DIRS["coord_hm"],
                                    "IGH_IGK_coords_sweep.csv"),
                       index=False)
    plot_igh_igk_heatmap(sweep_coord, "coords", DIRS["coord_hm"], overall_b)

    print(f"\nDone! All outputs are in: {BASE_OUT}")


if __name__ == "__main__":
    run()
