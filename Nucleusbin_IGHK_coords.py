"""
===========================================================================
NUCLEUS-STRICT version BCR QC analysis (IGH / IGK / IGL)
===========================================================================
Analysis on cell bins, but when counting unique IGH/IGK/IGL targetSeq /
unique coords, coordinates must satisfy BOTH:
  1) Fall within the cell bin polygon
  2) Fall on the nucleus mask foreground (mask > 0)

For each chain:
  1. Per-chain unique targetSeq frequency + spatial binning + dominant composition
  2. Per-chain unique coords    frequency + spatial binning + dominant composition
Additionally:
  3. Dominant-B heatmap for IGH x IGK unique targetSeq
  4. Dominant-B heatmap for IGH x IGK unique coords

RCTD convention (same as the cell-bin version):
  - Use RCTD_celltype_proportion_wide.csv
  - Sum the B_RELATED columns -> "B cell"
  - Keep the rest -> dominant_celltype = argmax
  - dominant_celltype == "B cell" is treated as B-dominant

Output directory structure:
  BASE_OUT/
    nucleus_IGH/
    nucleus_IGK/
    nucleus_IGL/
    nucleus_sequence_heatmap/
    nucleus_coordinates_heatmap/
===========================================================================
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.path import Path
from itertools import product
import anndata as ad
import tifffile


# =========================================================================
# User configuration
# =========================================================================
H5AD_PATH  = "/data/scratch/projects/punim1236/spleen_stomics/processed/spleen_out/outs/visualization/C04278C6.cellbin_1.0.adjusted.h5ad"
META_PATH  = "/data/scratch/projects/punim1236/AAA_chen_here/AAA克隆_XCR/C6_read1.fq.gz.meta.gz"
RCTD_PATH  = "/data/scratch/projects/punim1236/AAA_chen_here/RCTD_celltype_proportion_wide.csv"
MASK_PATH  = "/data/scratch/projects/punim1236/spleen_stomics/processed/spleen_out/outs/image/C04278C6_DAPI_mask.tif"

BASE_OUT = "/data/scratch/projects/punim1236/AAA_chen_here/AAA克隆_XCR/扫描高占比B/Nucleus_bin扫描高占比"

EPS = 0.01

CHAINS = ["IGH", "IGK", "IGL"]
SEQ_COL = "targetSequences"

B_RELATED = [
    "Activated", "Atypical", "B ASC", "B B1", "B follicular", "B MZ", "B1",
    "DZd", "DZp", "Follicular Naive", "GC", "Immature_Transitional",
    "LZ", "LZ_DZ_intermediate", "Memory", "MZ", "PB", "preMem", "preMem intermediate"
]

# Heatmap axis range (consistent with the bin20 version; adjust later if needed)
SEQ_AXIS_MAX_IGH   = 8
SEQ_AXIS_MAX_IGK   = 10
COORD_AXIS_MAX_IGH = 6
COORD_AXIS_MAX_IGK = 12

# Gold-border highlight range
TARGET_B_LOW  = 0.80
TARGET_B_HIGH = 0.90

# Subfolders
DIRS = {
    "IGH":      os.path.join(BASE_OUT, "nucleus_IGH"),
    "IGK":      os.path.join(BASE_OUT, "nucleus_IGK"),
    "IGL":      os.path.join(BASE_OUT, "nucleus_IGL"),
    "seq_hm":   os.path.join(BASE_OUT, "nucleus_sequence_heatmap"),
    "coord_hm": os.path.join(BASE_OUT, "nucleus_coordinates_heatmap"),
}
for d in DIRS.values():
    os.makedirs(d, exist_ok=True)


# =========================================================================
# Step 1: Load cell bin + nucleus mask
# =========================================================================

def get_cell_geometry(adata_cellbin):
    if {"x", "y"}.issubset(adata_cellbin.obs.columns):
        C = adata_cellbin.obs[["x", "y"]].to_numpy()
    else:
        C = adata_cellbin.obsm["spatial"]
    B = adata_cellbin.obsm["cell_border"]
    return C, B


def load_nucleus_mask_bool(mask_path):
    mask = tifffile.imread(mask_path)
    if mask.ndim == 3:
        mask = mask[..., 0] if mask.shape[-1] <= 4 else mask[0]
    nuc_bool = mask > 0
    print(f"  nucleus mask shape: {nuc_bool.shape}, "
          f"foreground px: {int(nuc_bool.sum())}")
    return nuc_bool


# =========================================================================
# Step 2: Per-cell statistics function for nucleus-strict
#         -> Returns both unique targetSeq count and unique coords count
# =========================================================================

def count_chain_nucleus_strict(df, chain, C, B, nuc_bool, eps=EPS):
    """
    Returns a DataFrame with one row per cell bin:
      cell_idx, cx, cy,
      n_unique_{chain}_targetSeq,
      n_unique_{chain}_coords
    Only counts reads falling inside (cell polygon AND nucleus mask foreground).
    """
    sub = df[df["topChains"].astype(str).str.strip() == chain].copy()
    sub = sub.dropna(subset=["x", "y"]).copy()

    # Clean the targetSeq column
    if SEQ_COL in sub.columns:
        sub[SEQ_COL] = sub[SEQ_COL].astype(str).str.strip()
        bad = {"", "nan", "None", "none", "NA"}
        sub["_has_seq"] = ~sub[SEQ_COL].isin(bad)
    else:
        sub["_has_seq"] = False

    px = sub["x"].to_numpy(dtype=np.float32)
    py = sub["y"].to_numpy(dtype=np.float32)
    seqs = sub[SEQ_COL].to_numpy() if SEQ_COL in sub.columns else np.array([""] * len(sub))
    has_seq = sub["_has_seq"].to_numpy()
    points = np.column_stack([px, py])

    # Global: which BCR points fall on the nucleus mask foreground
    H, W = nuc_bool.shape
    px_int = np.round(px).astype(np.int32)
    py_int = np.round(py).astype(np.int32)
    in_bounds = (
        (px_int >= 0) & (px_int < W) &
        (py_int >= 0) & (py_int < H)
    )
    point_in_nuc = np.zeros(points.shape[0], dtype=bool)
    point_in_nuc[in_bounds] = nuc_bool[py_int[in_bounds], px_int[in_bounds]]

    records = []
    for i in range(B.shape[0]):
        pts = B[i]
        valid = ~(pts == 32767).any(axis=1)
        pts = pts[valid].astype(np.float32)
        if pts.shape[0] < 3:
            records.append({
                "cell_idx": i,
                "cx": C[i, 0],
                "cy": C[i, 1],
                f"n_unique_{chain}_targetSeq": 0,
                f"n_unique_{chain}_coords": 0,
            })
            continue

        poly_xy = C[i] + pts[:, [0, 1]]
        xmin, xmax = poly_xy[:, 0].min(), poly_xy[:, 0].max()
        ymin, ymax = poly_xy[:, 1].min(), poly_xy[:, 1].max()

        mask = (
            (px >= xmin) & (px <= xmax) &
            (py >= ymin) & (py <= ymax) &
            point_in_nuc
        )
        if not np.any(mask):
            n_seq = 0
            n_coord = 0
        else:
            inside = Path(poly_xy).contains_points(points[mask], radius=-eps)
            if not np.any(inside):
                n_seq = 0
                n_coord = 0
            else:
                sel_idx = np.where(mask)[0][inside]
                # unique coords
                n_coord = int(np.unique(points[sel_idx], axis=0).shape[0])
                # unique targetSeq (only for entries with has_seq)
                sel_has_seq = has_seq[sel_idx]
                if sel_has_seq.any():
                    n_seq = int(np.unique(seqs[sel_idx][sel_has_seq]).shape[0])
                else:
                    n_seq = 0

        records.append({
            "cell_idx": i,
            "cx": C[i, 0],
            "cy": C[i, 1],
            f"n_unique_{chain}_targetSeq": n_seq,
            f"n_unique_{chain}_coords": n_coord,
        })
    return pd.DataFrame(records)


# =========================================================================
# Step 3: RCTD + dominant_celltype + is_B_dominant
# =========================================================================

def prepare_rctd_with_dominant(rctd_path, b_related):
    print(f"Loading RCTD: {rctd_path}")
    rctd = pd.read_csv(rctd_path)
    rctd["x"] = np.round(rctd["x"]).astype(int)
    rctd["y"] = np.round(rctd["y"]).astype(int)

    meta = {"spot_id", "x", "y"}
    ct_cols = [c for c in rctd.columns
               if c not in meta and pd.api.types.is_numeric_dtype(rctd[c])]
    b_cols = [c for c in b_related if c in rctd.columns]
    non_b = [c for c in ct_cols if c not in b_cols]

    collapsed = rctd[["spot_id", "x", "y"]].copy()
    collapsed["B cell"] = rctd[b_cols].sum(axis=1)
    for c in non_b:
        collapsed[c] = rctd[c]

    collapsed_ct = ["B cell"] + non_b
    collapsed["dominant_celltype"] = collapsed[collapsed_ct].idxmax(axis=1)
    collapsed["is_B_dominant"] = collapsed["dominant_celltype"] == "B cell"
    print(f"  RCTD spots: {len(collapsed)}, "
          f"B-dominant: {int(collapsed['is_B_dominant'].sum())}")
    return collapsed[["spot_id", "x", "y", "dominant_celltype", "is_B_dominant"]]


# =========================================================================
# Step 4: Three per-chain plots (frequency / spatial / dominant composition)
# =========================================================================

def plot_chain_frequency(stats, chain, value_col, xlabel_kind, outdir):
    freq = stats[value_col].value_counts().sort_index().reset_index()
    freq.columns = [value_col, "frequency"]
    freq.to_csv(os.path.join(outdir, f"{chain}_{xlabel_kind}_frequency.csv"),
                index=False)

    plot_freq = freq[freq[value_col] <= 30].copy()
    total = plot_freq["frequency"].sum()

    plt.figure(figsize=(12, 10))
    plt.barh(plot_freq[value_col], plot_freq["frequency"], height=0.8)
    plt.ylabel(f"Number of unique {chain} {xlabel_kind} per cell bin (nucleus-only)",
               fontsize=13)
    plt.xlabel("Number of cell bins", fontsize=13)
    plt.title(
        f"Frequency distribution of unique {chain} {xlabel_kind} per cell bin\n"
        f"(nucleus-strict)",
        fontsize=15
    )
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
    df_plot = stats[["cx", "cy", value_col]].copy()
    df_plot["bin_cat"] = df_plot[value_col].clip(upper=4)

    color_map = {0: "lightgrey", 1: "#4DBBD5", 2: "#00A087",
                 3: "#E64B35",   4: "#3C5488"}
    label_map = {0: "0", 1: "1", 2: "2", 3: "3", 4: ">=4"}

    plt.figure(figsize=(10, 10))
    for k in [0, 1, 2, 3, 4]:
        sub = df_plot[df_plot["bin_cat"] == k]
        plt.scatter(sub["cx"], sub["cy"], s=3, c=color_map[k],
                    label=f"{label_map[k]} (n={len(sub)})",
                    linewidths=0, alpha=0.9)
    plt.gca().set_aspect("equal", "box")
    plt.gca().invert_yaxis()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(
        f"Spatial distribution of {chain} {xlabel_kind} count per cell bin "
        f"(nucleus-strict)"
    )
    plt.legend(frameon=False, markerscale=3)
    plt.tight_layout()
    png = os.path.join(outdir, f"{chain}_{xlabel_kind}_spatial_bins.png")
    plt.savefig(png, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"    -> {png}")


def plot_bin_dominant_composition(stats, rctd, chain, value_col,
                                  xlabel_kind, outdir):
    """
    For each cell bin grouped by bin_cat (0/1/2/3/>=4), after matching with RCTD:
      - Compute the composition of their dominant_celltype
      - Merge all "B-prefixed" types into "B cells"
        (consistent with the bin20 version)
    """
    # Round the cell bin centroids from `stats` and merge with RCTD (x, y)
    stats_q = stats.copy()
    stats_q["x"] = np.round(stats_q["cx"]).astype(int)
    stats_q["y"] = np.round(stats_q["cy"]).astype(int)

    merged = pd.merge(
        stats_q[["x", "y", value_col]],
        rctd[["x", "y", "dominant_celltype"]],
        on=["x", "y"], how="left"
    )
    matched = merged[merged["dominant_celltype"].notna()].copy()
    matched["bin_cat"] = matched[value_col].clip(upper=4)
    bin_label_map = {0: "0", 1: "1", 2: "2", 3: "3", 4: ">=4"}
    matched["bin_cat_label"] = matched["bin_cat"].map(bin_label_map)

    # Merge B-prefixed types.
    # Note: in prepare_rctd_with_dominant we've already collapsed B subtypes
    # into "B cell", so in practice we only see the single "B cell" category
    # here. Keeping this function improves label readability.
    def _merge_b_types(ct):
        s = str(ct).strip()
        if s.startswith("B ") or s == "B" or s == "B1" or s == "B cell":
            return "B cells"
        return s
    matched["cell_type_merged"] = matched["dominant_celltype"].apply(_merge_b_types)

    all_comp = []
    for b in ["0", "1", "2", "3", ">=4"]:
        sub = matched[matched["bin_cat_label"] == b]
        if len(sub) == 0:
            continue
        count_df = (sub["cell_type_merged"].value_counts()
                    .rename_axis("cell_type").reset_index(name="count"))
        count_df["fraction"] = count_df["count"] / count_df["count"].sum()
        count_df["bin_label"] = b
        all_comp.append(count_df)

        # Individual barplot
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
        colors = ["#E64B35" if ct == "B cells" else "#4C72B0"
                  for ct in main_df["cell_type"]]
        plt.barh(main_df["cell_type"], main_df["fraction"],
                 height=0.8, color=colors)
        plt.ylabel("Dominant cell type (B-types merged)", fontsize=13)
        plt.xlabel("Fraction", fontsize=13)
        plt.title(
            f"Dominant cell-type distribution in cell bins with {b} "
            f"unique {chain} {xlabel_kind}\n(nucleus-strict)",
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
        png = os.path.join(
            outdir,
            f"{chain}_{xlabel_kind}_bin{safe_b}_dominant_composition.png"
        )
        plt.savefig(png, dpi=200, bbox_inches="tight")
        plt.close()

    if all_comp:
        out_csv = os.path.join(
            outdir, f"{chain}_{xlabel_kind}_dominant_composition_all.csv"
        )
        pd.concat(all_comp, ignore_index=True).to_csv(out_csv, index=False)
        print(f"    -> {out_csv}")


def run_chain_analysis(chain, stats, rctd, outdir):
    print(f"\n--- {chain} single-chain analysis (nucleus-strict) ---")
    # targetSeq
    plot_chain_frequency(stats, chain,
                         f"n_unique_{chain}_targetSeq", "targetSeq", outdir)
    plot_chain_spatial_bins(stats, chain,
                            f"n_unique_{chain}_targetSeq", "targetSeq", outdir)
    plot_bin_dominant_composition(stats, rctd, chain,
                                  f"n_unique_{chain}_targetSeq",
                                  "targetSeq", outdir)
    # coords
    plot_chain_frequency(stats, chain,
                         f"n_unique_{chain}_coords", "coords", outdir)
    plot_chain_spatial_bins(stats, chain,
                            f"n_unique_{chain}_coords", "coords", outdir)
    plot_bin_dominant_composition(stats, rctd, chain,
                                  f"n_unique_{chain}_coords",
                                  "coords", outdir)

    out_csv = os.path.join(outdir, f"{chain}_nucleus_stats.csv")
    stats.to_csv(out_csv, index=False)
    print(f"    -> {out_csv}")


# =========================================================================
# Step 5: IGH x IGK heatmap
# =========================================================================

def compute_igh_igk_heatmap(stats_igh, stats_igk, rctd,
                            col_igh, col_igk,
                            igh_max, igk_max):
    joint = pd.merge(
        stats_igh[["cell_idx", "cx", "cy", col_igh]],
        stats_igk[["cell_idx", col_igk]],
        on="cell_idx", how="outer"
    )
    joint[col_igh] = joint[col_igh].fillna(0).astype(int)
    joint[col_igk] = joint[col_igk].fillna(0).astype(int)

    joint["x"] = np.round(joint["cx"]).astype("Int64")
    joint["y"] = np.round(joint["cy"]).astype("Int64")

    joint = pd.merge(
        joint, rctd[["x", "y", "spot_id", "dominant_celltype", "is_B_dominant"]],
        on=["x", "y"], how="left"
    )

    results = []
    for igh_n, igk_n in product(range(igh_max + 1), range(igk_max + 1)):
        sub = joint[(joint[col_igh] == igh_n) & (joint[col_igk] == igk_n)]
        n_total = len(sub)
        matched = sub[sub["spot_id"].notna()]
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
                ax.text(j, i, "–", ha="center", va="center",
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
    ax.set_xlabel(f"Unique IGK {kind} per cell bin (nucleus-only)", fontsize=12)
    ax.set_ylabel(f"Unique IGH {kind} per cell bin (nucleus-only)", fontsize=12)
    ax.set_title(
        f"Dominant-B proportion by (IGH, IGK) {kind}  -  nucleus-strict\n"
        f"(overall dominant-B among matched bins = {overall_b_frac:.1%})",
        fontsize=13
    )
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Fraction of cell bins with dominant = B cell", fontsize=11)

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

    # 1. h5ad + mask
    print("Loading h5ad ...")
    adata_cellbin = ad.read_h5ad(H5AD_PATH)
    print(f"  cell bins: {adata_cellbin.n_obs}")
    C, B = get_cell_geometry(adata_cellbin)

    print("Loading nucleus mask ...")
    nuc_bool = load_nucleus_mask_bool(MASK_PATH)
    print()

    # 2. meta
    print(f"Loading meta.gz: {META_PATH}")
    df = pd.read_csv(META_PATH, sep=",", compression="gzip")
    print(f"  records: {len(df)}\n")

    # 3. RCTD
    rctd = prepare_rctd_with_dominant(RCTD_PATH, B_RELATED)
    print()

    # 4. Per-chain analysis (nucleus-strict)
    stats_by_chain = {}
    for chain in CHAINS:
        print(f"=== {chain} (nucleus-strict statistics) ===")
        stats = count_chain_nucleus_strict(df, chain, C, B, nuc_bool, EPS)
        stats_by_chain[chain] = stats
        run_chain_analysis(chain, stats, rctd, DIRS[chain])

    # 5. Overall dominant-B proportion (across all cell bins matched to RCTD)
    any_stats = stats_by_chain[CHAINS[0]][["cell_idx", "cx", "cy"]].copy()
    any_stats["x"] = np.round(any_stats["cx"]).astype(int)
    any_stats["y"] = np.round(any_stats["cy"]).astype(int)
    merged_all = pd.merge(
        any_stats, rctd[["x", "y", "dominant_celltype", "is_B_dominant"]],
        on=["x", "y"], how="left"
    )
    matched_all = merged_all[merged_all["dominant_celltype"].notna()]
    overall_b = (matched_all["is_B_dominant"].sum() / len(matched_all)
                 if len(matched_all) else np.nan)
    print(f"\nOverall dominant-B (among cell bins matched to RCTD): "
          f"{overall_b:.1%}  ({int(matched_all['is_B_dominant'].sum())}"
          f"/{len(matched_all)})")

    # 6. IGH x IGK heatmaps
    stats_igh = stats_by_chain["IGH"]
    stats_igk = stats_by_chain["IGK"]

    print(f"\n=== IGH x IGK heatmap: unique targetSeq (nucleus-strict) ===")
    sweep_seq = compute_igh_igk_heatmap(
        stats_igh, stats_igk, rctd,
        "n_unique_IGH_targetSeq", "n_unique_IGK_targetSeq",
        SEQ_AXIS_MAX_IGH, SEQ_AXIS_MAX_IGK
    )
    sweep_seq.to_csv(os.path.join(DIRS["seq_hm"],
                                  "IGH_IGK_targetSeq_sweep.csv"),
                     index=False)
    plot_igh_igk_heatmap(sweep_seq, "targetSeq", DIRS["seq_hm"], overall_b)

    print(f"\n=== IGH x IGK heatmap: unique coords (nucleus-strict) ===")
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
