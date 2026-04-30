"""
Rep2 Step 1:
  对 Rep2 的 cellbin / bin20 / nucleus 各跑一次"不区分 IGHKL"的分析
  完全复用 Rep1 的逻辑, 只换路径
输出:
  Rep2/cellbin_BCR_QC/
  Rep2/bin20_BCR_QC/
  Rep2/nucleus_BCR_QC/
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
import anndata as ad
import tifffile

# ====== Rep2 路径 ======
H5AD_CELLBIN = "/data/scratch/projects/punim1236/spleen_stomics_2/processed/STOmics_mice_spleen_2_brighter_morecores_omejustliker1/outs/visualization/C04278A3.cellbin_1.0.adjusted.h5ad"
H5AD_BIN20   = "/data/scratch/projects/punim1236/spleen_stomics_2/processed/STOmics_mice_spleen_2_brighter_morecores_omejustliker1/outs/visualization/C04278A3.bin20_1.0.h5ad"
META_PATH    = "/data/scratch/projects/punim1236/AAA_chen_here/AAA克隆_XCR/A3_read1.fq.gz.meta.gz"
RCTD_WIDE    = "/data/scratch/projects/punim1236/AAA_chen_here/AAA克隆_XCR/Rep2/RCTD基本数据/RCTD_celltype_proportion_wide_rep2.csv"
# 如果 Rep2 也有 nucleus mask 就填; 否则 nucleus 分析自动跳过
MASK_PATH    = ""   # ← 如果有, 填进来; 没有就保持空字符串

REP2_DIR     = "/data/scratch/projects/punim1236/AAA_chen_here/AAA克隆_XCR/Rep2"
OUT_CELLBIN  = os.path.join(REP2_DIR, "cellbin_BCR_QC")
OUT_BIN20    = os.path.join(REP2_DIR, "bin20_BCR_QC")
OUT_NUCLEUS  = os.path.join(REP2_DIR, "nucleus_BCR_QC")
for d in (OUT_CELLBIN, OUT_BIN20, OUT_NUCLEUS):
    os.makedirs(d, exist_ok=True)

EPS         = 0.01
BCR_CHAINS  = {"IGH", "IGK", "IGL"}
LINE_X_MAX  = 15
MIN_N_LINE  = 5

B_RELATED = [
    "Activated", "Atypical", "B ASC", "B B1", "B follicular", "B MZ", "B1",
    "DZd", "DZp", "Follicular Naive", "GC", "Immature_Transitional",
    "LZ", "LZ_DZ_intermediate", "Memory", "MZ", "PB", "preMem",
    "preMem intermediate", "B GC", "B naive",
]

# ====== 公共: RCTD wide -> dominant_celltype ======
def prepare_rctd_wide(path, b_related):
    print(f"读取 RCTD wide: {path}")
    rctd = pd.read_csv(path)
    rctd["x"] = np.round(rctd["x"]).astype(int)
    rctd["y"] = np.round(rctd["y"]).astype(int)

    meta = {"spot_id", "x", "y"}
    ct_cols = [c for c in rctd.columns
               if c not in meta and pd.api.types.is_numeric_dtype(rctd[c])]
    b_cols  = [c for c in b_related if c in rctd.columns]
    non_b   = [c for c in ct_cols if c not in b_cols]
    print(f"  B-related columns: {len(b_cols)}, non-B: {len(non_b)}")

    coll = rctd[["spot_id", "x", "y"]].copy()
    coll["B cell"] = rctd[b_cols].sum(axis=1)
    for c in non_b:
        coll[c] = rctd[c]

    ct_list = ["B cell"] + non_b
    coll["dominant_celltype"] = coll[ct_list].idxmax(axis=1)
    coll["is_B_dominant"]     = coll["dominant_celltype"] == "B cell"
    print(f"  spots: {len(coll)}, B-dominant: {int(coll['is_B_dominant'].sum())}")
    return coll[["spot_id", "x", "y", "dominant_celltype", "is_B_dominant"]]


# ====== 公共: 折线图 + 5 barplot ======
def _plot_line(stats, rctd, outdir, n_unit, x_max=LINE_X_MAX, min_n=MIN_N_LINE,
               coord_col="cx", coord_col2="cy"):
    sq = stats.copy()
    sq["x"] = np.round(sq[coord_col]).astype(int)
    sq["y"] = np.round(sq[coord_col2]).astype(int)
    merged = pd.merge(sq[["x","y","n_unique_BCR_coords"]],
                      rctd[["x","y","is_B_dominant"]],
                      on=["x","y"], how="left")
    matched = merged[merged["is_B_dominant"].notna()].copy()

    rows = []
    for k in range(x_max + 1):
        sub = matched[matched["n_unique_BCR_coords"] == k]
        if len(sub) == 0:
            rows.append({"n_coords": k, "n_matched": 0,
                         "n_B": 0, "B_proportion": np.nan})
        else:
            n_b = int(sub["is_B_dominant"].sum())
            rows.append({"n_coords": k, "n_matched": len(sub),
                         "n_B": n_b, "B_proportion": n_b/len(sub)})
    line_df = pd.DataFrame(rows).rename(columns={
        "n_matched": f"n_{n_unit}_matched",
        "n_B": "n_B_dominant",
    })
    line_df.to_csv(os.path.join(outdir,"B_proportion_by_coords_line.csv"),
                   index=False)
    print(f"\n折线图数据:\n{line_df}\n")

    overall = (matched["is_B_dominant"].sum()/len(matched)
               if len(matched) else np.nan)

    fig, ax = plt.subplots(figsize=(13, 8))
    enough = line_df[f"n_{n_unit}_matched"] >= min_n
    ax.plot(line_df["n_coords"], line_df["B_proportion"],
            color="#999", lw=1.2, alpha=0.5, ls="--", zorder=2)
    sub_ok = line_df[enough]
    ax.plot(sub_ok["n_coords"], sub_ok["B_proportion"],
            color="#E64B35", lw=2.8, zorder=3, label="B-cell proportion")
    sizes = np.clip(line_df[f"n_{n_unit}_matched"].fillna(0).to_numpy()**0.5*8, 20, 350)
    ax.scatter(line_df["n_coords"], line_df["B_proportion"],
               s=sizes,
               c=["#E64B35" if e else "#999" for e in enough],
               edgecolors="black", linewidths=0.8, zorder=4)
    for _, r in line_df.iterrows():
        if pd.notna(r["B_proportion"]):
            ax.text(r["n_coords"], r["B_proportion"]+0.025,
                    f"n={int(r[f'n_{n_unit}_matched'])}",
                    ha="center", va="bottom", fontsize=9, color="#333")
    ax.axhline(overall, color="#4DBBD5", ls=":", lw=1.8,
               label=f"overall = {overall:.1%}", zorder=1)
    ax.set_xlim(-0.5, x_max+0.5); ax.set_ylim(0, 1.05)
    ax.set_xticks(range(0, x_max+1))
    ax.set_xlabel(f"Number of unique BCR coordinates per {n_unit}", fontsize=13)
    ax.set_ylabel("Fraction with dominant = B cell", fontsize=13)
    ax.set_title(f"Rep2: B proportion vs BCR coords ({n_unit})", fontsize=14)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.grid(True, axis="y", ls="--", alpha=0.3)
    ax.legend(frameon=False, fontsize=11, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir,"B_proportion_by_coords_line.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_barplots(stats, rctd, outdir,
                   coord_col="cx", coord_col2="cy"):
    sq = stats.copy()
    sq["x"] = np.round(sq[coord_col]).astype(int)
    sq["y"] = np.round(sq[coord_col2]).astype(int)
    merged = pd.merge(sq[["x","y","n_unique_BCR_coords"]],
                      rctd[["x","y","dominant_celltype"]],
                      on=["x","y"], how="left")
    m = merged[merged["dominant_celltype"].notna()].copy()
    m["bin_cat"] = m["n_unique_BCR_coords"].clip(upper=4)
    label_map = {0:"0",1:"1",2:"2",3:"3",4:">=4"}
    m["lbl"] = m["bin_cat"].map(label_map)
    def merge_b(c):
        s=str(c).strip()
        if s.startswith("B ") or s in ("B","B1","B cell"): return "B cells"
        return s
    m["ct"] = m["dominant_celltype"].apply(merge_b)
    all_comp=[]
    for b in ["0","1","2","3",">=4"]:
        sub = m[m["lbl"]==b]
        if len(sub)==0: continue
        cd = (sub["ct"].value_counts()
              .rename_axis("cell_type").reset_index(name="count"))
        cd["fraction"] = cd["count"]/cd["count"].sum()
        cd["bin_label"] = b
        all_comp.append(cd)
        main = cd[cd["fraction"]>0.01].copy()
        of = cd.loc[cd["fraction"]<=0.01,"fraction"].sum()
        oc = cd.loc[cd["fraction"]<=0.01,"count"].sum()
        if of>0:
            main = pd.concat([main, pd.DataFrame({
                "cell_type":["Others"],"count":[oc],
                "fraction":[of],"bin_label":[b]})], ignore_index=True)
        main = main.sort_values("fraction")
        plt.figure(figsize=(12,10))
        cols = ["#E64B35" if c=="B cells" else "#4C72B0"
                for c in main["cell_type"]]
        plt.barh(main["cell_type"], main["fraction"], color=cols, height=0.8)
        plt.title(f"Rep2 dominant celltype: BCR coords = {b} (n={len(sub)})",
                  fontsize=14)
        plt.xlabel("Fraction"); plt.ylabel("Dominant cell type")
        for y_,x_,c_ in zip(main["cell_type"], main["fraction"], main["count"]):
            plt.text(x_, y_, f"  {c_} ({x_*100:.2f}%)",
                     ha="left", va="center", fontsize=11)
        ax=plt.gca()
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        plt.tight_layout()
        sb = b.replace(">=","ge")
        plt.savefig(os.path.join(outdir,
            f"dominant_composition_bin{sb}_BCRcoords.png"),
            dpi=200, bbox_inches="tight")
        plt.close()
    if all_comp:
        pd.concat(all_comp, ignore_index=True).to_csv(
            os.path.join(outdir,"dominant_composition_all_bins_BCRcoords.csv"),
            index=False)


# =====================================================================
# A. CELLBIN
# =====================================================================
def run_cellbin():
    print("\n" + "="*70 + "\nRep2 CELLBIN\n" + "="*70)
    adata = ad.read_h5ad(H5AD_CELLBIN)
    if {"x","y"}.issubset(adata.obs.columns):
        C = adata.obs[["x","y"]].to_numpy()
    else:
        C = adata.obsm["spatial"]
    B = adata.obsm["cell_border"]
    print(f"  cell bins: {len(C)}")

    df = pd.read_csv(META_PATH, sep=",", compression="gzip")
    print(f"  meta records: {len(df)}")

    sub = df[df["topChains"].astype(str).str.strip().isin(BCR_CHAINS)].copy()
    sub = sub.dropna(subset=["x","y"])
    print(f"  BCR records: {len(sub)}")
    px = sub["x"].to_numpy(np.float32); py = sub["y"].to_numpy(np.float32)
    points = np.column_stack([px,py])

    records=[]
    for i in range(B.shape[0]):
        pts = B[i]; valid = ~(pts==32767).any(axis=1)
        pts = pts[valid].astype(np.float32)
        if pts.shape[0]<3:
            records.append({"cell_idx":i,"cx":C[i,0],"cy":C[i,1],
                            "n_unique_BCR_coords":0}); continue
        poly = C[i] + pts[:,[0,1]]
        bbox = ((px>=poly[:,0].min())&(px<=poly[:,0].max())&
                (py>=poly[:,1].min())&(py<=poly[:,1].max()))
        if not bbox.any():
            n=0
        else:
            inside = MplPath(poly).contains_points(points[bbox], radius=-EPS)
            if not inside.any():
                n=0
            else:
                sel = np.where(bbox)[0][inside]
                n = int(np.unique(points[sel],axis=0).shape[0])
        records.append({"cell_idx":i,"cx":C[i,0],"cy":C[i,1],
                        "n_unique_BCR_coords":n})
    stats = pd.DataFrame(records)
    stats.to_csv(os.path.join(OUT_CELLBIN,"BCR_combined_cellbin_stats.csv"),
                 index=False)

    rctd = prepare_rctd_wide(RCTD_WIDE, B_RELATED)
    _plot_line(stats, rctd, OUT_CELLBIN, n_unit="cells")
    _plot_barplots(stats, rctd, OUT_CELLBIN)


# =====================================================================
# B. BIN20
# =====================================================================
def run_bin20():
    print("\n" + "="*70 + "\nRep2 BIN20\n" + "="*70)
    adata = ad.read_h5ad(H5AD_BIN20)
    obs=adata.obs
    if {"x","y"}.issubset(obs.columns):
        bx=obs["x"].to_numpy(float); by=obs["y"].to_numpy(float)
    else:
        bx=adata.obsm["spatial"][:,0].astype(float)
        by=adata.obsm["spatial"][:,1].astype(float)
    ux=np.sort(np.unique(bx)); uy=np.sort(np.unique(by))
    bs = float(round((np.min(np.diff(ux))+np.min(np.diff(uy)))/2))
    fx=(bx/bs)-np.floor(bx/bs)
    conv = "corner" if np.mean((fx<0.05)|(fx>0.95)) >= np.mean(np.abs(fx-0.5)<0.05) else "center"
    print(f"  bin_size={bs}, convention={conv}")
    bin_df = pd.DataFrame({"bx":bx,"by":by})

    def p2b(px,py):
        if conv=="corner":
            return np.floor(px/bs)*bs, np.floor(py/bs)*bs
        return np.floor(px/bs)*bs+bs/2, np.floor(py/bs)*bs+bs/2

    df = pd.read_csv(META_PATH, sep=",", compression="gzip")
    sub = df[df["topChains"].astype(str).str.strip().isin(BCR_CHAINS)].copy()
    sub = sub.dropna(subset=["x","y"])
    rx,ry = p2b(sub["x"].to_numpy(float), sub["y"].to_numpy(float))
    sub["bx"]=rx; sub["by"]=ry
    sub = pd.merge(sub, bin_df.drop_duplicates(), on=["bx","by"], how="inner")
    coord_u = sub[["bx","by","x","y"]].drop_duplicates()
    n_coord = (coord_u.groupby(["bx","by"]).size()
               .reset_index(name="n_unique_BCR_coords"))
    stats = bin_df.drop_duplicates().copy()
    stats = pd.merge(stats, n_coord, on=["bx","by"], how="left")
    stats["n_unique_BCR_coords"] = stats["n_unique_BCR_coords"].fillna(0).astype(int)
    stats.to_csv(os.path.join(OUT_BIN20,"BCR_combined_bin20_stats.csv"), index=False)

    # 把 RCTD spot 坐标也映射到 bin 参考点 (这样和 stats 能 merge 上)
    rctd = prepare_rctd_wide(RCTD_WIDE, B_RELATED)
    rrx, rry = p2b(rctd["x"].to_numpy(float), rctd["y"].to_numpy(float))
    rctd["x"] = rrx.astype(int); rctd["y"] = rry.astype(int)

    # bin20 stats 用 bx/by 当 (x,y); 改成统一字段名给 _plot_line/_plot_barplots
    stats_renamed = stats.rename(columns={"bx":"cx","by":"cy"})
    _plot_line(stats_renamed, rctd, OUT_BIN20, n_unit="bins")
    _plot_barplots(stats_renamed, rctd, OUT_BIN20)


# =====================================================================
# C. NUCLEUS-strict
# =====================================================================
def run_nucleus():
    if not MASK_PATH or not os.path.exists(MASK_PATH):
        print("\n>>> Skipping nucleus-strict: MASK_PATH 未提供或文件不存在")
        return
    print("\n" + "="*70 + "\nRep2 NUCLEUS-STRICT\n" + "="*70)

    adata = ad.read_h5ad(H5AD_CELLBIN)
    if {"x","y"}.issubset(adata.obs.columns):
        C = adata.obs[["x","y"]].to_numpy()
    else:
        C = adata.obsm["spatial"]
    B = adata.obsm["cell_border"]

    mask = tifffile.imread(MASK_PATH)
    if mask.ndim==3:
        mask = mask[...,0] if mask.shape[-1]<=4 else mask[0]
    nuc = mask>0
    H, W = nuc.shape
    print(f"  nucleus mask: {nuc.shape}, foreground: {int(nuc.sum())}")

    df = pd.read_csv(META_PATH, sep=",", compression="gzip")
    sub = df[df["topChains"].astype(str).str.strip().isin(BCR_CHAINS)].copy()
    sub = sub.dropna(subset=["x","y"])
    px = sub["x"].to_numpy(np.float32); py = sub["y"].to_numpy(np.float32)
    points = np.column_stack([px,py])
    pxi = np.round(px).astype(int); pyi = np.round(py).astype(int)
    inb = (pxi>=0)&(pxi<W)&(pyi>=0)&(pyi<H)
    in_nuc = np.zeros(len(sub), bool)
    in_nuc[inb] = nuc[pyi[inb], pxi[inb]]

    records=[]
    for i in range(B.shape[0]):
        pts = B[i]; valid = ~(pts==32767).any(axis=1)
        pts = pts[valid].astype(np.float32)
        if pts.shape[0]<3:
            records.append({"cell_idx":i,"cx":C[i,0],"cy":C[i,1],
                            "n_unique_BCR_coords":0}); continue
        poly = C[i] + pts[:,[0,1]]
        bbox = ((px>=poly[:,0].min())&(px<=poly[:,0].max())&
                (py>=poly[:,1].min())&(py<=poly[:,1].max())& in_nuc)
        if not bbox.any():
            n=0
        else:
            inside = MplPath(poly).contains_points(points[bbox], radius=-EPS)
            if not inside.any():
                n=0
            else:
                sel = np.where(bbox)[0][inside]
                n = int(np.unique(points[sel],axis=0).shape[0])
        records.append({"cell_idx":i,"cx":C[i,0],"cy":C[i,1],
                        "n_unique_BCR_coords":n})
    stats = pd.DataFrame(records)
    stats.to_csv(os.path.join(OUT_NUCLEUS,"BCR_combined_nucleus_stats.csv"),
                 index=False)

    rctd = prepare_rctd_wide(RCTD_WIDE, B_RELATED)
    _plot_line(stats, rctd, OUT_NUCLEUS, n_unit="cells")
    _plot_barplots(stats, rctd, OUT_NUCLEUS)


if __name__ == "__main__":
    run_cellbin()
    run_bin20()
    run_nucleus()
    print("\nStep 1 done.")