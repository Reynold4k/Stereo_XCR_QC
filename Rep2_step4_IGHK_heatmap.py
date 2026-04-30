"""
Rep2 Step 4: cellbin 上 IGH×IGK 双链热图
  - targetSeq 版本: y=IGH unique targetSeq, x=IGK unique targetSeq
  - coords 版本: y=IGH unique coords, x=IGK unique coords
  cell-color = dominant=B 的比例
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.path import Path as MplPath
from itertools import product

H5AD_CELLBIN = "/data/scratch/projects/punim1236/spleen_stomics_2/processed/STOmics_mice_spleen_2_brighter_morecores_omejustliker1/outs/visualization/C04278A3.cellbin_1.0.adjusted.h5ad"
META_PATH    = "/data/scratch/projects/punim1236/AAA_chen_here/AAA克隆_XCR/A3_read1.fq.gz.meta.gz"
RCTD_WIDE    = "/data/scratch/projects/punim1236/AAA_chen_here/AAA克隆_XCR/Rep2/RCTD基本数据/RCTD_celltype_proportion_wide_rep2.csv"
OUT_DIR      = "/data/scratch/projects/punim1236/AAA_chen_here/AAA克隆_XCR/Rep2/IGH_IGK_heatmap_cellbin"
os.makedirs(OUT_DIR, exist_ok=True)

EPS = 0.01
SEQ_AXIS_MAX_IGH, SEQ_AXIS_MAX_IGK = 8, 10
COORD_AXIS_MAX_IGH, COORD_AXIS_MAX_IGK = 6, 12
TARGET_B_LOW, TARGET_B_HIGH = 0.80, 0.90
B_RELATED = [
    "Activated","Atypical","B ASC","B B1","B follicular","B MZ","B1",
    "DZd","DZp","Follicular Naive","GC","Immature_Transitional",
    "LZ","LZ_DZ_intermediate","Memory","MZ","PB","preMem",
    "preMem intermediate","B GC","B naive",
]
SEQ_COL = "targetSequences"

# ====== 1. cellbin ======
print("Loading cellbin ...")
adata = ad.read_h5ad(H5AD_CELLBIN)
if {"x","y"}.issubset(adata.obs.columns):
    C = adata.obs[["x","y"]].to_numpy()
else:
    C = adata.obsm["spatial"]
B = adata.obsm["cell_border"]

# ====== 2. RCTD wide -> dominant_celltype ======
print("Loading RCTD ...")
rctd = pd.read_csv(RCTD_WIDE)
rctd["x"] = np.round(rctd["x"]).astype(int)
rctd["y"] = np.round(rctd["y"]).astype(int)
meta = {"spot_id","x","y"}
ct_cols = [c for c in rctd.columns if c not in meta and pd.api.types.is_numeric_dtype(rctd[c])]
b_cols = [c for c in B_RELATED if c in rctd.columns]
non_b  = [c for c in ct_cols if c not in b_cols]
coll = rctd[["spot_id","x","y"]].copy()
coll["B cell"] = rctd[b_cols].sum(axis=1)
for c in non_b: coll[c] = rctd[c]
coll["dominant_celltype"] = coll[["B cell"]+non_b].idxmax(axis=1)
coll["is_B_dominant"] = coll["dominant_celltype"]=="B cell"
rctd_d = coll[["spot_id","x","y","dominant_celltype","is_B_dominant"]]

# ====== 3. per-cell IGH/IGK 计数 ======
df = pd.read_csv(META_PATH, sep=",", compression="gzip")
df["topChains"] = df["topChains"].astype(str).str.strip()
if SEQ_COL in df.columns:
    df[SEQ_COL] = df[SEQ_COL].astype(str).str.strip()

def count_chain(df, chain):
    sub = df[df["topChains"]==chain].dropna(subset=["x","y"]).copy()
    if SEQ_COL in sub.columns:
        bad = {"","nan","None","none","NA"}
        sub["_has_seq"] = ~sub[SEQ_COL].isin(bad)
    else:
        sub["_has_seq"] = False
    px = sub["x"].to_numpy(np.float32); py = sub["y"].to_numpy(np.float32)
    seqs = sub[SEQ_COL].to_numpy() if SEQ_COL in sub.columns else np.array([""]*len(sub))
    has_seq = sub["_has_seq"].to_numpy()
    points = np.column_stack([px,py])

    rec=[]
    for i in range(B.shape[0]):
        pts=B[i]; v=~(pts==32767).any(axis=1); pts=pts[v].astype(np.float32)
        if pts.shape[0]<3:
            rec.append({"cell_idx":i,"cx":C[i,0],"cy":C[i,1],
                        f"n_{chain}_seq":0, f"n_{chain}_coord":0}); continue
        poly = C[i]+pts[:,[0,1]]
        bbox = ((px>=poly[:,0].min())&(px<=poly[:,0].max())&
                (py>=poly[:,1].min())&(py<=poly[:,1].max()))
        if not bbox.any():
            n_s=n_c=0
        else:
            inside = MplPath(poly).contains_points(points[bbox], radius=-EPS)
            if not inside.any():
                n_s=n_c=0
            else:
                sel = np.where(bbox)[0][inside]
                n_c = int(np.unique(points[sel],axis=0).shape[0])
                hs = has_seq[sel]
                n_s = int(np.unique(seqs[sel][hs]).shape[0]) if hs.any() else 0
        rec.append({"cell_idx":i,"cx":C[i,0],"cy":C[i,1],
                    f"n_{chain}_seq":n_s, f"n_{chain}_coord":n_c})
    return pd.DataFrame(rec)

print("Counting IGH ..."); s_h = count_chain(df, "IGH")
print("Counting IGK ..."); s_k = count_chain(df, "IGK")

joint = pd.merge(s_h, s_k.drop(columns=["cx","cy"]), on="cell_idx", how="outer")
for c in ["n_IGH_seq","n_IGH_coord","n_IGK_seq","n_IGK_coord"]:
    joint[c] = joint[c].fillna(0).astype(int)
joint["x"] = np.round(joint["cx"]).astype("Int64")
joint["y"] = np.round(joint["cy"]).astype("Int64")
joint = pd.merge(joint, rctd_d[["x","y","spot_id","dominant_celltype","is_B_dominant"]],
                 on=["x","y"], how="left")

# ====== 4. heatmap 计算 ======
def sweep(joint, col_h, col_k, h_max, k_max):
    rows=[]
    for hi, ki in product(range(h_max+1), range(k_max+1)):
        sub = joint[(joint[col_h]==hi) & (joint[col_k]==ki)]
        n_total=len(sub)
        m = sub[sub["spot_id"].notna()]
        n_match = len(m)
        n_b = int(m["is_B_dominant"].sum()) if n_match else 0
        rows.append({"IGH_n":hi,"IGK_n":ki,"n":n_total,
                     "n_matched":n_match,"n_B":n_b,
                     "B_fraction": n_b/n_match if n_match else np.nan})
    return pd.DataFrame(rows)

def plot_heatmap(sweep_df, kind, outdir, overall_b):
    pv_b   = sweep_df.pivot(index="IGH_n", columns="IGK_n", values="B_fraction").sort_index()
    pv_m   = sweep_df.pivot(index="IGH_n", columns="IGK_n", values="n_matched").sort_index()
    pv_dom = sweep_df.pivot(index="IGH_n", columns="IGK_n", values="n_B").sort_index()

    fig, ax = plt.subplots(figsize=(16,9))
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "bcell", ["#D9D9D9","#4DBBD5","#00A087","#E64B35"], N=256)
    norm = mcolors.TwoSlopeNorm(vmin=0.3, vcenter=0.8, vmax=1.0)
    im = ax.imshow(pv_b.values, cmap=cmap, norm=norm, aspect="auto", origin="lower")

    for i in range(pv_b.shape[0]):
        for j in range(pv_b.shape[1]):
            v = pv_b.values[i,j]; m = pv_m.values[i,j]; d = pv_dom.values[i,j]
            if np.isnan(v) or m==0:
                ax.text(j, i, "–", ha="center", va="center", fontsize=8, color="gray")
                continue
            color = "white" if v>0.75 else "black"
            weight= "bold" if TARGET_B_LOW<=v<=TARGET_B_HIGH else "normal"
            ax.text(j, i, f"{v:.1%}\n{int(d)}/{int(m)}",
                    ha="center", va="center", fontsize=7,
                    color=color, fontweight=weight)

    ax.set_xticks(range(len(pv_b.columns))); ax.set_xticklabels(pv_b.columns.astype(int))
    ax.set_yticks(range(len(pv_b.index)));   ax.set_yticklabels(pv_b.index.astype(int))
    ax.set_xlabel(f"Unique IGK {kind} per cell bin", fontsize=12)
    ax.set_ylabel(f"Unique IGH {kind} per cell bin", fontsize=12)
    ax.set_title(f"Rep2 cellbin: dominant-B by (IGH, IGK) {kind}\n"
                 f"(overall B = {overall_b:.1%})", fontsize=13)
    cb = plt.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("Fraction with dominant = B cell", fontsize=11)
    for i in range(pv_b.shape[0]):
        for j in range(pv_b.shape[1]):
            v = pv_b.values[i,j]
            if not np.isnan(v) and TARGET_B_LOW<=v<=TARGET_B_HIGH:
                ax.add_patch(plt.Rectangle((j-0.5,i-0.5),1,1,
                              fill=False, edgecolor="gold", linewidth=2.5))
    plt.tight_layout()
    p = os.path.join(outdir, f"IGH_IGK_{kind}_heatmap_cellbin.png")
    plt.savefig(p, dpi=200, bbox_inches="tight"); plt.close()
    print("  saved:", p)

# 整体 B 比例
matched_all = joint[joint["spot_id"].notna()]
overall_b = matched_all["is_B_dominant"].sum()/len(matched_all) if len(matched_all) else np.nan
print(f"  overall dominant-B = {overall_b:.1%}")

# targetSeq 版
sw_seq = sweep(joint, "n_IGH_seq", "n_IGK_seq",
                SEQ_AXIS_MAX_IGH, SEQ_AXIS_MAX_IGK)
sw_seq.to_csv(os.path.join(OUT_DIR,"IGH_IGK_targetSeq_sweep.csv"), index=False)
plot_heatmap(sw_seq, "targetSeq", OUT_DIR, overall_b)

# coords 版
sw_co = sweep(joint, "n_IGH_coord", "n_IGK_coord",
              COORD_AXIS_MAX_IGH, COORD_AXIS_MAX_IGK)
sw_co.to_csv(os.path.join(OUT_DIR,"IGH_IGK_coords_sweep.csv"), index=False)
plot_heatmap(sw_co, "coords", OUT_DIR, overall_b)

print("Step 4 done.")