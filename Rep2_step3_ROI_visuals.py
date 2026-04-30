"""
Rep2 Step 3:
  1) 6 个随机 ROI 的 cell segmentation QC 图 + 全局位置标注
  2) 中心 1 个 ROI 的 BCR 叠加图 (DAPI+cellbin / DAPI+nucleus)
  3) BCR in/out 三口径判断图
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
import tifffile
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.path import Path as MplPath
from skimage.segmentation import find_boundaries

# ====== Rep2 路径 (和 step1 保持一致) ======
H5AD_CELLBIN = "/data/scratch/projects/punim1236/spleen_stomics_2/processed/STOmics_mice_spleen_2_brighter_morecores_omejustliker1/outs/visualization/C04278A3.cellbin_1.0.adjusted.h5ad"
H5AD_BIN20   = "/data/scratch/projects/punim1236/spleen_stomics_2/processed/STOmics_mice_spleen_2_brighter_morecores_omejustliker1/outs/visualization/C04278A3.bin20_1.0.h5ad"
META_PATH    = "/data/scratch/projects/punim1236/AAA_chen_here/AAA克隆_XCR/A3_read1.fq.gz.meta.gz"
REGIST_PATH  = "/data/scratch/projects/punim1236/spleen_stomics_2/processed/STOmics_mice_spleen_2_brighter_morecores_omejustliker1/outs/image/C04278A3_DAPI_regist.tif"   # ← 如果 Rep2 有 DAPI registered tif 就填这里
MASK_PATH    = "/data/scratch/projects/punim1236/spleen_stomics_2/processed/STOmics_mice_spleen_2_brighter_morecores_omejustliker1/outs/image/C04278A3_DAPI_mask.tif"   # ← 如果 Rep2 有 nucleus mask 就填这里

REP2_DIR = "/data/scratch/projects/punim1236/AAA_chen_here/AAA克隆_XCR/Rep2"
OUT_SEG    = os.path.join(REP2_DIR, "cell_segmentation_QC")
OUT_BCR    = os.path.join(REP2_DIR, "BCR_overlay_center")
OUT_INOUT  = os.path.join(REP2_DIR, "BCR_inout_judgement")
for d in (OUT_SEG, OUT_BCR, OUT_INOUT):
    os.makedirs(d, exist_ok=True)

ROI_HALF = 200          # 半边长 (总 ROI = 400×400)
N_RANDOM = 6
RNG_SEED = 42
EPS = 0.01
BCR_CHAINS = {"IGH","IGK","IGL"}

# ====== 加载 ======
print("Loading cellbin ...")
adata = ad.read_h5ad(H5AD_CELLBIN)
if {"x","y"}.issubset(adata.obs.columns):
    C = adata.obs[["x","y"]].to_numpy()
else:
    C = adata.obsm["spatial"]
B = adata.obsm["cell_border"]
print(f"  cell bins: {len(C)}")

print("Loading bin20 ...")
adata_b = ad.read_h5ad(H5AD_BIN20)
if {"x","y"}.issubset(adata_b.obs.columns):
    bx = adata_b.obs["x"].to_numpy(float); by = adata_b.obs["y"].to_numpy(float)
else:
    bx = adata_b.obsm["spatial"][:,0].astype(float)
    by = adata_b.obsm["spatial"][:,1].astype(float)
ux=np.sort(np.unique(bx)); uy=np.sort(np.unique(by))
BIN_SIZE = float(round((np.min(np.diff(ux))+np.min(np.diff(uy)))/2))
fx=(bx/BIN_SIZE)-np.floor(bx/BIN_SIZE)
CONV = "corner" if np.mean((fx<0.05)|(fx>0.95)) >= np.mean(np.abs(fx-0.5)<0.05) else "center"
print(f"  bin20 bin_size={BIN_SIZE}, convention={CONV}")

# DAPI / mask 是可选的
img_show = None; H_img=W_img=None; nuc_bool=None; nuc_bound=None
if REGIST_PATH and os.path.exists(REGIST_PATH):
    print("Loading DAPI ...")
    img = tifffile.imread(REGIST_PATH)
    if img.ndim==3:
        img_show = img[...,0] if img.shape[-1] in [3,4] else img[0]
    else:
        img_show = img
    img_show = img_show.astype(np.float32)
    p1,p99 = np.percentile(img_show,[1,99])
    img_show = np.clip((img_show-p1)/(p99-p1), 0, 1)**1.8
    H_img, W_img = img_show.shape

if MASK_PATH and os.path.exists(MASK_PATH):
    print("Loading mask ...")
    m = tifffile.imread(MASK_PATH)
    if m.ndim==3:
        m = m[...,0] if m.shape[-1]<=4 else m[0]
    nuc_bool = m>0
    nuc_bound = find_boundaries(nuc_bool, mode="outer")

print("Loading BCR coords ...")
df_xcr = pd.read_csv(META_PATH, sep=",", compression="gzip")
df_xcr["topChains"] = df_xcr["topChains"].astype(str).str.strip()
df_bcr = df_xcr[df_xcr["topChains"].isin(BCR_CHAINS)].dropna(subset=["x","y"])
bcr_xy = df_bcr[["x","y"]].drop_duplicates().to_numpy(float)
print(f"  unique BCR coords: {len(bcr_xy)}")

# ====== 选 ROI: 6 个随机 + 1 个中心 ======
rng = np.random.default_rng(RNG_SEED)
# 选取 cellbin 质心分布的中心 80% 区域内, 避免靠边
qx = np.quantile(C[:,0],[0.1,0.9])
qy = np.quantile(C[:,1],[0.1,0.9])
center_x = float(np.mean(C[:,0])); center_y = float(np.mean(C[:,1]))
print(f"  cellbin center: ({center_x:.0f}, {center_y:.0f})")

# 6 个 random ROI
random_rois = []
for k in range(N_RANDOM):
    cx = rng.uniform(qx[0]+ROI_HALF, qx[1]-ROI_HALF)
    cy = rng.uniform(qy[0]+ROI_HALF, qy[1]-ROI_HALF)
    random_rois.append((int(cx), int(cy)))
print(f"  random ROIs: {random_rois}")

# 中心 ROI
center_roi = (int(center_x), int(center_y))
print(f"  center ROI:  {center_roi}")

# ====== 通用函数 ======
def make_roi_extent(cx, cy):
    return cx-ROI_HALF, cx+ROI_HALF, cy-ROI_HALF, cy+ROI_HALF   # x_min,x_max,y_min,y_max

def make_fig_black(figsize=10, dpi=400):
    fig = plt.figure(figsize=(figsize,figsize), dpi=dpi, facecolor="black")
    ax = fig.add_axes([0,0,1,1]); ax.set_facecolor("black")
    ax.set_aspect("equal","box"); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    return fig, ax

def draw_dapi_in_roi(ax, x_min,x_max,y_min,y_max):
    if img_show is None: return
    xs,xe=max(0,x_min),min(W_img,x_max); ys,ye=max(0,y_min),min(H_img,y_max)
    img_roi = img_show[ys:ye, xs:xe]
    ax.imshow(img_roi, cmap="gray",
              extent=[x_min,x_max,y_max,y_min],
              alpha=1, zorder=0, interpolation="nearest")

def draw_cell_borders_in_roi(ax, x_min,x_max,y_min,y_max,
                              color="#00E5FF", lw=1.2):
    cells = np.where(
        (C[:,0]>=x_min-100)&(C[:,0]<=x_max+100)&
        (C[:,1]>=y_min-100)&(C[:,1]<=y_max+100))[0]
    segs=[]
    for ci in cells:
        pts=B[ci]; v=~(pts==32767).any(axis=1)
        pts=pts[v].astype(np.float32)
        if pts.shape[0]<3: continue
        poly = C[ci]+pts[:,[0,1]]
        segs.append(np.vstack([poly, poly[0]]))
    if segs:
        ax.add_collection(LineCollection(segs, colors=color, linewidths=lw,
                                          antialiased=True, zorder=2))

def draw_nuc_boundary_in_roi(ax, x_min,x_max,y_min,y_max):
    if nuc_bound is None: return
    xs,xe=max(0,x_min),min(W_img,x_max); ys,ye=max(0,y_min),min(H_img,y_max)
    sub = nuc_bound[ys:ye, xs:xe]
    masked = np.ma.masked_where(~sub, sub)
    ax.imshow(masked, cmap="winter",
              extent=[x_min,x_max,y_max,y_min],
              alpha=0.95, zorder=1, interpolation="nearest")

def setup_roi_axes(ax, x_min,x_max,y_min,y_max):
    ax.set_xlim(x_min,x_max); ax.set_ylim(y_max,y_min)

# =====================================================================
# A. 6 个随机 ROI + 全局位置图
# =====================================================================
print("\n" + "="*70 + "\n  A. cell segmentation QC: 6 random ROIs\n" + "="*70)

# A1. 全局图带 6 个框
fig = plt.figure(figsize=(14,14), dpi=250, facecolor="black")
ax = fig.add_axes([0,0,1,1]); ax.set_facecolor("black")
if img_show is not None:
    ax.imshow(img_show, cmap="gray",
              extent=[0,W_img,H_img,0], alpha=1, zorder=0,
              interpolation="nearest")
else:
    # 没 DAPI 就画 cellbin 质心当背景
    ax.scatter(C[:,0], C[:,1], s=0.3, c="#444444", alpha=0.5)
for i,(cx,cy) in enumerate(random_rois, 1):
    x_min,x_max,y_min,y_max = make_roi_extent(cx,cy)
    rect = Rectangle((x_min,y_min), x_max-x_min, y_max-y_min,
                     fill=False, edgecolor="#FFEB3B", linewidth=3, zorder=5)
    ax.add_patch(rect)
    ax.text(x_min, y_min-30, f"#{i}", color="#FFEB3B", fontsize=18,
            fontweight="bold", ha="left", va="bottom",
            bbox=dict(facecolor="black", edgecolor="#FFEB3B",
                      linewidth=1.2, boxstyle="round,pad=0.3"))
ax.set_aspect("equal","box"); ax.set_xticks([]); ax.set_yticks([])
if img_show is not None:
    ax.set_xlim(0,W_img); ax.set_ylim(H_img,0)
else:
    ax.set_xlim(C[:,0].min(), C[:,0].max())
    ax.set_ylim(C[:,1].max(), C[:,1].min())
for s in ax.spines.values(): s.set_visible(False)
fig.savefig(os.path.join(OUT_SEG, "0_overview_6_random_ROIs.png"),
            dpi=250, facecolor="black", pad_inches=0)
plt.close(fig); print("  saved overview")

# A2. 每个 ROI 单独的 cell seg 图 (DAPI + cell border + nucleus)
for i,(cx,cy) in enumerate(random_rois, 1):
    x_min,x_max,y_min,y_max = make_roi_extent(cx,cy)
    fig, ax = make_fig_black()
    draw_dapi_in_roi(ax, x_min,x_max,y_min,y_max)
    draw_nuc_boundary_in_roi(ax, x_min,x_max,y_min,y_max)
    draw_cell_borders_in_roi(ax, x_min,x_max,y_min,y_max)
    setup_roi_axes(ax, x_min,x_max,y_min,y_max)
    p = os.path.join(OUT_SEG, f"ROI_{i}_x{cx}_y{cy}_cellseg_QC.png")
    fig.savefig(p, dpi=400, facecolor="black", pad_inches=0); plt.close(fig)
    print(f"  saved ROI {i}: {p}")

# =====================================================================
# B. 中心 ROI + BCR 叠加
# =====================================================================
print("\n" + "="*70 + "\n  B. BCR overlay (center ROI)\n" + "="*70)

cx,cy = center_roi
x_min,x_max,y_min,y_max = make_roi_extent(cx,cy)

# B0. 全局图: 中心 ROI 框
fig = plt.figure(figsize=(14,14), dpi=250, facecolor="black")
ax = fig.add_axes([0,0,1,1]); ax.set_facecolor("black")
if img_show is not None:
    ax.imshow(img_show, cmap="gray",
              extent=[0,W_img,H_img,0], alpha=1, zorder=0,
              interpolation="nearest")
else:
    ax.scatter(C[:,0], C[:,1], s=0.3, c="#444444", alpha=0.5)
ax.add_patch(Rectangle((x_min,y_min), x_max-x_min, y_max-y_min,
                        fill=False, edgecolor="#FFEB3B", lw=3, zorder=5))
ax.text(x_min, y_min-50, f"Center ROI\nx:[{x_min},{x_max}]\ny:[{y_min},{y_max}]",
        color="#FFEB3B", fontsize=14, fontweight="bold",
        bbox=dict(facecolor="black", edgecolor="#FFEB3B",
                  linewidth=1.5, boxstyle="round,pad=0.4"))
ax.set_aspect("equal","box"); ax.set_xticks([]); ax.set_yticks([])
if img_show is not None:
    ax.set_xlim(0,W_img); ax.set_ylim(H_img,0)
else:
    ax.set_xlim(C[:,0].min(), C[:,0].max())
    ax.set_ylim(C[:,1].max(), C[:,1].min())
for s in ax.spines.values(): s.set_visible(False)
fig.savefig(os.path.join(OUT_BCR, "0_overview_center_ROI.png"),
            dpi=250, facecolor="black", pad_inches=0); plt.close(fig)
print("  saved overview")

# 取 ROI 内 BCR
in_roi = ((bcr_xy[:,0]>=x_min)&(bcr_xy[:,0]<=x_max)&
          (bcr_xy[:,1]>=y_min)&(bcr_xy[:,1]<=y_max))
bcr_roi = bcr_xy[in_roi]
print(f"  BCR in center ROI: {len(bcr_roi)}")

def draw_bcr(ax, color="#FF00FF", s=8):
    if len(bcr_roi)==0: return
    ax.scatter(bcr_roi[:,0], bcr_roi[:,1],
               s=s, c=color, linewidths=0, alpha=1.0, zorder=10)

# B1. DAPI + cellbin + BCR
fig, ax = make_fig_black()
draw_dapi_in_roi(ax, x_min,x_max,y_min,y_max)
draw_cell_borders_in_roi(ax, x_min,x_max,y_min,y_max)
draw_bcr(ax)
setup_roi_axes(ax, x_min,x_max,y_min,y_max)
fig.savefig(os.path.join(OUT_BCR,"1_DAPI_cellbin_BCR.png"),
            dpi=400, facecolor="black", pad_inches=0); plt.close(fig)
print("  saved DAPI+cellbin+BCR")

# B2. DAPI + nucleus + BCR
fig, ax = make_fig_black()
draw_dapi_in_roi(ax, x_min,x_max,y_min,y_max)
draw_nuc_boundary_in_roi(ax, x_min,x_max,y_min,y_max)
draw_bcr(ax)
setup_roi_axes(ax, x_min,x_max,y_min,y_max)
fig.savefig(os.path.join(OUT_BCR,"2_DAPI_nucleus_BCR.png"),
            dpi=400, facecolor="black", pad_inches=0); plt.close(fig)
print("  saved DAPI+nucleus+BCR")

# B3. DAPI + nucleus + cellbin + BCR (合体)
fig, ax = make_fig_black()
draw_dapi_in_roi(ax, x_min,x_max,y_min,y_max)
draw_nuc_boundary_in_roi(ax, x_min,x_max,y_min,y_max)
draw_cell_borders_in_roi(ax, x_min,x_max,y_min,y_max)
draw_bcr(ax)
setup_roi_axes(ax, x_min,x_max,y_min,y_max)
fig.savefig(os.path.join(OUT_BCR,"3_DAPI_nucleus_cellbin_BCR.png"),
            dpi=400, facecolor="black", pad_inches=0); plt.close(fig)
print("  saved DAPI+nucleus+cellbin+BCR")


# =====================================================================
# C. BCR in/out 三口径判断 (中心 ROI)
# =====================================================================
print("\n" + "="*70 + "\n  C. BCR in/out three-criteria\n" + "="*70)

# bin20 in/out
def p2b(px,py):
    if CONV=="corner":
        return np.floor(px/BIN_SIZE)*BIN_SIZE, np.floor(py/BIN_SIZE)*BIN_SIZE
    return np.floor(px/BIN_SIZE)*BIN_SIZE+BIN_SIZE/2, np.floor(py/BIN_SIZE)*BIN_SIZE+BIN_SIZE/2

bxr, byr = p2b(bcr_roi[:,0], bcr_roi[:,1])
df_check = pd.DataFrame({"bx":bxr,"by":byr}).merge(
    pd.DataFrame({"bx":bx,"by":by, "valid":True}).drop_duplicates(),
    on=["bx","by"], how="left"
)
is_in_bin20 = df_check["valid"].fillna(False).to_numpy(bool)

# cell in/out
cells_in_roi = np.where(
    (C[:,0]>=x_min-100)&(C[:,0]<=x_max+100)&
    (C[:,1]>=y_min-100)&(C[:,1]<=y_max+100))[0]
is_in_cell = np.zeros(len(bcr_roi), bool)
for ci in cells_in_roi:
    pts=B[ci]; v=~(pts==32767).any(axis=1); pts=pts[v].astype(np.float32)
    if pts.shape[0]<3: continue
    poly = C[ci]+pts[:,[0,1]]
    bbox = ((bcr_roi[:,0]>=poly[:,0].min())&(bcr_roi[:,0]<=poly[:,0].max())&
            (bcr_roi[:,1]>=poly[:,1].min())&(bcr_roi[:,1]<=poly[:,1].max())&
            ~is_in_cell)
    if not bbox.any(): continue
    inside = MplPath(poly).contains_points(bcr_roi[bbox], radius=-EPS)
    is_in_cell[np.where(bbox)[0][inside]] = True

# nucleus in/out
is_in_nuc = np.zeros(len(bcr_roi), bool)
if nuc_bool is not None:
    pxi=np.round(bcr_roi[:,0]).astype(int); pyi=np.round(bcr_roi[:,1]).astype(int)
    inb=(pxi>=0)&(pxi<W_img)&(pyi>=0)&(pyi<H_img)
    is_in_nuc[inb] = nuc_bool[pyi[inb], pxi[inb]]

print(f"  in bin20:   {is_in_bin20.sum()}/{len(bcr_roi)}")
print(f"  in cell:    {is_in_cell.sum()}/{len(bcr_roi)}")
print(f"  in nucleus: {is_in_nuc.sum()}/{len(bcr_roi)}")

def draw_inout(ax, is_in):
    pin  = bcr_roi[is_in]; pout = bcr_roi[~is_in]
    if len(pin):
        ax.scatter(pin[:,0],pin[:,1], s=18, c="#FF00FF", marker="o",
                   linewidths=0, alpha=1, zorder=10)
    if len(pout):
        ax.scatter(pout[:,0],pout[:,1], s=75, c="#FF1744", marker="x",
                   linewidths=2.0, alpha=1, zorder=11)

# C1. bin20
fig, ax = make_fig_black()
draw_dapi_in_roi(ax, x_min,x_max,y_min,y_max)
# 画 bin20 网格
in_view = ((bx>=x_min-BIN_SIZE)&(bx<=x_max+BIN_SIZE)&
           (by>=y_min-BIN_SIZE)&(by<=y_max+BIN_SIZE))
bxv,byv = bx[in_view], by[in_view]
if CONV=="corner":
    cgx,cgy = bxv, byv
else:
    cgx,cgy = bxv-BIN_SIZE/2, byv-BIN_SIZE/2
rects = [Rectangle((cgx[i],cgy[i]), BIN_SIZE, BIN_SIZE) for i in range(len(cgx))]
ax.add_collection(PatchCollection(rects, facecolor="none",
                                   edgecolor="#00E5FF", linewidth=0.6,
                                   alpha=0.85, zorder=2))
draw_inout(ax, is_in_bin20)
setup_roi_axes(ax, x_min,x_max,y_min,y_max)
fig.savefig(os.path.join(OUT_INOUT,"1_bin20_DAPI_BCR_inout.png"),
            dpi=400, facecolor="black", pad_inches=0); plt.close(fig)

# C2. cellbin
fig, ax = make_fig_black()
draw_dapi_in_roi(ax, x_min,x_max,y_min,y_max)
draw_cell_borders_in_roi(ax, x_min,x_max,y_min,y_max)
draw_inout(ax, is_in_cell)
setup_roi_axes(ax, x_min,x_max,y_min,y_max)
fig.savefig(os.path.join(OUT_INOUT,"2_cellbin_DAPI_BCR_inout.png"),
            dpi=400, facecolor="black", pad_inches=0); plt.close(fig)

# C3. nucleus
fig, ax = make_fig_black()
draw_dapi_in_roi(ax, x_min,x_max,y_min,y_max)
draw_nuc_boundary_in_roi(ax, x_min,x_max,y_min,y_max)
draw_inout(ax, is_in_nuc)
setup_roi_axes(ax, x_min,x_max,y_min,y_max)
fig.savefig(os.path.join(OUT_INOUT,"3_nucleus_DAPI_BCR_inout.png"),
            dpi=400, facecolor="black", pad_inches=0); plt.close(fig)

print("\nStep 3 done.")