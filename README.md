# 🧬 STOmics BCR × Cell Type Analysis Pipeline

> A modular pipeline for **Stereo-seq** spatial transcriptomics data, combining BCR (IGH ∪ IGK ∪ IGL) mapping with **RCTD** cell-type deconvolution to quantify the relationship between BCR signal density and B-cell enrichment across three spatial-binning strategies.

[![Python](https://img.shields.io/badge/Python-3.9+-blue.svg)](https://www.python.org/)
[![Status](https://img.shields.io/badge/Status-Production-success.svg)]()
[![License](https://img.shields.io/badge/License-Internal-lightgrey.svg)]()

---

## ✨ Highlights

- 📊 **Three spatial-binning strategies** compared in one go: `bin20` / `cell bin` / `nucleus-strict`
- 🎯 **Replicate-agnostic** — change paths once, pipeline runs end-to-end
- 🖼️ **Publication-ready figures** with consistent black-background styling
- 🔁 **Reproducible** — fixed random seeds, deterministic outputs
- 💾 **Cached intermediates** — re-run downstream steps without recomputing heavy intersections

---

## 📚 Table of Contents

- [Pipeline Overview](#-pipeline-overview)
- [Required Inputs](#-required-inputs)
- [Configuration](#️-configuration)
- [Usage](#️-usage)
- [Outputs](#-outputs)
- [Methodology](#-methodology)
- [Visual Style Guide](#-visual-style-guide)
- [Implementation Notes](#-implementation-notes)
- [Dependencies](#-dependencies)
- [Adapting to a New Replicate](#-adapting-to-a-new-replicate)

---

## 🗂 Pipeline Overview

| # | Script | Purpose |
|:-:|---|---|
| 🚀 | `run_all.py` | One-click entry — runs steps 1 → 4 in order |
| 1️⃣ | `step1_three_QC.py` | BCR QC for **cellbin / bin20 / nucleus-strict** (no IGH/K/L distinction) → per-cell stats + B-proportion line plot + 5 dominant-celltype bar plots |
| 2️⃣ | `step2_integrated_line.py` | Combines the three line plots from step 1 into **one comparison figure** |
| 3️⃣ | `step3_ROI_visuals.py` | (a) 6 random ROIs for cell segmentation QC<br>(b) Center ROI with BCR overlay (DAPI + cellbin / nucleus)<br>(c) BCR in/out judgement under three criteria |
| 4️⃣ | `step4_IGHK_heatmap.py` | **IGH × IGK** heatmaps on cellbin (targetSeq + coords versions), colored by dominant-B fraction |

---

## 📥 Required Inputs

For each replicate, prepare the following files:

| 📄 File | Format | Purpose |
|---|---|---|
| `<sample>.cellbin_*.h5ad` | AnnData | Cell polygons (`obsm['cell_border']`) + centroids |
| `<sample>.bin20_*.h5ad` | AnnData | bin20 grid (used to detect bin size + valid bins) |
| `<sample>_DAPI_regist.tif` | TIFF | Registered DAPI image (background) |
| `<sample>_DAPI_mask.tif` | TIFF | Nucleus segmentation mask (binary or labeled) |
| `<sample>_read1.fq.gz.meta.gz` | gz CSV | XCR meta with `topChains`, `targetSequences`, `x`, `y` |
| `RCTD_celltype_proportion_wide.csv` | CSV | RCTD output (wide format: one column per celltype) — generate via the companion R script |

> **📝 RCTD CSV format** must contain the following columns:
> ```
> spot_id, x, y, <celltype_1>, <celltype_2>, ..., <celltype_N>
> ```
> Each celltype column = **proportion** of that celltype at the spot.
> Use the companion R script to export this from a `spacexr` RCTD `.rds` file.

---

## ⚙️ Configuration

Each script has a **path-config block** at the top. To switch replicates, edit:

```python
H5AD_CELLBIN = "/.../<sample>.cellbin_1.0.adjusted.h5ad"
H5AD_BIN20   = "/.../<sample>.bin20_1.0.h5ad"
REGIST_PATH  = "/.../<sample>_DAPI_regist.tif"
MASK_PATH    = "/.../<sample>_DAPI_mask.tif"
META_PATH    = "/.../<sample>_read1.fq.gz.meta.gz"
RCTD_WIDE    = "/.../RCTD_celltype_proportion_wide.csv"
OUT_DIR      = "/.../<replicate_root>"
```

The `B_RELATED` list at the top defines which RCTD columns are merged into a single `"B cell"` category. **Add or remove entries** if your RCTD reference has different B-cell subtype names.

---

## ▶️ Usage

### 🟢 Option A — Run everything

```bash
python run_all.py
```

This runs steps 1 → 4 sequentially. If any step fails, the pipeline stops.

### 🟡 Option B — Run individual steps

```bash
python step1_three_QC.py
python step2_integrated_line.py
python step3_ROI_visuals.py
python step4_IGHK_heatmap.py
```

> ⏱️ **Step 1 is the slowest** (per-cell polygon × point intersection over the whole tissue). Steps 2–4 reuse step 1's CSVs where applicable.

### 🎛 Optional flags

| Scenario | What to do |
|---|---|
| 🚫 No nucleus mask | Leave `MASK_PATH = ""` in step 1 → nucleus-strict analysis is skipped automatically |
| 🚫 No DAPI image | Leave `REGIST_PATH = ""` in step 3 → backgrounds fall back to cell centroids (other panels still work) |

---

## 📤 Outputs

```
<OUT_DIR>/
│
├── 📁 cellbin_BCR_QC/                 ← step 1: cellbin
│   ├── BCR_combined_cellbin_stats.csv
│   ├── B_proportion_by_coords_line.{png,csv}
│   └── dominant_composition_bin{0,1,2,3,ge4}_BCRcoords.png
│
├── 📁 bin20_BCR_QC/                   ← step 1: bin20
│   └── (same structure)
│
├── 📁 nucleus_BCR_QC/                 ← step 1: nucleus-strict (only if MASK_PATH given)
│   └── (same structure)
│
├── 📁 Integrated_line_plot/           ← step 2
│   └── integrated_B_proportion.{png,pdf}
│
├── 📁 cell_segmentation_QC/           ← step 3a
│   ├── 0_overview_6_random_ROIs.png
│   └── ROI_<i>_x*_y*_cellseg_QC.png   (×6)
│
├── 📁 BCR_overlay_center/             ← step 3b
│   ├── 0_overview_center_ROI.png
│   ├── 1_DAPI_cellbin_BCR.png
│   ├── 2_DAPI_nucleus_BCR.png
│   └── 3_DAPI_nucleus_cellbin_BCR.png
│
├── 📁 BCR_inout_judgement/            ← step 3c
│   ├── 1_bin20_DAPI_BCR_inout.png
│   ├── 2_cellbin_DAPI_BCR_inout.png
│   └── 3_nucleus_DAPI_BCR_inout.png
│
└── 📁 IGH_IGK_heatmap_cellbin/        ← step 4
    ├── IGH_IGK_targetSeq_heatmap_cellbin.png
    ├── IGH_IGK_coords_heatmap_cellbin.png
    └── IGH_IGK_*_sweep.csv
```

---

## 🧪 Methodology

### 🩸 BCR Coordinates

> **BCR coordinates** = unique `(x, y)` records where `topChains ∈ {IGH, IGK, IGL}`.

By default **IGH / IGK / IGL are merged** ("any BCR"). Step 4 separates them only for the dual-chain heatmap.

### 📐 Three Spatial-Binning Criteria

Used throughout the pipeline for in/out judgement:

| Criterion | Definition |
|---|---|
| 🟦 **`bin20`** | BCR coordinate falls inside **any valid bin** in the bin20 AnnData |
| 🟢 **`cell bin`** | BCR coordinate falls inside **any cell polygon** (`obsm['cell_border']`) |
| 🔴 **`nucleus-strict`** | BCR coordinate falls inside **both a cell polygon AND nucleus mask foreground** (`mask > 0`) |

### 🎯 Dominant-B Definition

RCTD per-spot proportions are merged across all `B_RELATED` columns into a single `"B cell"` category:

```
dominant_celltype = argmax(proportions)
```

A spot is **B-dominant** if its `dominant_celltype == "B cell"`.

---

## 🎨 Visual Style Guide

| Element | Color | Notes |
|---|---|---|
| 🖤 Background | Black | All single-tissue images use black for inline-DAPI compatibility |
| 🟦 Boundaries / grids | Cyan `#00E5FF` | Avoids clashing with red/pink cell-status palettes |
| 💗 BCR — in-boundary | Magenta `#FF00FF` | Filled circle marker |
| ❌ BCR — out-of-boundary | Red `×` `#FF1744` | "x" marker for emphasis |
| 🥇 Heatmap highlight | Gold border | Cells in 80–90% B-fraction range |

---

## 📝 Implementation Notes

- 🔍 **bin20 coordinate convention** (corner vs. center) is **auto-detected** from the AnnData by checking whether `bx` values are integer-multiples of `bin_size`.
- 💾 **Step 1 caches** per-cell stats to CSV, so steps 2–4 can be re-run quickly without recomputing the polygon × point intersection.
- 🖨️ All figures use `bbox_inches="tight"` and `facecolor="black"` for **clean PPT pasting** (no white bleed).

---

## 📦 Dependencies

### 🐍 Python (≥ 3.9)

```
numpy            pandas         scipy
matplotlib       scikit-image
anndata          scanpy
tifffile
```

### 📊 R (only for generating `RCTD_celltype_proportion_wide.csv`)

```
spacexr          dplyr          tidyr
readr            tibble         Matrix
```

---

## 🔁 Adapting to a New Replicate

> ✅ **No logic changes required** — just paths.

1. 📁 Copy the entire script folder to a new working location (or just re-edit in place)
2. ✏️ Update the **path block** at the top of each `step*.py` to point to the new replicate's files
3. 🧬 If RCTD column names differ, update `B_RELATED` accordingly
4. 🚀 Run `python run_all.py`

That's it — same pipeline, same logic, new data.

---

<div align="center">

**🧬 Built for reproducible spatial BCR analysis 🧬**

</div>
