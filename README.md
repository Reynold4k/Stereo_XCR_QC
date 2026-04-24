[README.md](https://github.com/user-attachments/files/27032211/README.md)
# Spatial BCR Quality Control Analysis

Three companion scripts for quality-controlling **B-cell receptor (BCR)** reads in **spatial transcriptomics** data (STOmics / Stereo-seq, spleen tissue). They evaluate how well BCR chain assignments (IGH / IGK / IGL) reflect true B-cell localization, at three progressively stricter spatial resolutions.

---

## What this code does

For every spatial unit (a square bin, a segmented cell, or the nuclear portion of a cell), each script:

1. **Counts BCR signal** per unit, for each chain (IGH, IGK, IGL), using two independent measures:
   - **Unique `targetSequences`** — how many distinct BCR sequences are captured
   - **Unique `(x, y)` coordinates** — how many distinct read positions fall inside
2. **Cross-references** each unit against **RCTD cell-type deconvolution** to determine whether the unit is dominated by a B-cell type.
3. **Produces:**
   - Frequency histograms of unique BCR counts per unit (per chain)
   - Spatial maps of BCR density
   - Dominant cell-type composition bar charts for each count bucket (0, 1, 2, 3, ≥4)
   - An **IGH × IGK joint heatmap** showing the fraction of units dominated by B cells at every (IGH_n, IGK_n) combination, with gold borders marking the 80–90% B-purity sweet spot

### The three scripts

| Script | Spatial unit | Constraint | Purpose |
|---|---|---|---|
| `bin20_BCR_QC_analysis.py` | Fixed **20 × 20** square bins | Coordinate must fall inside the bin | Coarse grid baseline, no cell segmentation assumptions |
| `cellbin_BCR_QC_analysis.py` | **Cell segmentation polygon** | Read must fall inside the cell border | Cell-level, permissive — includes cytoplasmic reads |
| `nucleus_strict_BCR_QC_analysis.py` | Cell polygon **AND** nucleus mask (`mask > 0`) | Read must fall inside *both* the polygon and the DAPI nucleus mask | Most stringent — nuclear reads only |

All three produce the same outputs structure, so results can be compared side-by-side.

---

## Why this QC matters

BCR chains (IGH, IGK, IGL) are biologically restricted to the **B-cell lineage**. In healthy tissue, a spatial unit densely populated by unique BCR sequences should, almost by definition, be a B cell. Any deviation is informative:

1. **Detects diffusion and contamination.** Transcripts in spatial transcriptomics can migrate during tissue permeabilization, landing in neighboring bins. If bins with high unique-IGH counts are *not* dominated by B cells in RCTD, that signals meaningful lateral diffusion or ambient-RNA noise.
2. **Validates the segmentation pipeline.** Going from `bin20` → `cellbin` → `nucleus-strict` should monotonically **increase B-cell purity** at the cost of sensitivity. If stricter constraints don't improve purity, the cell segmentation or nucleus mask may be misaligned with the transcript coordinates.
3. **Produces defensible thresholds for downstream clonotype work.** The IGH × IGK heatmap reveals concrete rules like *"units with ≥3 unique IGH and ≥4 unique IGK sequences are ≥85% dominated by B cells"*. Those thresholds become the confidence gates for clonotype assignment, repertoire diversity, and somatic-hypermutation analyses — gates that can be cited and reproduced.
4. **Surfaces technical artifacts early.** A low overall dominant-B fraction in a chain, or suspiciously high BCR counts in T-cell regions, flags problems (library prep bias, primer mis-priming, cross-chain contamination) *before* they propagate into biological conclusions.
5. **Cross-chain consistency check.** Because IGH pairs with either IGK or IGL in a given B cell, the joint IGH × IGK heatmap implicitly checks that chain co-occurrence patterns make biological sense.

In short: if these three scripts don't agree that your high-BCR units are B cells, **nothing downstream is trustworthy**.

---

## Inputs

Each script expects the following, configured in the `User configuration` block at the top of the file:

| File | Variable | Description | Used by |
|---|---|---|---|
| `*.bin20_1.0.h5ad` | `H5AD_PATH` | AnnData with bin20 coordinates in `.obs[['x','y']]` or `.obsm['spatial']` | bin20 |
| `*.cellbin_1.0.adjusted.h5ad` | `H5AD_PATH` | AnnData with cell centroids and `.obsm['cell_border']` (polygon vertices, `32767` = padding) | cellbin, nucleus-strict |
| `*.meta.gz` | `META_PATH` | Per-read BCR metadata; must contain `topChains`, `targetSequences`, `x`, `y` | all three |
| `GCsurrounding.csv` or `RCTD_celltype_proportion_wide.csv` | `RCTD_PATH` | RCTD deconvolution output. The bin20 script uses the `first_type` column; the cellbin and nucleus scripts use wide-format proportions and collapse B-related columns into a single `B cell` category. | all three |
| `*_DAPI_mask.tif` | `MASK_PATH` | DAPI nucleus segmentation mask (2-D, foreground `> 0`) | nucleus-strict only |

The list of B-related RCTD labels is hard-coded in `B_RELATED` at the top of each script — edit it if your cell-type naming differs.

---

## Installation

```bash
pip install numpy pandas matplotlib anndata tifffile
```

Python ≥ 3.8 recommended.

---

## Usage

1. Open the script and edit the `User configuration` block — paths, the output root `BASE_OUT`, and (optionally) the heatmap axis caps `SEQ_AXIS_MAX_*` / `COORD_AXIS_MAX_*`.
2. Run:

   ```bash
   python bin20_BCR_QC_analysis.py
   python cellbin_BCR_QC_analysis.py
   python nucleus_strict_BCR_QC_analysis.py
   ```

3. Inspect the console log for overall dominant-B percentages and per-file output paths.

No command-line arguments — all configuration is at the top of each file so one run is fully reproducible from the script itself.

---

## Outputs

Each script creates five subfolders under `BASE_OUT`:

```
BASE_OUT/
├── {prefix}_IGH/                    # per-chain plots + stats CSV
├── {prefix}_IGK/
├── {prefix}_IGL/
├── {prefix}_sequence_heatmap/       # IGH × IGK heatmap for unique targetSeq
└── {prefix}_coordinates_heatmap/    # IGH × IGK heatmap for unique coords
```

where `{prefix}` is `bin20`, `cellbin`, or `nucleus`.

Inside each per-chain folder you get, for both `targetSeq` and `coords`:

- `{chain}_{kind}_frequency.csv` + `.png` — how many units have N unique sequences/coordinates
- `{chain}_{kind}_spatial_bins.png` — spatial map colored by count bucket (0, 1, 2, 3, ≥4)
- `{chain}_{kind}_bin{0,1,2,3,ge4}_dominant_composition.png` — cell-type breakdown for each bucket (B subtypes merged into "B cells")
- `{chain}_{kind}_dominant_composition_all.csv` — combined table of the above
- `{chain}_{bin20|cellbin|nucleus}_stats.csv` — raw per-unit counts, joinable by `bx, by` (bin20) or `cell_idx` (cell-bin variants)

The heatmap folders contain:

- `IGH_IGK_{kind}_sweep.csv` — full sweep with `n_bins`, `n_matched`, `n_dominant_B`, `B_fraction` per (IGH_n, IGK_n)
- `IGH_IGK_{kind}_heatmap_dominant.png` — the annotated heatmap; cells with B-fraction in **[0.80, 0.90]** are gold-bordered to highlight the recommended confidence band

---

## Interpretation cheat sheet

| You see… | Likely meaning |
|---|---|
| Dominant-B fraction rises steadily with IGH_n and IGK_n | Normal — more BCR signal → cleaner B cells. Expected. |
| Dominant-B fraction is flat or noisy | Low signal, high ambient noise, or misaligned RCTD spots. Check alignment first. |
| `nucleus-strict` and `cellbin` give nearly identical results | Most BCR transcripts are nuclear / nucleus mask covers most of the cell. Fine. |
| `nucleus-strict` gives much higher B purity but very low counts | Expected for cytoplasm-biased transcripts — use `cellbin` for sensitivity, `nucleus-strict` for specificity. |
| Lots of high-IGH bins with non-B dominant type | Diffusion artifact or RCTD cell-type label mismatch — verify `B_RELATED` matches your RCTD labels. |
| IGL shows very different patterns from IGK | Biological (κ:λ ratio is tissue-dependent) or library-prep bias; compare to bulk expectations. |

---

## Notes & caveats

- **Path strings contain non-ASCII directory names** (e.g. `AAA克隆_XCR`, `扫描高占比B`). These are the author's original directory names; change them to your own layout when deploying.
- **RCTD label mismatch is the #1 source of confusion.** If the B-fractions look wrong everywhere, print `rctd["first_type"].unique()` (bin20 version) or the columns of the wide table (cellbin/nucleus versions) and make sure every B-related label is captured in `B_RELATED`.
- **Coordinate convention for bin20** is auto-detected (lower-left corner vs. bin center). Check the console log for `Bin coordinate convention: corner/center` on first run.
- **Cell-border padding** in `.obsm['cell_border']` uses the sentinel value `32767`; the code strips these before building each polygon.
- The scripts are **single-threaded** and I/O-bound on the h5ad read. A full run on one spleen sample typically finishes in a few minutes.

---

## File list

```
bin20_BCR_QC_analysis.py            # bin20-level QC
cellbin_BCR_QC_analysis.py          # cell-polygon QC
nucleus_strict_BCR_QC_analysis.py   # cell-polygon + nucleus mask QC
README.md                           # this file
```
