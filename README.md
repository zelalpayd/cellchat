# CellChat Analysis of Cell–Cell Communication in the Developing Human Brain

## Background

Brain organoid work by Yentür et al. (2025) showed higher secretome levels and higher **ADAMTS7** and **MMP15** expression at day 35 than at day 50. This project tests whether similar communication and proteolytic patterns appear in the developing human brain. It does this by running CellChat on the single-cell RNA-seq data from Braun et al. (2023).

## Data

- **Source:** CZ CELLxGENE Discover, collection *"Comprehensive cell atlas of the first-trimester developing human brain"* (Braun et al., *Science* 2023)
- **Raw data:** EGA accession `EGAD00001006049`
- **Portal filters:** Organism: *Homo sapiens*; Disease: normal; Developmental stage: Embryonic human (0–56 days); Tissue: brain
- **Full collection:** 1,665,937 cells from 26 specimens, 5–14 post-conception weeks (pcw)

## Requirements

- **Python:** `scanpy`, `seaborn`, `matplotlib`
- **R:** `Seurat`, `SeuratObject`, `zellkonverter`, `CellChat`, `patchwork`, `ggplot2`

## Workflow

Each brain region is processed and analyzed independently. In the file names below, `<brainregion>` stands for the region name.

### 1. Preprocessing

| Script | Description |
| `1_gene_expression.py` | Normalizes counts to 10,000 per cell and log-transforms them . Plots MMP15 and ADAMTS7 expression across stages to choose the time points: **9 pcw** (early), **12 pcw** (midpoint), **15 pcw** (late). |
| `2_subset_data.py` | Subsets the data to the selected stages, removes genes expressed in fewer than 10 cells, and replaces Ensembl IDs with gene symbols (`adata.var_names = adata.var["Genes"]`). Output: `<brainregion>_subset.h5ad` |
| `3_h5ad_to_rds.R` | Loads the `.h5ad` file, converts it to a Seurat object, re-normalizes it, and saves it as `.rds`. |


### 2. CellChat Analysis

This part follows the CellChat tutorial *"Comparison Analysis for Multiple Datasets using CellChat"* (S. Jin, updated 14 Feb 2025). Two signaling databases are used: **Cell–Cell Contact** and **ECM–Receptor**.

| File | Description |
|---|---|
| `4_brainpart_cellchat_interestedDB.R` | General template of the analysis. It shows the main steps for any brain region and database combination. |
| `4E_create_cell_chat_object.R` | Defines `create_cellchat()`, which builds a CellChat object for one region, stage, and signaling database (see steps below). |
| `4.1_<brainregion>_cellchat_cell_cell.ipynb` | Full analysis with the **Cell–Cell Contact** database, including plots and results. |
| `4.2_<brainregion>_cellchat_ECM.ipynb` | Full analysis with the **ECM–Receptor** database, including plots and results. |


The tutorial's comparison functions are built for two datasets, but this project has three stages. So the data is first split by pcw and a CellChat object is made for each stage. Cell types that are not present at every time point are then removed, and the objects are merged into:

- `cellchat_9_12`
- `cellchat_12_15`
- `cellchat_merged`

**Steps inside `create_cellchat()`:**

1. Remove cell groups with fewer than 10 cells and drop unused factor levels (`droplevels`).
2. Keep only the cells that appear in both the expression matrix and the metadata.
3. Remove genes that are not in the CellChat database.
4. Create the CellChat object (`createCellChat`) and attach the chosen database.
5. Preprocess the data (`subsetData`, `identifyOverExpressedGenes`, `identifyOverExpressedInteractions`) and compute communication probabilities (`computeCommunProb`).
6. Compute pathway-level signaling (`computeCommunProbPathway`) and summarize the network (`aggregateNet`).

## Usage

1. Run `1_gene_expression.py`, then `2_subset_data.py`, then `3_h5ad_to_rds.R` for each brain region.
2. To see the general workflow, read `4_brainpart_cellchat_interestedDB.R`.
3. To reproduce the results and plots, open the matching `4.1` or `4.2` notebook for the region and database you want.


## References

- Braun E. et al. (2023). Comprehensive cell atlas of the first-trimester developing human brain. *Science*.
- Yentür et al. (2025). Human dorsal forebrain organoids show differentiation-state-specific protein secretion *iScience*
- Jin S. et al. (2021). Inference and analysis of cell-cell communication using CellChat. *Nature Communications*.

