#CellChat Analysis

Data Preparation for CellChat Analysis

This repository contains scripts for preprocessing and analyzing single-cell RNA-seq data from three brain regions across developmental stages.



Preprocessing

The preprocessing pipeline is divided into three steps:
1_gene_expression.py
Data normalization and log transformation were performed using Scanpy. Expression of MMP15 and ADAMTS7 was used to select developmental stages (9, 12, 15 pcw).
2_subset_data.py
Data was subsetted to selected stages, lowly expressed genes (less than 10 cells) were removed, and gene names were updated.
3_h5ad_to_rds.R
Data was converted from .h5ad to .rds format using zellkonverter and Seurat for downstream R analysis.



CellChat Analysis

Cell-cell communication analysis was performed using CellChat with three signaling categories:
Cell-cell contact
ECM-receptor interaction
A custom function (create_cellchat_object.R) was used to generate CellChat objects for each developmental stage. Since the dataset includes three time points, data were split by stage, processed independently, and then merged for comparison.



Repository

All preprocessing and simplified analysis scripts are included in this repository.
