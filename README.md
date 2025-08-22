scRNA-seq Analysis of CD83+ Cells in Mouse Colitis Models
Project Overview

This project analyzes publicly available single-cell RNA-seq (scRNA-seq) data from mouse models of colitis (Healthy Control, Acute Colitis, Chronic Colitis) to investigate CD83+ immune cell subsets. The aim is to identify and characterize CD83+ cells, particularly in T cell, stromal, and cycling/progenitor populations, using computational methods.

Data set source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264408


Analysis Workflow:

Data Loading and Preparation- Raw count matrices and metadata were loaded and merged into a Seurat object.

QC and Normalisation- Cells with low feature counts or high mitochondrial content were filtered; data were normalised for downstream analysis.

Clustering and Marker Identification- Seurat clustering identified 22 cell populations; presto was used to identify top markers per cluster.

CD83 Analysis- CD83 expression was quantified across cell types and conditions; statistical tests (Chi-squared and Fisher’s Exact Test) were performed to evaluate differences.

Visualisation- DimPlots, FeaturePlots, and bar plots illustrate cell type distribution and CD83+ fraction.

Reproducibility

R version: 4.5.x

Required packages: Seurat, presto, dplyr, data.table, ggplot2

To reproduce the analysis, run the scripts in order:

source("scripts/01_data_loading_preparation.R")
source("scripts/02_qc_normalization_analysis.R")

Limitations: 
-Analysis was performed on a laptop with 16GB RAM, requiring downsampling (~50%) of the full dataset. This may reduce the detection of very rare cell populations.
-The study uses publicly available data and a limited number of samples (n = 2 per condition), so results are intended to demonstrate computational analysis competency rather than to generate novel biological conclusions.
