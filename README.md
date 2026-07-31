# 🧬 Single-Cell RNA Sequencing Analysis of Human PBMCs

This repository contains an end-to-end bioinformatics pipeline in **R** for the quality control, processing, and downstream analysis of **Human (*Homo sapiens*) Peripheral Blood Mononuclear Cells (PBMCs)** using single-cell RNA sequencing (scRNA-seq) data.

Additionally, this repository includes a comprehensive slide deck (`.pdf`) covering both the single-cell workflow and an integrated **Bulk RNA-seq analysis**.

---

## 📌 Overview

Single-cell RNA sequencing allows for the characterization of transcriptomic heterogeneity at individual cell resolution. This workflow covers the standard Seurat-based processing pipeline to identify distinct immune cell populations from human PBMC samples.

### Key Analysis Steps
1. **Quality Control (QC) & Filtering:** Removal of low-quality cells, doublets, and damaged cells based on mitochondrial gene expression (`percent.mt`), UMIs, and detected gene counts.
2. **Normalization:** Log-normalization.
3. **Dimensionality Reduction:** Linear (PCA) and non-linear (UMAP/t-SNE) dimensionality reduction.
4. **Clustering & Data Visualization:** Unsupervised graph-based clustering (Louvain algorithm) and data visualization via Heatmap, Feature plot and Dot plot
5. **Annotation:** Cell-type annotation using canonical cell markers.
