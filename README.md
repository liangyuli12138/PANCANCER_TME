# Pan-Cancer Project Analysis Readme

This repository documents the key datasets, analysis workflows, scripts, and output files used in the pan-cancer single-cell and spatial transcriptomics project.

All paths below are shown as **relative project paths** and start directly from the numbered analysis directory.

---

## Contents

- [1. Core Data Resources](#1-core-data-resources)
- [2. Single-Cell Analysis](#2-single-cell-analysis)
  - [2.1 Data Integration](#21-data-integration)
  - [2.2 Doublet Removal](#22-doublet-removal)
  - [2.3 Primary Clustering](#23-primary-clustering)
  - [2.4 Cell-Type Annotation and Differential Expression](#24-cell-type-annotation-and-differential-expression)
  - [2.5 Secondary Clustering by Major Cell Type](#25-secondary-clustering-by-major-cell-type)
  - [2.6 Malignant Epithelial Cell Identification by CNV](#26-malignant-epithelial-cell-identification-by-cnv)
  - [2.7 Malignant Cell State Identification by NMF](#27-malignant-cell-state-identification-by-nmf)
  - [2.8 ICM Analysis](#28-icm-analysis)
- [3. Spatial Transcriptomics Analysis](#3-spatial-transcriptomics-analysis)
  - [3.1 Raw Spatial Data and Cell Segmentation](#31-raw-spatial-data-and-cell-segmentation)
  - [3.2 CellBin Quality Filtering](#32-cellbin-quality-filtering)
  - [3.3 CellBin Expansion and Correction](#33-cellbin-expansion-and-correction)
  - [3.4 Cell2location Annotation](#34-cell2location-annotation)
  - [3.5 UCell Signature Scoring](#35-ucell-signature-scoring)
  - [3.6 Distance to Tumor Boundary](#36-distance-to-tumor-boundary)
  - [3.7 ICM Spatial Distance Analysis](#37-icm-spatial-distance-analysis)
  - [3.8 Cell-Type Distance Analysis in LYM123](#38-cell-type-distance-analysis-in-lym123)
  - [3.9 Cell Neighborhood Detection](#39-cell-neighborhood-detection)
  - [3.10 Cell Neighborhood Clustering](#310-cell-neighborhood-clustering)
  - [3.11 Figure 5 and Figure 6 Analyses](#311-figure-5-and-figure-6-analyses)
  - [3.12 ICAR Contour Expansion Analysis](#312-icar-contour-expansion-analysis)
- [4. Key Notes](#4-key-notes)

---

# 1. Core Data Resources

**All Data havehave been deposited in the OMIX, China National Center for Bioinformation / Beijing Institute of Genomics, Chinese Academy of Sciences (https://ngdc.cncb.ac.cn/omix: accession no. OMIX006904).**
## Single-cell data

### Single-cell expression matrix list

```text
18.upload/SC.matrix.list
```

### Single-cell matrices compressed by cancer type

```text
18.upload/SC.matrix.tar.gz.file
```

### Final clustered and annotated single-cell objects

```text
09.cluster_show/01.all.cluster/oridata_allcancer/pancancer.ref.raw.0917.h5ad
```

```text
09.cluster_show/01.all.cluster/oridata_allcancer/h5ad2rds/pancancer.ref.raw.0917.rds
```

## Spatial data

### Spatial tissue GEM list

```text
18.upload/SC.tissue.gem.list
```

### CellBin segmentation and Cell2location annotation data

```text
18.upload/ST.cellbin.h5ad.list
```

```text
18.upload/ST.cellbin.list
```

---

# 2. Single-Cell Analysis

## 2.1 Data Integration

All single-cell samples were first merged into a unified pan-cancer object.

### Script

```text
05.merge_ori_maxtir/01.merge_ori_matrix/join.ori.matrix.py
```

### Integrated dataset

```text
05.merge_ori_maxtir/01.merge_ori_matrix/pancancer.concatenate.oridata.h5ad
```

---

## 2.2 Doublet Removal

Doublets were identified and removed before downstream clustering.

### Analysis directory

```text
02.DoubletFind
```

### Script

```text
02.DoubletFind/DoubletFind.R
```

---

## 2.3 Primary Clustering

Primary clustering was performed on the integrated single-cell dataset.

### Initial clustering and filtering

Cells with high expression of neural-related genes were filtered during the initial clustering step.

```text
03.cluster_primary/cluster.py
```

### Final clustering used in the manuscript

After manually removing clusters that could not be confidently annotated, the following script was used to generate the final clustering result shown in the manuscript.

```text
09.cluster_show/01.all.cluster/cluster.py
```

### Differential expression between major cell populations

```text
08.cluster_show/02.diff_gene/findmarker.py
```

---

## 2.4 Cell-Type Annotation and Differential Expression

### Marker dot plot

```text
09.cluster_show/04.marker_plot/plot.py
```

### h5ad objects for individual major cell populations

```text
09.cluster_show/02.split
```

### Differential expression analysis for individual cell populations

```text
09.cluster_show/03.diff
```

These outputs were used for further annotation and downstream subclustering of major cell populations.

---

## 2.5 Secondary Clustering by Major Cell Type

The following directories contain the original and manually revised secondary clustering results.

| Cell type | Original clustering | Revised clustering |
|---|---|---|
| B cells | `07.cluster_secondary/B/1000_15_0.5/` | `05.cluster_secondary_filter/07.cluster_new/B` |
| CD4 T cells | `05.cluster_secondary_filter/05.cluster/CD4_pca_filter/CD4_pca/2000_20_0.6/` | `05.cluster_secondary_filter/07.cluster_new/CD4` |
| CD8 T cells | `05.cluster_secondary_filter/05.cluster/CD8_pca/CD8_pca/2000_50_1.5/` | `05.cluster_secondary_filter/07.cluster_new/CD8` |
| Endothelial cells | `05.cluster_secondary_filter/05.cluster/EC_pca/EC_pca/2000_30_1.5/` | `05.cluster_secondary_filter/07.cluster_new/EC` |
| Fibroblasts | `07.cluster_secondary/Fibroblast/more/Fibroblast/1500_8_0.5_0.5/` | `05.cluster_secondary_filter/07.cluster_new/Fibroblast` |
| Mural cells | `05.cluster_secondary/Mural_cell/2000_50_0.3/` | `05.cluster_secondary_filter/07.cluster_new/Mural_cell` |
| Myeloid cells | `07.cluster_secondary/Myeloid/1500_25_1` | `05.cluster_secondary_filter/07.cluster_new/Myeloid` |
| NK cells | `05.cluster_secondary_filter/05.cluster/NK/2000_50_1.5` | `05.cluster_secondary_filter/07.cluster_new/NK` |

> Epithelial cells were not annotated through the same secondary clustering workflow. Malignant epithelial cells were identified using **CNV analysis followed by NMF-based state analysis**.

---

## 2.6 Malignant Epithelial Cell Identification by CNV

CNV analysis was performed separately for the tumor sample of each patient.

### Analysis directory

```text
08.filter_by_gene/09.cancer_cnv/05.find_cancer_donor
```

### Script

```text
run_infercnvpy.py
```

Malignant clusters were manually identified based on CNV scores shown on the UMAP embedding.

### Extracted malignant epithelial cells

High-confidence malignant populations:

```text
08.filter_by_gene/09.cancer_cnv/11.cnmf/00.data_core
```

Secondary malignant populations:

```text
08.filter_by_gene/09.cancer_cnv/11.cnmf/01.data_sub
```

The directory naming convention is:

- `core`: clusters with clear CNV score differences
- `sub`: clusters with weaker or secondary CNV score differences

Both groups were retained for downstream NMF analysis.

---

## 2.7 Malignant Cell State Identification by NMF

NMF was performed independently for malignant epithelial cells from each sample.

### NMF analysis for core malignant populations

```text
08.filter_by_gene/09.cancer_cnv/11.cnmf/02.cnmf.core
```

### NMF analysis for secondary malignant populations

```text
08.filter_by_gene/09.cancer_cnv/11.cnmf/03.cnmf.sub
```

### Merge NMF results across samples

```text
08.filter_by_gene/09.cancer_cnv/11.cnmf/04.merge
```

### Identify consensus NMF modules

```text
08.filter_by_gene/09.cancer_cnv/11.cnmf/04.merge/hypergeometric_consensus
```

### Final sample-by-module score matrix

```text
08.filter_by_gene/09.cancer_cnv/11.cnmf/04.merge/hypergeometric_consensus/final_modules_monify/final.0616.malig.gene.modules.xls
```

### Signature genes for each NMF module

```text
08.filter_by_gene/09.cancer_cnv/11.cnmf/04.merge/hypergeometric_consensus/final_modules_monify/modules_filter_gene.csv
```

### Cell-level module scoring

Each malignant epithelial cell was scored against the final modules to assign malignant cell states.

```text
08.filter_by_gene/09.cancer_cnv/11.cnmf/05.cell.scoring
```

---

## 2.8 ICM Analysis

### Normalized proportions of single-cell populations

```text
10.co_existence/pancancer.ref.0905.final.obs.csv.from.windows.zscore.xls
```

### Pearson correlation analysis and clustering

```text
10.co_existence/filter_20240317
```

### Differentially expressed genes for the three immune groups

```text
10.co_existence/marker_immune
```

---

# 3. Spatial Transcriptomics Analysis

## 3.1 Raw Spatial Data and Cell Segmentation

### CellBin segmentation TIFF downloaded from SAP

```text
12.spa/04.mask.tif.from.qingdao/mask.tif
```

### Registered ssDNA TIFF files

```text
12.spa/05.cell_correct/regist_tif/
```

---

## 3.2 CellBin Quality Filtering

Potential false-positive CellBins were filtered based on ssDNA signal intensity.

```text
13.spa/04.filter_bright/02.filter
```

---

## 3.3 CellBin Expansion and Correction

The original CellBin segmentation was expanded and corrected before downstream analysis.

### Expansion/correction scripts

```text
13.spa/05.cell_correct/shell
```

### Final corrected h5ad files

```text
13.spa/05.cell_correct/result/*/*.final.h5ad
```

---

## 3.4 Cell2location Annotation

Cell2location analysis was performed on a GPU-enabled computing environment in Hainan.

### Analysis directory

```text
12.spa/10.cell2location/run
```

### Main script

```text
cel2loc.py
```

For each spatial chip, CellBins were randomly divided into batches of approximately **10,000 cells** for computation.

---

## 3.5 UCell Signature Scoring

UCell scoring was performed for each spatial cell.

The scoring included:

- four tumor microenvironment subtypes
- Hallmark gene signatures
- additional Reactome terms

### Tumor microenvironment subtype signatures

```text
16.microenvironment/02.subtypes.all
```

### Hallmark signatures

```text
16.microenvironment/04.hallmarker.all
```

### Additional Reactome terms

```text
16.microenvironment/06.fig2
```

Results stored under `scale` directories were Z-score normalized across samples from different cancer types for each signature.

---

## 3.6 Distance to Tumor Boundary

Spatial distances between individual cells and the tumor boundary were calculated.

### Analysis directory

```text
13.spa/14.icm/02.border
```

### Distance calculation script

```text
13.spa/14.icm/02.border/get.dist.pl
```

### Distance definition

Each cell was first classified as either inside or outside the tumor region.

- For a cell inside the tumor, the distance to the nearest cell outside the tumor was used as its distance to the tumor boundary.
- For a cell outside the tumor, the distance to the nearest cell inside the tumor was used as its distance to the tumor boundary.

### Plotting script

```text
13.spa/14.icm/02.border/plot/all.plot.py
```

---

## 3.7 ICM Spatial Distance Analysis

Distances between ICM cells belonging to the same or different ICM classes were analyzed.

### Analysis directories

```text
13.spa/14.icm/03.distance
```

```text
13.spa/14.icm/03.distance/merge
```

---

## 3.8 Cell-Type Distance Analysis in LYM123

Spatial distances between different cell types were calculated within LYM123.

### Analysis directory

```text
13.spa/14.icm/05.distance
```

### Script

```text
13.spa/14.icm/05.distance/get.all.dis.pl
```

For each cell, the mean distance to the **five nearest cells** of either the same or a different cell type was calculated.

### Output files

```text
13.spa/14.icm/05.distance/Lymphoid_1_2.distance.stat.csv
```

```text
13.spa/14.icm/05.distance/Lymphoid3.distance.stat.csv
```

---

## 3.9 Cell Neighborhood Detection

### Neighborhood detection script

```text
13.spa/09.cluster_immune/03.net/network.py
```

### Neighborhood detection output

```text
13.spa/09.cluster_immune/03.net/out
```

### Neighborhood visualization

```text
13.spa/10.cluster_immune/04.plot
```

---

## 3.10 Cell Neighborhood Clustering

### Clustering script

```text
13.spa/10.cluster_immune/06.more_celltype_cluster/03.cluster_test/cluster.py
```

The clustering radius used in the final analysis was:

```text
r = 1.5
```

Some clusters were manually merged after the initial clustering.

### Comparison between lymphoid and myeloid neighborhoods

```text
13.spa/10.cluster_immune/06.more_celltype_cluster/05.l2m
```

---

## 3.11 Figure 5 and Figure 6 Analyses

### Figure 5-related analyses

```text
13.spa/10.cluster_immune/07.fig5.stage
```

### Figure 6-related analyses

```text
13.spa/15.fig6
```

---

## 3.12 ICAR Contour Expansion Analysis

ICAR regions were expanded by **100, 300, and 500 bins** for downstream analysis.

### Expansion script

```text
13.spa/18.contour/zhangzhao.py
```

### Downstream analysis directory

```text
13.spa/19.contour.after
```

---

# 4. Key Notes

1. Malignant epithelial cells were identified using a **CNV + NMF** workflow rather than conventional secondary clustering.
2. CNV analysis was performed independently for each patient's tumor sample.
3. Malignant clusters were manually selected based on CNV scores visualized on UMAP.
4. Both `core` and `sub` malignant populations were included in downstream NMF analysis.
5. Cell2location was run in a GPU-enabled environment.
6. UCell results under `scale` directories were Z-score normalized across cancer types.
