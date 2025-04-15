import scanpy as sc
import os
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/03.cluster_test/cell")
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/03.cluster_test/out/immune.cluster.11.r1.5.h5ad")

sc.pl.umap(adata_concat,color=["Lymphoid_CD4_Tn"],frameon=False, na_color="grey",save=".Lymphoid_CD4_Tn.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_CD4_Tcm"],frameon=False, na_color="grey",save=".Lymphoid_CD4_Tcm.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_CD4_Tfh"],frameon=False, na_color="grey",save=".Lymphoid_CD4_Tfh.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_CD4_Th17"],frameon=False, na_color="grey",save=".Lymphoid_CD4_Th17.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_CD4_CTL"],frameon=False, na_color="grey",save=".Lymphoid_CD4_CTL.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_CD4_Treg"],frameon=False, na_color="grey",save=".Lymphoid_CD4_Treg.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_CD4_Tstr"],frameon=False, na_color="grey",save=".Lymphoid_CD4_Tstr.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_CD8_Tn"],frameon=False, na_color="grey",save=".Lymphoid_CD8_Tn.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_CD8_Tm"],frameon=False, na_color="grey",save=".Lymphoid_CD8_Tm.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_CD8_Teff"],frameon=False, na_color="grey",save=".Lymphoid_CD8_Teff.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_CD8_Tex"],frameon=False, na_color="grey",save=".Lymphoid_CD8_Tex.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_CD8_Tstr"],frameon=False, na_color="grey",save=".Lymphoid_CD8_Tstr.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_CD8_Tisg"],frameon=False, na_color="grey",save=".Lymphoid_CD8_Tisg.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_ILC"],frameon=False, na_color="grey",save=".Lymphoid_ILC.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_MAIT"],frameon=False, na_color="grey",save=".Lymphoid_MAIT.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_NK_CD16."],frameon=False, na_color="grey",save=".Lymphoid_NK_CD16..cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_NK_CD56."],frameon=False, na_color="grey",save=".Lymphoid_NK_CD56..cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_NKT"],frameon=False, na_color="grey",save=".Lymphoid_NKT.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_B_naive"],frameon=False, na_color="grey",save=".Lymphoid_B_naive.cell.png")
sc.pl.umap(adata_concat,color=["Lymphoid_B_memory"],frameon=False, na_color="grey",save=".Lymphoid_B_memory.cell.png")
sc.pl.umap(adata_concat,color=["Myeloid_Mono"],frameon=False, na_color="grey",save=".Myeloid_Mono.cell.png")
sc.pl.umap(adata_concat,color=["Myeloid_cDC1"],frameon=False, na_color="grey",save=".Myeloid_cDC1.cell.png")
sc.pl.umap(adata_concat,color=["Myeloid_cDC2"],frameon=False, na_color="grey",save=".Myeloid_cDC2.cell.png")
sc.pl.umap(adata_concat,color=["Myeloid_cDC3"],frameon=False, na_color="grey",save=".Myeloid_cDC3.cell.png")
sc.pl.umap(adata_concat,color=["Myeloid_pDC"],frameon=False, na_color="grey",save=".Myeloid_pDC.cell.png")
sc.pl.umap(adata_concat,color=["Myeloid_Marco_C1QC"],frameon=False, na_color="grey",save=".Myeloid_Marco_C1QC.cell.png")
sc.pl.umap(adata_concat,color=["Myeloid_Marco_SPP1"],frameon=False, na_color="grey",save=".Myeloid_Marco_SPP1.cell.png")
sc.pl.umap(adata_concat,color=["Myeloid_Marco_LYVE1"],frameon=False, na_color="grey",save=".Myeloid_Marco_LYVE1.cell.png")
