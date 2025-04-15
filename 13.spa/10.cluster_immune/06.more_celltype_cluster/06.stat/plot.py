import scanpy as sc
import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

adata = sc.read("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/06.stat/immune.cluster.11.r1.5.new.h5ad")
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/06.stat/cell.merge.at",index_col=0)
adata.obs = adata.obs.join(atlist)

adata_new = sc.AnnData(obs=adata.obs)
adata_new.obsm['X_umap'] = adata.obsm['X_umap']

mg=['Lymphoid_B_naive', 'Lymphoid_B_memory', 'Lymphoid_T', 'Lymphoid_T_inhibitory',
'Myeloid_cDC', 'Myeloid_Marco_C1QC', 'Myeloid_Marco_LYVE1', 'Myeloid_Marco_SPP1']
#marker_genes_dict = {'4','5,','0','6','1','2','3'}
sc.pl.stacked_violin(adata_new,mg,groupby='new_groups', rotation=90, swap_axes=False, 
categories_order=['Myeloid8','Myeloid7','Myeloid6','Myeloid5','Myeloid4','Lymphoid3','Lymphoid2','Lymphoid1','Lymphoid0'],
vmax=0.5,
save='.group.cell.pdf')

adata_concat=adata_new
sc.pl.umap(adata_concat,color=["Lymphoid_B_naive"],frameon=False, na_color="grey",save=".Lymphoid_B_naive.cell.pdf")
sc.pl.umap(adata_concat,color=["Lymphoid_B_memory"],frameon=False, na_color="grey",save=".Lymphoid_B_memory.cell.pdf")
sc.pl.umap(adata_concat,color=["Lymphoid_T"],frameon=False, na_color="grey",save=".Lymphoid_T.cell.pdf")
sc.pl.umap(adata_concat,color=["Lymphoid_T_inhibitory"],frameon=False, na_color="grey",save=".Lymphoid_T_inhibitory.cell.pdf")
sc.pl.umap(adata_concat,color=["Myeloid_cDC"],frameon=False, na_color="grey",save=".Myeloid_cDC.cell.pdf")
sc.pl.umap(adata_concat,color=["Myeloid_Marco_C1QC"],frameon=False, na_color="grey",save=".Myeloid_Marco_C1QC.cell.pdf")
sc.pl.umap(adata_concat,color=["Myeloid_Marco_SPP1"],frameon=False, na_color="grey",save=".Myeloid_Marco_SPP1.cell.pdf")
sc.pl.umap(adata_concat,color=["Myeloid_Marco_LYVE1"],frameon=False, na_color="grey",save=".Myeloid_Marco_LYVE1.cell.pdf")

