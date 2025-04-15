import scanpy as sc
import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

adata = sc.read("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/05.stat/leiden3/immune.cluster.r1.h5ad")
#obslist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/05.stat/leiden3/new.obs",index_col=0)
#adata.obs = adata.obs.join(obslist)
adata.obs["leiden"] = adata.obs["leiden"].astype("category")
mg=["Lymphoid_B_naive","Lymphoid_B_memory","Lymphoid_Plamsa_IGKC","Lymphoid_Plamsa_IGLC","Lymphoid_CD4_Tn","Lymphoid_CD4_Tcm","Lymphoid_CD4_Th17",
"Lymphoid_CD4_CTL","Lymphoid_CD4_Tstr","Lymphoid_CD4_Tfh","Lymphoid_CD4_Treg","Lymphoid_CD8_Tn","Lymphoid_CD8_Tm","Lymphoid_CD8_Teff","Lymphoid_CD8_Tisg",
"Lymphoid_CD8_Tstr","Lymphoid_CD8_Tex","Lymphoid_ILC","Lymphoid_NK_CD16.","Lymphoid_NK_CD56.","Lymphoid_NKT","Myeloid_cDC1","Myeloid_cDC2","Myeloid_cDC3",
"Myeloid_Marco_LYVE1","Myeloid_Marco_C1QC","Myeloid_Marco_SPP1","Myeloid_Mono","Myeloid_pDC"]
#marker_genes_dict = {'4','5,','0','6','1','2','3'}
sc.pl.stacked_violin(adata,mg,groupby='leiden', rotation=90, swap_axes=False, 
categories_order=['9','6','0','5','3','1','7','8','4','2'],
#vmax=0.5,
#colorbar_title='Median proportion in gruop',
save='.group.cell.zscor.png')
