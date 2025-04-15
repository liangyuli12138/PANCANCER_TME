import scanpy as sc
import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

adata = sc.read("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/06.ex/leiden/immune.cluster.r0.5.h5ad")
mg=["EC_Angiogenic","EC_Artery","EC_Lymph","EC_Vein","Epithelium_Malig_Basal","Epithelium_Malig_Cycle","Epithelium_Malig_Glandular","Epithelium_Malig_Interferon","Epithelium_Malig_Migration","Epithelium_Malig_Stress","Epithelium_Normal","Fibroblast_apCAF_CD74","Fibroblast_iCAF_IGFBP6","Fibroblast_iCAF_KCNN3","Fibroblast_mCAF_POSTN","Fibroblast_mCAF_WNT5A","Lymphoid_Plamsa_IGKC","Lymphoid_Plamsa_IGLC","Mural_cell_Pericyte1","Mural_cell_Pericyte2","Mural_cell_SMC1","Mural_cell_SMC2","Unknown"]
#marker_genes_dict = {'4','5,','0','6','1','2','3'}
sc.pl.stacked_violin(adata,mg,groupby='leiden', rotation=90, swap_axes=False, 
#categories_order=['3','2','1','6','0','5','4'],
vmax=0.5,
colorbar_title='Median proportion in gruop',
save='.ex.group.cell.png')
