import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib as mpl

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/04.marker_plot")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Fibroblast.0905.h5ad")
marker_genes_dict = {
"Fibroblast_mCAF_POSTN":["BGN","SPARC","POSTN","CTHRC1","COL1A1","COL3A1","COL8A1","VCAN"],
"Fibroblast_mCAF_KRT19":["MYLK","KRT19","TCEAL2","TCEAL4","PEG3","ESR1"],
"Fibroblast_mCAF_WNT5A":["TBX3","WNT5A","PRDM1","TRPA1","MME","COL7A1"],
"Fibroblast_iCAF_IGFBP6":["IGFBP6","SCARA5","FBLN2","MFAP5"],
"Fibroblast_iCAF_IL6":["IL6","GEM","NAMPT","SOD2","SOCS3","CXCL2","SLC2A3"],
"Fibroblast_iCAF_KCNN3":["KCNN3","FXYD6","SCN7A","NDRG2","MATN2","MBP","SPARCL1"],
"Fibroblast_apCAF_CD74":["CD74","HLA-DRB1","HLA-DRB5","HLA-DRA"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.9, categories_order=["Fibroblast_mCAF_POSTN","Fibroblast_mCAF_KRT19","Fibroblast_mCAF_WNT5A","Fibroblast_iCAF_IGFBP6","Fibroblast_iCAF_IL6","Fibroblast_iCAF_KCNN3","Fibroblast_apCAF_CD74",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Fibroblast.pdf")


