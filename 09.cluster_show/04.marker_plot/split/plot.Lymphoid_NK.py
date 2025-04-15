import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib as mpl

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/04.marker_plot")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Lymphoid_NK.0905.h5ad")
marker_genes_dict = {
"Lymphoid_ILC":["IL4I1","LST1","KIT","FXYD5","NFKB1","AHR","NCOA7","AFF3","IL1R1"],
"Lymphoid_MAIT":["SLC4A10","TRAC","DPP4","SLAMF1","DUSP1","CD28","CITED2","TNF","EGR1"],
"Lymphoid_NK_CD16+":["FCGR3A","FGFBP2","NKG7","SPON2","ADGRG1","PRF1","GZMB","GNLY","CX3CR1","CST7"],
"Lymphoid_NK_CD56+":["IL2RB","FCER1G","CD7","IRF8","LAT2","PIK3R1","GSTP1","CD160"],
"Lymphoid_NKT":["GZMH","CD3D","CD8A","CD8B","CD52","S100A10","CD3G","PATL2","ITGB1"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.9, categories_order=["Lymphoid_ILC","Lymphoid_MAIT","Lymphoid_NK_CD16+","Lymphoid_NK_CD56+","Lymphoid_NKT",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Lymphoid_NK.pdf")


