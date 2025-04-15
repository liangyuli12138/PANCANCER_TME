import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib as mpl

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/04.marker_plot")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Lymphoid_CD8.0905.h5ad")
marker_genes_dict = {
"Lymphoid_CD8_Tn":["CCR7","TCF7","SELL"],
"Lymphoid_CD8_Tm":["MT2A","ITGA1","MT1E","KLF6"],
"Lymphoid_CD8_Teff":["CCL5","TRGC2","GZMK","CD69","DUSP2"],
"Lymphoid_CD8_Tex":["DUSP4","CTLA4","TNFRSF9","RBPJ","RGS1","PDCD1","TIGIT","CD82","CXCL13","LAG3"],
"Lymphoid_CD8_Tstr":["FOSB","ERCC1","JUN","EGR1","BRD2","IER2","JUNB","NR4A1","HEXIM1"],
"Lymphoid_CD8_Tisg":["IFIT3","IFIT2","LIPA","HERC5","IFI44","DDX58","IFIH1"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.9, categories_order=["Lymphoid_CD8_Tn","Lymphoid_CD8_Tm","Lymphoid_CD8_Teff","Lymphoid_CD8_Tex","Lymphoid_CD8_Tstr","Lymphoid_CD8_Tisg",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Lymphoid_CD8.pdf")


