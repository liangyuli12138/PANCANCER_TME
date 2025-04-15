import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib as mpl

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/04.marker_plot")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Lymphoid_CD4.0905.h5ad")
marker_genes_dict = {
"Lymphoid_CD4_Tn":["CCR7","LEF1","SELL"],
"Lymphoid_CD4_Tcm":["GPR183","CD55","YPEL5","PABPC1","EZR","IL7R"],
"Lymphoid_CD4_Tfh":["CXCL13","PDCD1","TOX","IL21","GNG4","CD200"],
"Lymphoid_CD4_Th17":["IL17A","KLRB1","CCL20","RORA","AHR","IL17F"],
"Lymphoid_CD4_CTL":["GZMA","GZMK","DNAJA1","BTG2"],
"Lymphoid_CD4_Treg":["FOXP3","IL2RA","CTLA4","TNFRSF4","IKZF2","TIGIT"],
"Lymphoid_CD4_Tstr":["FOS","JUN","TNF","FOSB","JUNB","DUSP1"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.9, categories_order=["Lymphoid_CD4_Tn","Lymphoid_CD4_Tcm","Lymphoid_CD4_Tfh","Lymphoid_CD4_Th17","Lymphoid_CD4_CTL","Lymphoid_CD4_Treg","Lymphoid_CD4_Tstr",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Lymphoid_CD4.pdf")


