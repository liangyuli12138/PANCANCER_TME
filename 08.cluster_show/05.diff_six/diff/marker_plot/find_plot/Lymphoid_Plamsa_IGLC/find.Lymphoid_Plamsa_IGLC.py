import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/marker_plot/find_plot/Lymphoid_Plamsa_IGLC");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_B.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4","MZB1","XBP1","JCHAIN","PRDM1","SDC1","IGLC2","IGLC3","IGLC7","IGLL5","ITM2C",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_Plamsa_IGLC",
                  save="Lymphoid_Plamsa_IGLC.markergene.0.raw.pdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_Plamsa_IGLC.markergene.0.xxx.pdf")
