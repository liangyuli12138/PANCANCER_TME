import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/marker_plot/find_plot/Fibroblast_mCAF_KRT19");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Fibroblast.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["MYLK","KRT19","TCEAL2","TCEAL4","PEG3","ESR1","ST13","MDK","CSRP1","RPS4X",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Fibroblast_mCAF_KRT19",
                  save="Fibroblast_mCAF_KRT19.markergene.0.raw.pdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Fibroblast_mCAF_KRT19.markergene.0.xxx.pdf")
