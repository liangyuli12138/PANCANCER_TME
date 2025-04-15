import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/marker_plot/find_plot/Lymphoid_B_naive");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_B.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["YBX3","TCL1A","IGHD","IL4R","FCER2","SELL","CD79A","CD79B","MS4A1","CD19","HLA-A","LAPTM5","CD37","CD74","CD22","FCMR",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_B_naive",
                  save="Lymphoid_B_naive.markergene.0.raw.pdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_B_naive.markergene.0.xxx.pdf")
