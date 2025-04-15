import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/marker_plot/find_plot/Myeloid_Mono");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Myeloid.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["S100A8","S100A9","FCN1","LST1","FCGR3A","LILRB2","CD68","CD14","CD300E","CLEC12A","LYZ","CD36","CYP1B1","SELL","CDKN1C","C5AR1","SERPINA1","SMIM25","EREG","SLC11A1","LTA4H",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_Mono",
                  save="Myeloid_Mono.markergene.0.raw.pdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_Mono.markergene.0.xxx.pdf")
