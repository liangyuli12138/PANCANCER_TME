import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/marker_plot/find_plot/Epithelium_Malig_Interferon");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Epithelium.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["B2M","STAT1","UBE2L6","PSMB8","HLA-E","IFITM3","XAF1","HLA-A","IFI6","CD74",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Malig_Interferon",
                  save="Epithelium_Malig_Interferon.markergene.0.raw.pdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Malig_Interferon.markergene.0.xxx.pdf")
