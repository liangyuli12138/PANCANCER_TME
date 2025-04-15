import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/marker_plot/find_plot/Myeloid_Mast");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Myeloid.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["TPSB2","TPSAB1","KIT","CLU","VWA5A","CPA3","MS4A2","IL1RL1","HDC","HPGDS",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_Mast",
                  save="Myeloid_Mast.markergene.0.raw.pdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_Mast.markergene.0.xxx.pdf")
