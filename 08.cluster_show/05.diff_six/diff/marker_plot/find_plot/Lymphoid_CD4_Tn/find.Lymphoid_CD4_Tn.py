import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/marker_plot/find_plot/Lymphoid_CD4_Tn");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_CD4.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["CCR7","LEF1","SELL","TCF7","RPL31","RPL21","FHIT","TCEA3","CDC25B","MAL","WDR74","PGAP1","AP3M2","TXK","CXCR5","GIMAP7","GIMAP4","IKZF3","TRAF3IP3","GIMAP1","TC2N","NUCB2",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD4_Tn",
                  save="Lymphoid_CD4_Tn.markergene.0.raw.pdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD4_Tn.markergene.0.xxx.pdf")
