import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib as mpl

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/04.marker_plot")

adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/EC_split/pancancer.split.EC.0905.h5ad")
marker_genes_dict = {
"EC_Angiogenic":["INSR","COL4A1","LAMB1","COL4A2","VWA1","MCAM"],
"EC_Artery":["SEMA3G","IGF2","AIF1L","CXCL12","GLUL","SRGN","ASS1","KCTD12","FBLN2"],
"EC_Capillary":["CA4","SLCO2A1","RGCC"],
"EC_Vein":["ACKR1","CLU","SELE","SELP","OLFM1"],
"EC_Lymph":["MMRN1","TFF3","COLEC12","CCL21","PROX1","LYVE1"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.9, categories_order=["EC_Angiogenic","EC_Artery","EC_Capillary","EC_Vein","EC_Lymph",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.EC.filter.pdf")

