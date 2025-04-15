import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib as mpl

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/04.marker_plot")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Myeloid_Mono.0905.h5ad")
marker_genes_dict = {
"Myeloid_Mono":["S100A8","S100A9","SLC11A1","SERPINA1","FCN1"],
"Myeloid_cDC1":["CLEC9A","XCR1","WDFY4","SNX3","CPNE3"],
"Myeloid_cDC2":["CD1C","FCER1A","CLEC10A","JAML","LST1"],
"Myeloid_cDC3":["LAMP3","CCR7","FSCN1","CCL19","CCL17"],
"Myeloid_pDC":["TCF4","GZMB","LILRA4","PLAC8","IRF7"],
"Myeloid_Mast":["TPSAB1","TPSB2","MS4A2","CPA3","VWA5A","KIT"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.9, categories_order=["Myeloid_Mono","Myeloid_cDC1","Myeloid_cDC2","Myeloid_cDC3","Myeloid_pDC","Myeloid_Mast",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Myeloid_Mono.pdf")


