import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib as mpl

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/04.marker_plot")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Lymphoid_B.0905.h5ad")
marker_genes_dict = {
"Lymphoid_B_naive":["YBX3","TCL1A","IGHD","IL4R","FCER2","CD79A"],
"Lymphoid_B_memory":["GPR183","BCL2A1","CRIP1"],
"Lymphoid_Plamsa_IGLC":["IGHA1","IGHG1","MZB1","JCHAIN","IGLC2","IGLC3","IGLC7"],
"Lymphoid_Plamsa_IGKC":["IGHA1","IGHG1","MZB1","JCHAIN","IGKC"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.9, categories_order=["Lymphoid_B_naive","Lymphoid_B_memory","Lymphoid_Plamsa_IGLC","Lymphoid_Plamsa_IGKC",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Lymphoid_B.pdf")


