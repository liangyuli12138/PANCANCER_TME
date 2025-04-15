import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib as mpl

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/04.marker_plot")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Epithelium.0905.h5ad")
marker_genes_dict = {
"Epithelium_Malig_Migration":["MALAT1","NEAT1","NFAT5","ZSWIM6","FTX","SRRM2","SOX4","STK3"],
"Epithelium_Malig_Cycle":["TUBA1B","SNRPG","TUBB","STMN1","H2AFZ","SNRPD1","HMGN2","HMGB1"],
"Epithelium_Malig_cEMT":["TPM2","CAVIN1","COL17A1","CALD1","MYL9","ACTA2","MYLK","ACTN1","SPARCL1","TAGLN"],
"Epithelium_Malig_Interferon":["B2M","STAT1","UBE2L6","PSMB8","HLA-E","XAF1","HLA-A","IFI6","CD74"],
"Epithelium_Malig_Stress":["FOSB","JUN","ERCC1","IER2","HSPA1B","EGR1","HSP90AB1","DNAJA1","CLK1","MCL1"],
"Epithelium_Malig_Basal":["ANXA2","YWHAZ","KRT19","S100A6","KRT18"],
"Epithelium_Malig_Glandular":["TMED3","RABAC1","TTC3","TFF3","AGR2","WFDC2","CKAP4","MUC2","CKAP4"],
"Epithelium_Normal":["MYBPC1","MGP","AZGP1","TAT","SFTPC","AGR3","PHGR1","ERO1B","PIP","GFRA1"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.9, categories_order=["Epithelium_Malig_Migration","Epithelium_Malig_Cycle","Epithelium_Malig_cEMT","Epithelium_Malig_Interferon","Epithelium_Malig_Stress","Epithelium_Malig_Basal","Epithelium_Malig_Glandular","Epithelium_Normal",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Epithelium.pdf")


