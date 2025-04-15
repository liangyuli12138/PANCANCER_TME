import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Epithelium_Malig_Interferon");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Epithelium.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["B2M","STAT1","UBE2L6","PSMB8","PSMB9","HLA-E","BST2","PSME1","IFITM3","IFI27","HLA-C","XAF1","TAP1","HLA-B","HLA-A","PLSCR1","PARP9","PARP14","GUK1","PPIA","PSME2","ITM2B","APOL6","LY6E","AC040162.1","AL121594.1","ENO1","IFI44","PLAAT4","LGALS3BP","TMSB10","IFI16","TRMT112","NMI","PKM","GAPDH","RPL37","AP000350.4","SNHG6","FAU","SHISA5","GPX1","TPI1","RPS2","SUMO2","IFITM2","TAP2","GBP2","NDUFB4","RPS20",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Malig_Interferon",
                  save="Epithelium_Malig_Interferon.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Malig_Interferon.markergene.0.xxx.pdf")

marker_genes_dict=["PSMA4","SEM1","RNF213","RPS19","SNRPG","EPSTI1","OAS1","APOL3","PGK1","IFI6","TRIM22","RPS11","UBL5","ARPC3","MT-CO1","AL136295.5","RPLP0","MYL12B","S100A11","RPL19","OCIAD2","NME1-NME2","TAGLN2","IFI44L","SP100","PSMA3","CD74","TMA7","SNRPD2","RPL10","H3F3A","AC093512.2","RPL30","ATP5F1E","LDHA","NOP10","GBP4","SELENOH","OST4","RPS16","RPS7","CFL1","RPL18A","ELOB","UQCRH","AC010422.8","HNRNPA3","HLA-DRA","HINT1","ADAR",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Malig_Interferon",
                  save="Epithelium_Malig_Interferon.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Malig_Interferon.markergene.1.xxx.pdf")
