import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Myeloid_Mono");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Myeloid.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["SLC11A1","AQP9","SERPINA1","S100A9","S100A4","FLNA","S100A6","FCN1","LTA4H","UPP1","MCEMP1","TSPO","TREM1","TKT","FGR","SMIM25","ATP2B1","MARCO","NCF2","S100A8","ANPEP","FBP1","GSTO1","BCL2A1","ANXA2","STXBP2","PLIN2","FABP4","CD52","AC026369.3","RETN","EREG","SOD2","KYNU","OLR1","CD300E","S100A10","CYP27A1","CEBPB","DENND5A","EMP3","ALDH2","GLIPR2","PECAM1","OAZ1","AC007192.1","MYL6","LY6E","SDC4","OASL",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_Mono",
                  save="Myeloid_Mono.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_Mono.markergene.0.xxx.pdf")

marker_genes_dict=["CRIP1","CCRL2","PDLIM1","GNB2","PPARG","OSCAR","FTH1","SLC7A7","CD55","HCK","IL1RN","CXCL3","GMFG","CLEC12A","QSOX1","TIMP1","GPCPD1","SIRPB1","TCF7L2","SUN2","LRRFIP1","VCAN","PPIF","GABARAPL1","HBEGF","CD300LF","MIR3945HG","INHBA","HK3","GRN","CXCL2","HCAR2","PLBD1","CSTA","COTL1","NOP10","ADAM17","C5AR1","FPR1","ATP6V0D1","SNX10","PLAUR","CD44","LRPAP1","MYO1G","PILRA","EHD1","MGST3","LYZ","ADGRE5",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_Mono",
                  save="Myeloid_Mono.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_Mono.markergene.1.xxx.pdf")
