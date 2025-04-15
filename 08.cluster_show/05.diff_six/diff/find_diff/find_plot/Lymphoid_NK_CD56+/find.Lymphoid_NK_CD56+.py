import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_NK_CD56+");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_NK.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["IL2RB","TYROBP","FCER1G","CD7","LAT2","PIK3R1","IFITM2","PLCG2","CMC1","NCAM1","GSTP1","CD160","IRF8","SPRY2","SYK","STK17A","TXK","KLRC1","TRDC","KLRD1","SLFN13","CD63","KLRF1","IFITM3","KLRK1","REL","AOAH","XCL1","RAB11FIP1","ARHGAP9","HSPD1","NCR1","CLDND1","ATP8B4","MATK","CEMIP2","B3GNT7","CD38","TIPARP","PABPC1","XCL2","SH2D1B","ADGRG3","LDB2","MAP3K8","HSPE1","TIGIT","RGS1","EOMES","CCND2",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_NK_CD56+",
                  save="Lymphoid_NK_CD56+.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_NK_CD56+.markergene.0.xxx.pdf")

marker_genes_dict=["IFITM1","ZNF331","RHOQ","ICAM1","CLNK","SRGN","TTN","YPEL5","MCTP2","IFRD1","ELF1","PDCD4","MAFF","ABCB1","SPECC1","KLRC3","NR4A2","TNFSF14","TOX2","CDHR1","SH2D1A","KIR2DL4","IKZF3","FAM43A","MAPK1","GFOD1","DZIP3","NFKBIZ","PDE4B","AREG","HSPH1","ERGIC1","CD247","ADAMTS17","STRBP","ADGRE5","FYN","CDYL","CHORDC1","VPS37B","P2RY11","PDE7A","ACTN4","FOSL2","CCL3","PLCH2","GZMK","FAM177A1","LIF","SMC4",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_NK_CD56+",
                  save="Lymphoid_NK_CD56+.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_NK_CD56+.markergene.1.xxx.pdf")
