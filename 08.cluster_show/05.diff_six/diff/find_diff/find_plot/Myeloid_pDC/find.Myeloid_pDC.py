import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Myeloid_pDC");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Myeloid.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["BX571818.1","TRAF4","CXCR4","EZR","LILRA4","C12orf75","PPP1R16B","TCF4","RPS27","SPIB","PLD4","SEMA7A","GZMB","IRF4","BCL11A","PLAC8","IRF7","RPL23A","YPEL5","RAB11FIP1","PIM2","SEL1L3","CCDC50","RPS23","CYB561A3","RPS3A","MZB1","PTPRCAP","SEPTIN6","SEC61B","PTPRS","SLC7A5","ITM2C","NOP58","RPL3","RPS11","SLC38A1","SULF2","IRF8","SLC15A4","CYFIP2","RPS12","RPL13A","CARD11","RBM38","RPLP0","AC005912.1","NIBAN3","SELL","RPS6",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_pDC",
                  save="Myeloid_pDC.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_pDC.markergene.0.xxx.pdf")

marker_genes_dict=["RPS5","RPL18","HERPUD1","CHAF1A","CORO1A","RPS18","RPL10A","GPR183","IRF2BP2","RPS4X","TSPAN13","RPL13","CCDC69","LPIN1","RPS15A","RPS3","MKNK2","RPL6","SLAMF7","CYTIP","AC244205.1","RUBCN","RPS8","MALT1","FYTTD1","CD79A","DUSP5","NAPSB","DERL3","RPL18A","RPSA","AP001324.1","AC245884.11","LDLRAD4","POLB","RPL21","CLEC4C","HINT1","OGT","RPL37","RPS10-NUDT3","RIPOR2","AC024293.1","RPL19","FAU","TXNDC5","RPS7","CLIC3","RPS23P8","C12orf45",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_pDC",
                  save="Myeloid_pDC.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_pDC.markergene.1.xxx.pdf")
