import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_B_memory");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_B.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["HLA-DRA","CD83","HLA-DRB1","REL","HLA-DQA1","BANK1","HLA-DPB1","HLA-DQB1","MS4A1","TUBA1A","CD74","LAPTM5","CD37","CD82","HLA-DPA1","LTB","HLA-DRB5","CCR7","SWAP70","PIKFYVE","BCL2A1","HLA-DRB6","IRF8","SLC2A3","RILPL2","NAP1L1","RPS27","CMTM6","CXCR4","NAPSB","CD69","ZFP36L1","SPIB","NFKBID","CD52","GRASP","WSB1","CXCR5","COTL1","CD40","SEMA7A","NR4A2","RUNX3","LCP1","KDM6B","PTMA","GPR183","DEK","DDX21","SNN",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_B_memory",
                  save="Lymphoid_B_memory.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_B_memory.markergene.0.xxx.pdf")

marker_genes_dict=["TRIM22","RIPOR2","STAT6","RPL30","RPLP2","PEA15","YWHAZ","TMC8","HLA-DQA2","IDS","ZFAS1","RPS20","RPS29","LY9","CHORDC1","PNRC2","FCMR","RASGRP2","LRRFIP1","RHOG","SFPQ","VASP","ARL6IP1","RPS21","RPL31","CDKN2D","TAGLN2","MAP4K4","PARP15","RPS23","RUBCNL","WEE1","TRAF4","HLA-DMB","SH2B3","CD53","BLK","TOB2","SNED1","INPP5D","TMEM273","HSPA4","TTN","FAM102A","RPS12","SMIM14","RPL9","RPL38","TNFRSF1B","PMAIP1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_B_memory",
                  save="Lymphoid_B_memory.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_B_memory.markergene.1.xxx.pdf")
