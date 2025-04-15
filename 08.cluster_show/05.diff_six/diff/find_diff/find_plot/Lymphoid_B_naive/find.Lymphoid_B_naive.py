import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_B_naive");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_B.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["LAPTM5","MS4A1","HLA-DPB1","HLA-DRA","HLA-DRB1","CD37","IL4R","CD74","BACH2","HLA-DPA1","HLA-DMB","TCL1A","CD22","HLA-DQA1","HLA-DQB1","FCMR","NIBAN3","CD72","CD52","FCRL1","FCER2","HLA-DRB5","LINC00926","CD79B","CXCR4","RUBCNL","HLA-DOA","BCL11A","LIMD2","HVCN1","HLA-DRB6","FOXP1","AFF3","PAX5","TMEM123","CD83","SELL","PLEKHA2","RASGRP2","SNX2","PTPRC","BANK1","BTG1","RBM38","HLA-DMA","LTB","CIITA","RIPOR2","ETS1","ARRDC2",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_B_naive",
                  save="Lymphoid_B_naive.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_B_naive.markergene.0.xxx.pdf")

marker_genes_dict=["SWAP70","COL19A1","ZNF318","IRF8","CLEC2D","YBX3","CORO1A","GGA2","CYB561A3","EVL","GPR18","TNFRSF13C","LCP1","HLA-DQA2","SP100","SLC2A3","IRS2","SORL1","CD19","TRIM38","REL","STAT6","STX7","TMC8","CD69","SATB1","ZFP36L1","TRBC2","UCP2","AC025279.1","BCL7A","LGALS9","MAP4K4","TRIM22","STK17A","SMAD3","SLC44A2","RFTN1","PLAC8","CD53","PIK3IP1","VPREB3","VAV3","S1PR1","TRAF5","P2RX5","TMEM131L","CDCA7L","PIK3AP1","PLEKHG1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_B_naive",
                  save="Lymphoid_B_naive.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_B_naive.markergene.1.xxx.pdf")
