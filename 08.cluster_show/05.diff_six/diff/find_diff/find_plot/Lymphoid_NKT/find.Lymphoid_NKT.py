import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_NKT");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_NK.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["GZMH","FGFBP2","CD3D","S100A10","CD52","CD8B","CD8A","CD3G","NKG7","PATL2","ZEB2","ADGRG1","SYNE1","HLA-DPB1","GZMB","ARL4C","LITAF","FLNA","ITGB1","TRGC2","S1PR5","TGFBR3","LILRB1","EMP3","PXN","C12orf75","MIAT","SYNE2","HLA-DRB1","AC243829.2","S100A4","GNLY","ENC1","S100A6","RUNX3","LGALS1","PPP2R5C","CDC25B","TNFRSF1B","SH3BGRL3","MT1E","CX3CR1","CLEC2D","S1PR1","HLA-DRB5","C1orf21","LINC02446","MT2A","KLF3","CAPN2",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_NKT",
                  save="Lymphoid_NKT.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_NKT.markergene.0.xxx.pdf")

marker_genes_dict=["VCL","MT1X","ANXA1","OASL","PRSS23","AL136018.1","RIPOR2","GNG2","SSBP3","TRBC2","FGL2","DSTN","RAP1GAP2","TMBIM1","RASA3","CST7","AHNAK","SELPLG","STK38","RASGRP2","HLA-DQA1","PLEK","CRIP1","LYAR","LAIR2","YBX3","HLA-DPA1","RAB29","TTC16","HLA-DQA2","EFHD2","LRRFIP1","AGPAT4","MT1F","PROK2","FCGR3A","B3GAT1","SYT11","ORAI1","KLRG1","HLA-DRB6","S1PR4","ANXA2","RAP2A","LINC02384","MYO1G","PLEKHG3","SH3BP5","ARRB2","CD5",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_NKT",
                  save="Lymphoid_NKT.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_NKT.markergene.1.xxx.pdf")
