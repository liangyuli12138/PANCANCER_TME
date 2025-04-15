import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_CD4_Th17");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_CD4.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["Z94721.2","KLRB1","RPLP1","IL4I1","TPT1","CCL20","GPR183","VIM","TNFAIP3","GPR65","RORA","ANXA1","ANKRD28","RPLP0","IL7R","MYADM","MAP3K4","CLEC2B","ADGRE5","ABCB1","MGAT4A","ERN1","RPS19","PNP","AHR","MCL1","PDE4D","LINC01871","AQP3","EEF1A1","JAML","ARL4C","RPS12","FURIN","KLF6","PTGER4","PDE4B","FKBP11","S100A6","RUNX3","FLNA","NR1D2","RORA-AS1","SPTY2D1","ANXA2","MYBL1","HSPA8","TAGLN2","EEF1G","AMPD2",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD4_Th17",
                  save="Lymphoid_CD4_Th17.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD4_Th17.markergene.0.xxx.pdf")

marker_genes_dict=["RPL10","S100A10","PRNP","RPL13","CERK","RPS18","ELOVL5","GLIPR1","SNHG16","HECA","CRIP1","AC073610.2","EFHD2","AC116533.1","RPS27A","RNF149","CREM","FTH1","KIF5C","SERP1","TMEM273","RPL36A","AC008440.2","SRGN","RPS16","RPL31","RPL3","DPP4","RBMS1","RPS2","APOL3","AUTS2","RPS3","RPSA","RAP1B","LMNA","SLAMF1","RORC","BTG1","RPL34","CTSH","IL23R","RPLP2","CD69","RPS6KA3","TIFA","CDC42EP3","CAPG","RPL17-C18orf32","RPL9",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD4_Th17",
                  save="Lymphoid_CD4_Th17.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD4_Th17.markergene.1.xxx.pdf")
