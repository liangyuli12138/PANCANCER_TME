import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_CD8_Tm");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_CD8.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["MT2A","MT1X","LMNA","EZR","GPR183","IL7R","MT1E","MCL1","HSP90AB1","HSPA8","SLC30A1","ANXA1","VIM","HSPE1","YPEL5","MT1F","BTG1","MYADM","HSPD1","HSP90AA1","RUNX3","HECA","PRNP","FTH1","CD44","TPT1","MT2P1","S100A10","HSPH1","ELOVL5","HSPA1A","CDKN1A","PTMA","SLC7A5","TUBA1B","IDS","SENP3-EIF4A1","PABPC1","GABARAPL1","SKI","SRGN","DNAJB6","H2AFZ","TAGLN2","RGCC","YWHAZ","PNP","UBE2D3","CREM","PDE4B",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD8_Tm",
                  save="Lymphoid_CD8_Tm.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD8_Tm.markergene.0.xxx.pdf")

marker_genes_dict=["MKNK2","SQSTM1","GPBP1","MAP1LC3B","FOSL2","FKBP4","CACYBP","CRIP1","NDRG1","MT1G","TNFAIP3","ELL2","CHD2","AQP3","ARL4C","TOP1","ATP1A1","UBB","EIF1","TUBA1C","PTGER4","GCLM","MAT2A","KDM6B","KLF6","SERP1","SLC2A3","SOD1","TSPYL2","JOSD1","CNBP","CSRNP1","AUTS2","TUBB2A","SPOCK2","PTGES3","UBC","FLNA","TUBA4A","AC008440.2","HNRNPAB","TSC22D3","JUND","RORA","MT1XP1","HNRNPA2B1","AHNAK","ZNF165","JMJD1C","HSP90AA2P",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD8_Tm",
                  save="Lymphoid_CD8_Tm.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD8_Tm.markergene.1.xxx.pdf")
