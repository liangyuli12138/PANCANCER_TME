import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_CD4_Tcm");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_CD4.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["YPEL5","PABPC1","EZR","GPR183","CCR7","CD55","RPS12","SARAF","RPL9","RPL32","TPT1","CHD2","RPS3A","RPL21","RPS6","USP36","RPL30","RPL34","SLC2A3","RPL3","RPS20","RPLP2","RPS27","RPS13","KDM6B","FAM177A1","RPL11","RPL13","ATP1B3","RPS14","RPS3","FTH1","PTMA","RPS18","RPL38","SCML1","RPL13A","ELL2","RPL31","EIF1","RPS27A","HECA","RPS16","BTG1","ZBTB10","ARIH1","CREM","RPS2","MARCKSL1","RPS9",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD4_Tcm",
                  save="Lymphoid_CD4_Tcm.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD4_Tcm.markergene.0.xxx.pdf")

marker_genes_dict=["RPS23","FAU","PER1","RPL35A","RPL10","RPS4X","RPL18A","RPL4","SC5D","RPS21","RPS28","EEF1A1","RPL27","EEF1G","RPL19","VPS37B","ARL4C","RPL36","RPS15A","RPS29","RPL18","RPL37","GAS5","RPS25","RPS8","ZNF331","HNRNPU","RPL5","PNRC1","DDX21","RPL23A","RANBP2","MKNK2","IL7R","AC013394.1","RPL27A","NR4A2","EIF4A3","NFKB1","SELENOK","RPL7","SYAP1","EEF1B2","SFPQ","RPS5","GNL3","NAP1L1","ARID5A","GABARAPL1","RPS17",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD4_Tcm",
                  save="Lymphoid_CD4_Tcm.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD4_Tcm.markergene.1.xxx.pdf")
