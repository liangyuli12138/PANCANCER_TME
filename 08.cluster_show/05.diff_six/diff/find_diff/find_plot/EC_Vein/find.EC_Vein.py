import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/EC_Vein");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.EC.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["ACKR1","CLU","CADM3-AS1","CCL14","SELP","GAS5","RPS2","DUSP23","RPS4X","RPS18","RPS6","RPL34","SELE","OLFM1","RPS23","RPS3","RPL32","RPL3","ACTN1","RPL13","RPL10A","RPL12","ZNF385D","RPS3A","RPS14","RPL5","NFKBIA","RPLP0","ADIRF","RPL31","RPL30","RPL18","EEF1G","RPL10","RPL37","EEF1B2","RPS9","CST3","RPL7","RPL36","TXN","RPL9","MAP3K8","UGCG","TPT1","RPL35A","RPS12","NAMPT","RPS27A","RPL17-C18orf32",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Vein",
                  save="EC_Vein.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Vein.markergene.0.xxx.pdf")

marker_genes_dict=["NPM1","RPL7A","RPL8","RPS15A","RPS19","AC116533.1","NNMT","NACA","RPL24","SOD2","EEF1A1","RPS13","CDKN1A","CTSC","YBX3","RPS17","RACK1","RPS21","RPL11","RPS8","RPS5","RPS24","RPL4","IL1R1","RPS7","NDRG1","RPL41","PABPC1","MEOX1","RPL26","RPL13A","RPL36A","NOP53","CPXM2","AC244100.2","SLC25A25","TPD52L1","ICAM1","TMEM273","C1orf112","RPL19","RPS16","FBL","MAT2A","RPL15","SLC2A3","RPL18A","SNHG8","EIF3E","RPS25",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Vein",
                  save="EC_Vein.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Vein.markergene.1.xxx.pdf")
