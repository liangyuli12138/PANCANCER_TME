import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/EC_Sinusoidal");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.EC.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["DNASE1L3","CTSD","FCN3","CLEC4G","CTSL","OIT3","FCN2","SELENOP","MS4A6A","CLEC1B","TIMP1","MRC1","ACP5","LGMN","STAB2","CD4","STAB1","FCGR2B","CRHBP","CLEC4M","CD36","CXCL16","ST6GAL1","MAF","BMPER","RBP1","AKAP12","GPR182","NPL","ANPEP","CCL14","Z84466.1","GSTO1","FEZ1","FCGRT","CD14","F2R","EXOSC7","CLEC3B","SLC40A1","PPFIBP1","C1QTNF1","GYPC","ITSN1","CETP","SHC2","TMEM37","PLAC8","VMO1","MYO10",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Sinusoidal",
                  save="EC_Sinusoidal.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Sinusoidal.markergene.0.xxx.pdf")

marker_genes_dict=["COTL1","PDE2A","LDB2","BST2","AP2S1","RELN","ART4","NPY1R","ABCA9","GATA4","FCHSD2","NTN4","PLTP","HECW2","CTSB","PVALB","IL33","CYBA","SORBS1","SNX6","SLC9A9","NUPR1","SNORD3A","SLC7A8","AL078590.2","NRG3","MIR99AHG","RAB20","F8","ADD3","SNX8","DENND4C","CCDC152","LINC01197","VAMP8","ST8SIA6","IFITM10","MAN1C1","HIBCH","HIPK2","GPM6A","ST6GALNAC3","AC007686.2","LGALS3BP","RIPOR3","MASP1","ITGA9","RGL3","AP2A2","PRICKLE2",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Sinusoidal",
                  save="EC_Sinusoidal.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Sinusoidal.markergene.1.xxx.pdf")
