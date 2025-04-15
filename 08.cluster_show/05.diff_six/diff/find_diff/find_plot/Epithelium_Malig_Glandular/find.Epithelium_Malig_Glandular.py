import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Epithelium_Malig_Glandular");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Epithelium.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["TMED3","RABAC1","TMED9","TTC3","TFF3","ST6GALNAC1","SPINK4","LY6E","FCGBP","HPCAL1","AGR2","MUC2","CKAP4","WFDC2","IFITM3","SH3BGRL3","CTSD","SERF2","KCNMA1","ATP5F1E","SELENOM","PPIB","AP000350.4","SSR4","CD63","CIB1","IFITM2","SERPINA1","TCEA3","TMEM205","NDUFA13","GPX1","PRDX5","GALNT6","CYBA","HLA-A","SEC61G","RPL37","SMIM22","RPN2","TPI1","ATP6V0B","RPL36","TMED10","CHCHD2","PSAP","MUC5B","LDHA","BST2","GUK1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Malig_Glandular",
                  save="Epithelium_Malig_Glandular.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Malig_Glandular.markergene.0.xxx.pdf")

marker_genes_dict=["OSTC","KLK1","OCIAD2","SELENOH","REG4","OST4","METTL9","HSBP1","ROMO1","SNHG6","ATP6AP1","C3","PDIA3","ENO1","BCAP31","ARPC3","LGALS3BP","ITM2B","ITPR2","SPDEF","PKM","TMEM208","CNPY2","RCN1","CD24","NOP10","RAB2A","TMSB10","KRT18","IFI27","SLC40A1","EIF2S2","RPN1","ARPC1B","HM13","C4orf48","PGK1","PRDX4","UQCRQ","RPS11","TSTA3","RPS20","MYRIP","CD82","LMAN2","NPDC1","SEPTIN9","SELENOW","S100A13","GCHFR",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Malig_Glandular",
                  save="Epithelium_Malig_Glandular.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Malig_Glandular.markergene.1.xxx.pdf")
