import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Epithelium_Malig_Basal");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Epithelium.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["MYL6","ANXA2","PKM","YWHAZ","ENO1","RAC1","TAGLN2","AP000350.4","AC093512.2","TPI1","RPS20","MYL12B","KRT19","LDHA","C4orf3","PPIA","S100A6","CLIC1","ACTG1","H3F3A","ATP5F1E","CFL1","TMSB10","PGK1","RPS11","SDC1","NACA","CDC42","ANXA11","ARPC3","KRT18","CAST","S100A11","CHCHD2","RPL4","SELENOW","FTH1","GPX1","RAB11A","EIF1","PPA1","TMBIM1","S100A10","SUMO2","PGAM1","RPL8","RPL30","RPLP0","ROMO1","CD24",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Malig_Basal",
                  save="Epithelium_Malig_Basal.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Malig_Basal.markergene.0.xxx.pdf")

marker_genes_dict=["ELOB","VAMP8","LGALS3","SERF2","PABPC1","MYL12A","OST4","CTNNA1","EIF6","RAN","HSBP1","TMSB4X","UBA52","CLDN7","NOP10","ACTB","RPS13","JUP","ARPC2","RPL18","YWHAB","TXNDC17","SENP3-EIF4A1","PTTG1IP","HMGA1","TPT1","RPL28","VDAC1","PSMA7","SEM1","GAPDH","GUK1","SPINT2","SERBP1","B2M","SPINT1","RPSA","FAU","RPS21","CD63","EIF3E","RPL15","AC106886.5","RPS19","BCAP31","POMP","BTF3","RHOA","RPL37","LAMTOR5",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Malig_Basal",
                  save="Epithelium_Malig_Basal.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Malig_Basal.markergene.1.xxx.pdf")
