import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Mural_cell_SMC1");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Mural_cell.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["ADIRF","MUSTN1","RERGL","PLN","DSTN","SORBS2","AC006254.1","CRIP1","NET1","NDUFA4","CCDC3","CASQ2","RCAN2","SNCG","MYL9","TAGLN","CRIP2","CSRP2","NTRK2","SUSD5","MYH11","LBH","SOD3","MGST3","ADIRF-AS1","WTIP","NIBAN1","MFGE8","C12orf75","ARPC1A","SGCA","GGTA1P","KCNAB1","SBSPON","MCAM","PHLDA2","NRIP2","C11orf96","CX3CL1","S100A6","DLGAP4-AS1","MOB2","KCNA5","TBX2-AS1","GADD45B","HCFC1R1","CAVIN3","RBPMS2","NTRK3","LDB3",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Mural_cell_SMC1",
                  save="Mural_cell_SMC1.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Mural_cell_SMC1.markergene.0.xxx.pdf")

marker_genes_dict=["MAP3K7CL","NRGN","SYS1-DBNDD2","GPRC5C","ITIH3","SLC25A4","FHL5","ATP1A2","LDOC1","ITIH5","MT1M","AOC3","FRZB","HES4","MSRB3","ACKR3","MAPRE2","PGAM2","DACT3","PTP4A3","HACD1","AKAP1","ANAPC16","ENTPD3","CCN3","AC091182.1","TSC22D1","CD151","PLAC9","KLF2","SELENOW","CCDC107","AC113167.1","S100A4","NTN4","COX7A1","PRR34-AS1","MYOM1","PRPH","AIF1L","SNTA1","SPARCL1","PLS3","LSMEM2","DBI","CAP2","P2RX1","FOXC1","MT1E","ID1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Mural_cell_SMC1",
                  save="Mural_cell_SMC1.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Mural_cell_SMC1.markergene.1.xxx.pdf")
