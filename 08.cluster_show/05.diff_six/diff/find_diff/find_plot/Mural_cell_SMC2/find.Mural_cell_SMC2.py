import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Mural_cell_SMC2");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Mural_cell.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["RAMP1","PALLD","MYLK","TCEAL4","NCAM1","TIMP2","EEF1A1","COLEC12","PSAP","CNN1","TCEAL3","LTBP2","RPS24","CTSB","EPDR1","S100A10","FHL2","PSD","CLU","PBX1","DPYSL3","MDK","SFRP4","CES1","TGM2","HSD17B6","KCNMA1","GREM1","MFAP5","TUBA1A","SEMA3C","MYLK-AS1","CA12","PGR","CD81","PTGDS","CDC42EP3","PGM5","CYP1B1","FBLN1","ESR1","BEX3","CBR4","LAPTM4A","LRIG1","TGFBR3","MTCH1","OLFML3","RHOB","GJA1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Mural_cell_SMC2",
                  save="Mural_cell_SMC2.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Mural_cell_SMC2.markergene.0.xxx.pdf")

marker_genes_dict=["RPLP1","REXO2","ACTG2","COL12A1","COL27A1","RAI14","CSRP1","BAG2","FRMD6-AS2","TMED10","COL5A1","PRSS23","SFRP1","MORF4L2","CKB","NUCB2","ATP1B1","OGN","RPL5","SCARA3","RPS12","ANXA2","ITGA8","NRP2","OAZ1","DES","TCEAL1","MZT2B","LGALS1","TMEM123","MT-CYB","CLMP","NPTN","PPIB","EMP3","TLE5","PRNP","FN1","SEC61B","CETN2","NBL1","PCDH7","PAM","BNC2","PIM1","RPS4X","PDLIM7","TPM1","SRP14","NLRP1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Mural_cell_SMC2",
                  save="Mural_cell_SMC2.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Mural_cell_SMC2.markergene.1.xxx.pdf")
