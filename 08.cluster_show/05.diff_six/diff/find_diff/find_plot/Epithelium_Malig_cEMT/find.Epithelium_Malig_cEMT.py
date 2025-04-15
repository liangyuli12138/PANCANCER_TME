import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Epithelium_Malig_cEMT");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Epithelium.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["TPM2","CAVIN1","COL17A1","FLNA","CALD1","SPARCL1","MYL9","DST","ACTA2","MYLK","ACTN1","DKK3","CSRP1","C1R","SYNM","LTBP2","TAGLN","ITGB1","SERPING1","ANO1","PRNP","CNN1","AMIGO2","MSRB3","SCPEP1","GEM","NNMT","SPARC","NFIB","SEMA3C","SERPINF1","MYH11","ACTG2","NGFR","MYLK-AS1","VIM","EFEMP1","FAM126A","FBXO2","SFRP1","FSTL1","TPM1","CAV1","WLS","ECRG4","FBXO32","CDH3","STK17A","C1S","ARL4A",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Malig_cEMT",
                  save="Epithelium_Malig_cEMT.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Malig_cEMT.markergene.0.xxx.pdf")

marker_genes_dict=["IFI16","YBX3","KRT14","MMP2","BOC","CFH","GAS6","SSPN","PDLIM3","ARL4C","TNC","IFITM3","SPRY2","MME","APOE","LAMC1","STAMBPL1","MFGE8","ANXA1","BST2","AC083870.1","PIK3R1","HTRA1","LGALS1","DPYSL3","MXRA8","GLIPR1","LMOD1","LY6E","GPR87","TP63","FGFR1","CXCL14","AEBP1","ARID5B","TUBB6","ACTA2-AS1","IFITM2","TUBA1A","COL6A2","MIA","TSHZ2","NES","IRX1","SNAI2","EDNRA","SULF1","KCNMB1","STAC2","TIMP3",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Malig_cEMT",
                  save="Epithelium_Malig_cEMT.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Malig_cEMT.markergene.1.xxx.pdf")
