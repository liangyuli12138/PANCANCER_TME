import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Fibroblast_mCAF_KRT19");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Fibroblast.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["MYLK","TCEAL4","PEG3","ESR1","TCEAL2","ST13","KRT19","FHL2","MDK","CSRP1","SPINT2","NCAM1","PRLR","BEX2","RPS4X","RPL5","MYH11","TCEAL3","RBP1","OBSL1","PGR","RPL3","RPS3","HNRNPA1","MARCKSL1","RPL18","NRP2","TAGLN","GOLGA8B","DIRAS3","RPS8","TMOD1","AC006115.2","RHOB","LDLR","EMX2","MYLK-AS1","MT-CO1","TCEAL1","AR","SULF2","STAR","CSDC2","COL27A1","RAMP1","CNN1","LMOD1","FOXL2","CDH3","NPM1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Fibroblast_mCAF_KRT19",
                  save="Fibroblast_mCAF_KRT19.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Fibroblast_mCAF_KRT19.markergene.0.xxx.pdf")

marker_genes_dict=["ATP1B1","RPL15","TPM2","NAP1L1","GAS5","CAV1","MAOB","KLHDC8A","RPS3A","MASP1","CDC42EP3","ACTA2","SPEG","PALLD","TSPAN5","SYTL4","TGM2","GJA1","EIF4B","LONRF2","SYNPO2","EDNRA","SNCAIP","PPP1R14A","EMX2OS","SIGLEC11","KCNMA1","RPS6","EPDR1","SLC7A2","ALDH7A1","SCD5","DES","LRIG1","PTPRF","PMEPA1","ACSM3","PAWR","BTF3","LDHB","NACA","KCNK6","WT1","RPS7","PPM1K","FADS2","LEFTY2","RPL9","NEFM","GATA4",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Fibroblast_mCAF_KRT19",
                  save="Fibroblast_mCAF_KRT19.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Fibroblast_mCAF_KRT19.markergene.1.xxx.pdf")
