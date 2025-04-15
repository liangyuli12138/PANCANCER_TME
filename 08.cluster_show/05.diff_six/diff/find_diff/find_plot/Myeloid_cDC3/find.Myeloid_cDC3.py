import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Myeloid_cDC3");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Myeloid.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["LAMP3","BASP1","EEF1A1","CD83","DAPP1","RFTN1","CCR7","BIRC3","CFLAR","TRAF1","LSP1","BTG1","FSCN1","FNBP1","PTPN1","GPR157","KDM2B","TMSB10","DUSP5","TRABD2A","NUB1","ZFAS1","VOPP1","REL","SYNGR2","CLIC2","RNASEK-C17orf49","CDKN1A","RHOF","CERS6","RASSF4","STK4","NAP1L1","MARCKS","FAM49A","GSTP1","RPS19","PFN1","RELB","TNFAIP2","RPS27A","TSPAN33","EIF1","ACTG1","MARCKSL1","RPL35A","CD40","PNRC1","RPL10","MALAT1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_cDC3",
                  save="Myeloid_cDC3.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_cDC3.markergene.0.xxx.pdf")

marker_genes_dict=["HINT1","CYTH1","UVRAG","TRIP10","RPS2","ETV3","RUFY3","SLCO5A1","TRAFD1","USP12","ERICH1","BCL2A1","ID2","GRSF1","GPX4","LAD1","RAB8B","N4BP2L1","PABPC1","MREG","IRF1","IRF4","TBC1D4","IDO1","RPS7","RPL41","SNN","RPS27","JAK1","SINHCAF","CD80","MGLL","RPL9","RAB9A","SUB1","ANKRD33B","HLA-F","IL2RG","CD274","KIF2A","TNFAIP8","GSN","RPS29","BIRC2","TTYH2","SPECC1L-ADORA2A","CEP350","MCOLN2","LRRFIP2","GPBP1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_cDC3",
                  save="Myeloid_cDC3.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_cDC3.markergene.1.xxx.pdf")
