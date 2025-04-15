import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/EC_Artery");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.EC.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["SEMA3G","IGF2","AIF1L","SRP14","CXCL12","ASS1","SRGN","ARL15","KCTD12","GLUL","MECOM","PPA1","FBLN2","ICAM2","CLDN5","VEGFC","JAG2","GJA4","CAV1","PALMD","FAM107A","CLIC5","PCSK5","FBLN5","UTRN","TSPAN2","ANXA3","TIMP3","S100A6","JAG1","PLLP","SERPINE2","PPP1R14A","SSUH2","LTBP4","CRIP2","SLC6A6","PODXL","HEY1","SLC9A3R2","FN1","MAST4","PI16","DEPP1","MRPL33","SOX17","TAGLN2","NEBL","PPT2-EGFL8","LPCAT2",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Artery",
                  save="EC_Artery.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Artery.markergene.0.xxx.pdf")

marker_genes_dict=["MAP4","PLS3","CCDC3","SYN3","TMEM178A","DST","ADAM15","SYNPO","CLEC14A","F11R","A2M","CGNL1","ADAMTS6","CRIP1","EMP3","GUCY1A1","OCIAD2","ABI3BP","KCNN3","CAVIN1","AQP7","SYT11","ITPR2","FUT8","EDN1","PRDM16","CFH","MTUS1","PLCG2","ATP13A3","NES","TM6SF1","IL32","EFNB2","DHH","THSD7A","PKN3","STMN1","ERC1","ACTN4","PIK3R3","RAMP2","RHOC","STOM","IGFBP3","LGALS3","PDCD4","SYS1-DBNDD2","MMRN2","CRYBG3",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Artery",
                  save="EC_Artery.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Artery.markergene.1.xxx.pdf")
