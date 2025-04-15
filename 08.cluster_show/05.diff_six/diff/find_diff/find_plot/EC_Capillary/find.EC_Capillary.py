import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/EC_Capillary");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.EC.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["RGCC","CA4","PLVAP","EXOC3L2","ITIH5","TMEM88","GSN","EHD4","CEACAM1","FLT1","SLCO2A1","TMEM204","VWA1","IRX3","RAMP3","HSPG2","DNASE1L3","BEX5","LNX1","PLAU","PDGFD","RFLNB","BTNL9","DHRS3","PODXL","KCNIP4","ACE","RBP7","RARB","ITGA1","OSBPL1A","INMT","TMEM233","ATOH8","RBP5","SLC51B","ARHGAP18","NOTCH4","TM4SF18","SLC9A3R2","FABP5","JCAD","FAM167B","TMEM150C","ETS1","TCEAL2","GPRC5B","NCKIPSD","EDNRB","CD34",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Capillary",
                  save="EC_Capillary.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Capillary.markergene.0.xxx.pdf")

marker_genes_dict=["KDR","SEC14L1","CD320","GRB10","LIMCH1","AC119674.2","TIMP3","LIPE-AS1","PTPRG","FXYD6","CD36","DDIT4","MCF2L","NQO1","QRFPR","SORBS2","CFAP20","TIAM2","PDGFB","C11orf96","TTC9","PLPP3","EPAS1","AC104211.2","SPTBN1","CAVIN2","PREX1","TRIL","ARHGAP29","ISLR","NPY1R","OLFML2A","AC020659.1","RAPGEF2","CYP4B1","FHL1","TMEM94","SYN3","TACC1","EFNB1","RXRA","GAS6","EML1","CETP","CD300LG","UACA","PITPNM2","BTN3A2","FABP5P7","JUP",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Capillary",
                  save="EC_Capillary.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Capillary.markergene.1.xxx.pdf")
