import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/EC_Glomerular");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.EC.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["IGFBP5","EMCN","PLAT","EHD3","SOST","CRHBP","SYNE1","TGFBR2","PLPP1","PLPP3","ITGA8","AC119674.2","SLC14A1","ADGRF5","GJA5","MEIS2","TBX3","KDR","NTN4","SYNE2","PTPRB","AC007563.2","PBX1","HECW2","APP","HIPK2","ST6GAL1","CCND1","LINC00551","ANO6","MGP","LGALS3BP","TMEM204","SLC9A3R2","LIMCH1","MEG3","CEACAM1","PRX","GATA5","RXFP1","EFNB1","IL13RA2","MYO1C","SNED1","MYLK","RNASE1","RGL1","FGL2","BTNL9","SLC44A1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Glomerular",
                  save="EC_Glomerular.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Glomerular.markergene.0.xxx.pdf")

marker_genes_dict=["LDB2","IL6ST","EPAS1","CD74","NOSTRIN","HLA-DQA1","DNAJC1","RAPGEF4","GAS6","TEK","TP53I11","PIP4K2A","AJAP1","PDIA5","FKBP5","CHRM3","ANTXR2","SMAD6","FCN3","TMEM94","GIMAP8","RCSD1","SEC14L1","PPFIBP1","HLA-DQA2","EMID1","BEX5","SLC2A9","BST2","MTSS1","EXOSC7","TBC1D8","AC097459.1","F8","PSEN2","GRIP2","FXYD6","TMEM26","KHDRBS3","RCAN2","TNS3","RGS9","MEG8","S100A4","ZEB1","EML1","HLA-DRB6","FLT1","NCKAP5","ANKRD44",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Glomerular",
                  save="EC_Glomerular.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Glomerular.markergene.1.xxx.pdf")
