import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_CD4_Tn");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_CD4.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["GIMAP7","SNORD3A","GIMAP4","GIMAP1-GIMAP5","RN7SK","IKZF3","AC004585.1","DENND2D","TRAF3IP3","GIMAP6","TCF7","SORCS2","SMCO4","GIMAP1","IL6R","RBL2","TC2N","FKBP5","SLFN5","RN7SL1","GNG4","SFXN1","TOX2","TRG-AS1","AC109326.1","LINC00861","RN7SL4P","GIMAP8","CD3G","NUCB2","TTC14","DRAIC","RNU4-2","CXCR5","GRAP2","TRGC2","TRIM56","CHI3L2","ST8SIA1","IL16","ABLIM1","SCARNA2","GIMAP2","LINC00892","SDAD1","ADD3","TNFAIP8L2","CCR2","FZD3","DDX17",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD4_Tn",
                  save="Lymphoid_CD4_Tn.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD4_Tn.markergene.0.xxx.pdf")

marker_genes_dict=["MBP","C2orf68","S1PR1","SPN","MHENCR","CD84","SNHG22","DIDO1","PEBP1","FP671120.4","ZNF22","GIMAP5","SGCE","TSPOAP1-AS1","NUP43","RIPOR2","CAV1","NSUN5P2","ST6GAL1","AL355816.2","ARMH1","LRMP","ASCL2","SAMD3","ANKRD34B","PAQR8","PVALB","RNU1-67P","GZMA","STAG3","PTPN11","ZNF512","RAB27A","PRMT2","SESN1","ELK4","XXYLT1","GZMK","F2R","TXK","TMEM223","AC253572.2","GVINP1","AC243829.1","NKG7","GFOD1","AC135983.2","TMEM204","AC022239.1","GDF7",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD4_Tn",
                  save="Lymphoid_CD4_Tn.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD4_Tn.markergene.1.xxx.pdf")
