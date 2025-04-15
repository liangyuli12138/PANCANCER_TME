import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_CD8_Tisg");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_CD8.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["IFIT3","OASL","IFIT2","LIPA","HERC5","ISG15","IFIT1","ISG20","IFI44","PMAIP1","FGFBP2","DDX58","IFIH1","ZC3HAV1","OAS2","S1PR5","APOL6","FCGR3A","GBP4","HERC6","EPSTI1","MX1","GNLY","RNF149","IFI44L","MX2","CMPK2","OAS1","IFITM2","PARP14","USP18","TRIM26","RSAD2","XAF1","SH2B3","IFI6","TYROBP","PRSS23","AIM2","LPIN2","AC004551.1","IFIT5","OAS3","ZNFX1","RTP4","LY6E","DDX60","IFI16","PPM1K","GBP5",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD8_Tisg",
                  save="Lymphoid_CD8_Tisg.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD8_Tisg.markergene.0.xxx.pdf")

marker_genes_dict=["TGFBR3","RHEBL1","PARP9","PLAAT4","STAT1","EIF2AK2","TRANK1","RAB29","CD300A","SAMD9","GBP1","KLRF1","IFITM1","PLAC8","IGFBP3","ZBP1","SPON2","LINC02384","ITGAM","HELB","GNA13","CX3CR1","SAMD9L","SP140","GZMH","MIA3","B4GALT5","N4BP1","TNFSF10","ETV3","PLSCR1","LGALS9","DHX58","SP110","IFITM3","NEDD9","PLEK","HSH2D","C1orf21","S100A10","DRAP1","IFI35","DTX3L","APOL1","DDX60L","AC118553.2","LYSMD2","SAMHD1","EFHD2","TKTL1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD8_Tisg",
                  save="Lymphoid_CD8_Tisg.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD8_Tisg.markergene.1.xxx.pdf")
