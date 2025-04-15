import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_CD8_Tstr");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_CD8.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["RPL13P5","LRRC23","FOSB","ERCC1","JUN","EGR1","BRD2","IER2","JUNB","HEXIM1","DNAJA1","NR4A1","DNAJB1","UBE2S","SERTAD1","BTG2","HSP90AA1","H3F3B","GADD45B","DSTNP2","UBB","PPP1R10","HSPA1B","TUBB4B","DUSP1","FOS","EHD1","MALAT1","ZFP36","HSP90AB1","DDX5","SNHG5","NXF1","ATF4","SRSF7","DUSP2","HSPA8","SRSF3","DDIT3","SENP3-EIF4A1","CITED2","H2AFX","SRSF5","ZFAS1","TPI1","EIF1","HNRNPH1","TRA2B","HSPA1A","ID2",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD8_Tstr",
                  save="Lymphoid_CD8_Tstr.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD8_Tstr.markergene.0.xxx.pdf")

marker_genes_dict=["PPP1R15A","EIF4A2","SNHG12","NR4A2","EIF4A3","HSP90AA2P","HERPUD1","INTS6","HOOK2","HSPH1","JUND","SAT1","SRSF2","IER5","MYLIP","CD83","FUS","MAT2A","DDX3X","SAP18","EIF4A1","TUBA1A","SLC20A1","MORF4L2","MKNK2","MARCKSL1","SNHG1","AC022217.3","MYADM","RBBP6","TFRC","AMD1","EIF5","LBH","AC020916.1","IFRD1","CD69","RHOB","RSRP1","ZFP36L1","NASP","PCF11","ZNF331","POLR2A","ATP1A1","NUFIP2","MIDN","SFPQ","AC016876.2","TUBA1C",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD8_Tstr",
                  save="Lymphoid_CD8_Tstr.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD8_Tstr.markergene.1.xxx.pdf")
