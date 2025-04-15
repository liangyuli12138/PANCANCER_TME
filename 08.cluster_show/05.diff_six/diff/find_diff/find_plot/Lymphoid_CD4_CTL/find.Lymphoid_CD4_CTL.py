import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_CD4_CTL");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_CD4.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["CCL5","DNAJB1","FOSB","HSPA1B","JUN","HSPA1A","ERCC1","EGR1","HSP90AA1","ZFP36","BTG2","DUSP1","GZMA","NR4A1","CD69","PPP1R15A","DNAJA1","UBB","IER2","FOS","MYADM","DUSP2","HSPA8","HSP90AB1","GZMK","JUNB","H3F3B","SERTAD1","UBE2S","ID2","HSPH1","ANXA1","SRSF7","HSP90AA2P","CITED2","BRD2","AC015849.1","ZFP36L2","NFKBIA","CDKN1A","IER5","AC022217.3","TNFAIP3","GADD45B","NR4A2","NEU1","RGS2","IFNG","AC020916.1","TNF",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD4_CTL",
                  save="Lymphoid_CD4_CTL.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD4_CTL.markergene.0.xxx.pdf")

marker_genes_dict=["ZC3H12A","UBC","CD83","ATF3","MYLIP","HSPD1","HSPA5","EHD1","ANKRD37","AC008440.2","TRA2B","H2AFX","ZFAND2A","JUND","SRSF3","HSPA6","SLC20A1","SENP3-EIF4A1","TUBA1A","IFNG-AS1","MCL1","DUSP5","NKG7","MAFF","KLF2","PNP","AL592429.2","TUBB4B","PARP8","TUBA1C","MYBL1","CXCR4","SLAMF7","EIF4A3","ERN1","EIF4A2","TGFB1","SLC1A5","XBP1","ETV3","MIDN","EIF1","DDIT3","TUBB2A","SNHG12","ATF4","AC025259.3","NFKBIZ","JMJD6","SAMD3",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD4_CTL",
                  save="Lymphoid_CD4_CTL.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD4_CTL.markergene.1.xxx.pdf")
