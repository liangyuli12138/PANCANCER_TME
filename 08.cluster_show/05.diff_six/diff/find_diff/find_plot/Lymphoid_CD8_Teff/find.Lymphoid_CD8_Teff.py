import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_CD8_Teff");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_CD8.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["JUN","DUSP1","FOS","IER2","HSPA1B","FOSB","GADD45B","DNAJB1","EGR1","ERCC1","JUNB","BTG2","AC022217.3","AOAH","CCL5","PEBP1P3","NR4A1","DUSP2","FCRL6","BRD2","GZMK","CD69","RHOB","CMC1","CITED2","HSPA1A","ID2","TRGC2","TNF","TRDC","ITGAL","HOOK2","PLEK","SRSF7","SRSF3","IFNG","AL499604.1","ANKRD37","KLRG1","PPP1R15A","TRA2B","CD160","ZFP36","HIST1H2BG","HSPE1-MOB4","DNAJA1","IER3","AL592429.2","TRG-AS1","ARRDC3",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD8_Teff",
                  save="Lymphoid_CD8_Teff.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD8_Teff.markergene.0.xxx.pdf")

marker_genes_dict=["PPP1R10","AL691403.1","ZFAND2A","SNAI1","EOMES","SLAMF7","HCST","RASD1","MYL12A","MATK","MYO1F","AC015849.1","NFKBIA","HSPA1L","SLAMF6","HIST1H3J","EIF4A2","ARHGAP25","ILF3-DT","SOCS1","DNAJB4","CD244","RBKS","PYHIN1","GZMA","PTPRCAP","TXNIP","DDIT3","GZMH","GIMAP4","HSPA6","RESF1","AC007384.1","AP001160.1","IL32","F2R","NR4A2","DUSP6","AL031777.3","PDE4A","IKZF2","ICAM1","DTHD1","EGR2","BTN3A1","DDIT4","LINC01841","KLF2","TBC1D10C","RBL2",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD8_Teff",
                  save="Lymphoid_CD8_Teff.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD8_Teff.markergene.1.xxx.pdf")
