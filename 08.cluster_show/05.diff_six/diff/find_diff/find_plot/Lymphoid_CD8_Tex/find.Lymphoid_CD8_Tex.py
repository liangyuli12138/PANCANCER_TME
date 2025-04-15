import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_CD8_Tex");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_CD8.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["DUSP4","CTLA4","TNFRSF9","RBPJ","RGS1","GAPDH","TIGIT","CD82","TNFRSF1B","AKAP5","LAG3","ITGAE","GEM","PARK7","ICOS","NAB1","NKG7","HAVCR2","ENTPD1","RAB27A","CXCR6","MYO7A","LINC01871","HLA-DRB1","IKZF3","AL390957.1","SLA","LYST","MIR155HG","PDCD1","CHN1","AHI1","CD27","PHLDA1","ITM2A","PGAM1","COTL1","SNX9","CXCL13","LRRN3","PKM","PRDM1","CD3D","SAMSN1","FABP5","LAYN","APOBEC3G","APOBEC3C","DUSP16","BATF",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD8_Tex",
                  save="Lymphoid_CD8_Tex.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD8_Tex.markergene.0.xxx.pdf")

marker_genes_dict=["KRT86","CD8A","PTMS","TNIP3","RNF19A","IGFLR1","TOX","CLEC2D","CTSD","CD2","NR3C1","SEM1","NFIL3","SH3BGRL3","MAST4","VCAM1","ZEB2","ARID5B","CD74","RGS2","CARD16","GOLIM4","NDFIP2","CCL3","HLA-DRB5","TBC1D4","HLA-DPA1","FAM3C","ATP8B4","GZMB","DGKH","HLA-DPB1","CYTOR","INPP5F","HLA-DRA","ADGRG1","ZBED2","HLA-DMA","LSP1","MS4A6A","CCND2","HNRNPLL","SYNGR2","CRTAM","TNFRSF18","TNFSF4","TRAF5","ZC3H12C","PAM","PTPN22",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD8_Tex",
                  save="Lymphoid_CD8_Tex.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD8_Tex.markergene.1.xxx.pdf")
