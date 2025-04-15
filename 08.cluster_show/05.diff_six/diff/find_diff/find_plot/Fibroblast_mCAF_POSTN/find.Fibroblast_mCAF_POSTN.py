import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Fibroblast_mCAF_POSTN");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Fibroblast.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["BGN","SPARC","TAGLN","CTHRC1","COL1A1","COL3A1","COL8A1","COL10A1","VCAN","RAB31","INHBA","IGFBP7","MYL9","FN1","COL5A2","ACTN1","FAP","COL1A2","F2R","POSTN","ASPN","TPM2","ACTA2","COMP","LOXL1","TPM1","HOPX","LGALS1","COL5A1","PMEPA1","NREP","LTBP2","PALLD","SERPINH1","TPM4","SDC1","SULF1","MOXD1","ITGB1","MYL12A","LRRC15","STRA6","CCN2","MXRA5","CDH11","COL11A1","ANTXR1","MYL6","TNFRSF12A","MYH9",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Fibroblast_mCAF_POSTN",
                  save="Fibroblast_mCAF_POSTN.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Fibroblast_mCAF_POSTN.markergene.0.xxx.pdf")

marker_genes_dict=["MICAL2","DIO2","CRIP2","CCN4","ANO1","ATP5F1E","UNC5B","COL4A2","FNDC1","ARL4C","THY1","CFL1","PDLIM7","CNN2","MMP11","NOTCH3","C1QTNF6","TIMP1","ENC1","S100A16","LAMP5","LOXL2","MFAP2","ISG15","ITGA1","CYTOR","MSRB3","C4orf48","ARPC2","KIF26B","IL32","AC079328.2","GJB2","VCAN-AS1","COL8A2","RGS3","IFI27","HES4","ADAM12","EDNRA","COL4A1","PXDN","STARD4-AS1","NTM","PYCR1","TNC","SLC12A8","KDELR3","DKK3","SFRP4",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Fibroblast_mCAF_POSTN",
                  save="Fibroblast_mCAF_POSTN.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Fibroblast_mCAF_POSTN.markergene.1.xxx.pdf")
