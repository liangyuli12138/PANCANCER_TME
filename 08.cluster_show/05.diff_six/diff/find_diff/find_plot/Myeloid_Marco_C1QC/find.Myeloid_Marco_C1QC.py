import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Myeloid_Marco_C1QC");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Myeloid.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["C1QB","C1QA","AC016876.2","FCGR3A","C1QC","CD81","FTL","TYROBP","APOC1","SERPING1","CYBA","UBE2L6","B2M","GPX1","DNASE2","MARCO","CTSS","NPC2","LIPA","VAMP8","LGALS3BP","CTSD","CTSB","FCER1G","MGST3","GRN","SLC7A7","CD68","AC007192.1","STAT1","S100A11","LY6E","APOE","PSME2","GIMAP4","FCGR1A","SERPINA1","COMT","ALDH2","PLBD1","IGSF6","CYP27A1","CREG1","CTSA","GSTO1","FCGRT","PPT1","HLA-C","BCAP31","SAMHD1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_Marco_C1QC",
                  save="Myeloid_Marco_C1QC.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_Marco_C1QC.markergene.0.xxx.pdf")

marker_genes_dict=["BRI3","TMSB4X","ARPC1B","C2","PSAP","SCPEP1","ACP2","NOP10","AC004922.1","ASAH1","VSIG4","SCARB2","ATP6V1F","Z84488.2","GPNMB","AIF1","CTSC","CTSH","IFI27","MDH1","TRIM14","TNFSF12-TNFSF13","IFIT3","UQCR10","ELOB","GBP1","LAP3","LRPAP1","SQOR","SMIM30","SCD","VAMP5","OAS1","GAA","ATP5MC3","CALHM6","HLA-DMB","SLC31A2","GCHFR","LGALS3","LAMP2","UQCR11","NR1H3","CNDP2","SAMD9L","GBP4","IFI6","HEXB","APOL1","ACP5",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_Marco_C1QC",
                  save="Myeloid_Marco_C1QC.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_Marco_C1QC.markergene.1.xxx.pdf")
