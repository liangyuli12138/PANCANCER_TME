import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Fibroblast_apCAF_CD74");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Fibroblast.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["CD74","HLA-DRB1","HLA-DRB5","HLA-DRA","HLA-DRB6","ARHGDIB","IGFBP7","SRGN","LAPTM5","HLA-DPB1","HLA-B","LCP1","B2M","HLA-DPA1","HLA-A","TNFSF13B","PTPRC","HLA-F","GGT5","FMO2","NRP1","CD53","CXCR4","HLA-C","RAC2","CD48","IFITM2","HLA-DQB1","HLA-DQA1","TRBC2","RPS2","HLA-DMA","CD52","FYB1","CORO1A","AC116533.1","RPS18","RPL36A","RPS19","RUNX3","HCST","TNFRSF1B","ARHGAP15","RPS8","RPS29","ITGB2","RPSA","CD3D","PSME2","PLAAT4",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Fibroblast_apCAF_CD74",
                  save="Fibroblast_apCAF_CD74.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Fibroblast_apCAF_CD74.markergene.0.xxx.pdf")

marker_genes_dict=["BIRC3","HCLS1","TYROBP","LAP3","IL32","STK17B","RPL12","TAGAP","CD69","RPS12","DUSP2","CTSC","CYTIP","RPLP0","RPS21","BST2","RPL28","PTPRCAP","TMC8","RPS3","SAMSN1","CFI","NAP1L1","HLA-DMB","RBP5","PSMB9","CD2","RPS23","RPS27A","IL2RG","RPL3","STEAP4","PCED1B-AS1","RPS17","EZR","RPLP2","RPL8","EVI2B","NTRK2","RHOH","ACAP1","MARCKSL1","RPL41","SPOCK2","RPS4X","RPS6","RPL17-C18orf32","AC024293.1","UBD","RPL37",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Fibroblast_apCAF_CD74",
                  save="Fibroblast_apCAF_CD74.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Fibroblast_apCAF_CD74.markergene.1.xxx.pdf")
