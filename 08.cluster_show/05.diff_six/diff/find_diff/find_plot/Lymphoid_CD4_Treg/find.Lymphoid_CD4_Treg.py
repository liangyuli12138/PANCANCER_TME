import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_CD4_Treg");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_CD4.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["IKZF2","CARD16","TIGIT","CTLA4","IL32","FOXP3","TNFRSF9","BATF","TNFRSF18","DUSP4","TNFRSF1B","TBC1D4","STAM","DUSP16","RTKN2","IL2RA","LAYN","RGS1","NAMPT","CD27","PMAIP1","CTSC","MIR4435-2HG","HLA-A","IL2RB","ICOS","CYTOR","ARID5B","GBP5","CD74","LINC01943","GADD45A","GBP2","DNPH1","VAV3","HLA-DRB1","TNFRSF4","IKZF4","B2M","S100A4","HNRNPA1P21","PRDM1","LAIR2","CASP1","CCR8","PIM2","OAZ1","PHLDA1","NIBAN1","HTATIP2",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD4_Treg",
                  save="Lymphoid_CD4_Treg.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD4_Treg.markergene.0.xxx.pdf")

marker_genes_dict=["VDR","SPOCK2","PARK7","CCNG2","GLRX","ACOT9","SAT1","TYMP","AC104850.2","IL10RA","TTN","IL12RB2","UBE2B","IL21R","NOP58","PHTF2","IL1R2","FANK1","TRAC","GAPDH","CFLAR","LAPTM4B","IL2RG","GK","FAM53B","ZBTB38","BTG3","PKM","WHRN","UGP2","CALM3","CXCR6","F5","SAMSN1","CRADD","ANKRD10","BACH1","MAGEH1","PICALM","PELI1","HLA-DPB1","CACYBP","ENTPD1","MAST4","LINC02195","FOXO1","ARHGEF12","HLA-DPA1","PTPRJ","SNX9",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD4_Treg",
                  save="Lymphoid_CD4_Treg.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD4_Treg.markergene.1.xxx.pdf")
