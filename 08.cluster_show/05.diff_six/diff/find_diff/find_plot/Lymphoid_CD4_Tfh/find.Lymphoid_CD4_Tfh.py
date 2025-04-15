import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_CD4_Tfh");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_CD4.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["NR3C1","RNF19A","CD200","PPP1CC","ITM2A","MAF","FABP5","PDCD1","AHI1","DUSP4","DRAIC","IL6ST","COTL1","CXCL13","LIMS1","GNG4","SRGN","NMB","FKBP5","IGFBP4","TNFSF8","SH2D1A","RAB11FIP1","AC004585.1","TOX2","RBPJ","FAM107B","CNIH1","TOX","CD82","NAB1","PTPN11","IGFL2","GEM","CCDC50","REL","SPAG1","BTLA","IL21","CHN1","ARID5B","ICA1","GNAS","TSHZ2","PTPN13","OAZ1","ZNF281","YWHAQ","HSPE1","NFATC1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD4_Tfh",
                  save="Lymphoid_CD4_Tfh.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD4_Tfh.markergene.0.xxx.pdf")

marker_genes_dict=["MYO7A","RDH10","FZD3","SEC11A","PGM2L1","PGGHG","ST8SIA1","SMCO4","HSPA4","CPM","LRMP","IL6R","SIAH2","SGPP2","IQGAP1","THADA","CTSB","CHORDC1","AC008914.1","KIAA1324","RAB27A","HSPH1","NUDC","TP53INP1","ZEB2","PCNX2","TNFAIP8","HSPD1","TIGIT","DLEU2","FGFR1","TRIM8","FBLN7","HIF1A","ARMH1","AZIN1","NFATC2","MSI2","STX11","ETV6","AHSA1","ADD3","LINC01480","HBP1","MFHAS1","SESN3","AGFG1","CEMIP2","CNOT6L","GFOD1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD4_Tfh",
                  save="Lymphoid_CD4_Tfh.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD4_Tfh.markergene.1.xxx.pdf")
