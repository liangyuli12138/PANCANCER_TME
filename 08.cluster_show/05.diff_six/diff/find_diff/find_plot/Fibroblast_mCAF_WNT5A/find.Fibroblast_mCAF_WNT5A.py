import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Fibroblast_mCAF_WNT5A");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Fibroblast.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["FHL2","TGFBI","TBX3","NSG1","PRDM1","WNT5A","FBLIM1","ETS2","TRPA1","KLHDC8A","RBP1","IGFBP2","CDH11","MARCKSL1","MME","HSD17B2","PEG3","NRG1","FOXL2","GOLGA8B","ATP1B1","TSPAN5","COL7A1","SEMA5A","STAR","TCIM","PPP1R14A","BMP5","KRT18","LAMC3","SOX6","BMP7","RNASE1","MRVI1","DACH1","FARP1","LDLR","MYO10","PAG1","ITGA8","WNT4","CHCHD10","PCSK6","GATA4","MEST","WFDC1","DOCK5","PTCH1","NAP1L1","PMAIP1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Fibroblast_mCAF_WNT5A",
                  save="Fibroblast_mCAF_WNT5A.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Fibroblast_mCAF_WNT5A.markergene.0.xxx.pdf")

marker_genes_dict=["PWWP3B","MSI2","TSPAN33","ATRNL1","SIGLEC11","EMID1","COL27A1","P3H2","BST2","TANC1","DHRS2","FRY","PSD3","NR2F2","TBX2","EDNRA","NDNF","GOLGA8A","COLEC11","GREB1","ANK3","RYR2","SPINT2","PDGFC","GATM","NRP2","SLC14A1","MGAT5","LOXL2","ENC1","SEZ6L2","ADAMTS9","STMN1","FENDRR","MYLK","DMD","TSPAN6","ALCAM","SCD","PTGER2","TM4SF1","HSPA1B","SYNPO2","ALKAL2","MDK","ARX","COL4A5","TFAP2C","PXDN","ELK1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Fibroblast_mCAF_WNT5A",
                  save="Fibroblast_mCAF_WNT5A.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Fibroblast_mCAF_WNT5A.markergene.1.xxx.pdf")
