import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_ILC");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_NK.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["IL4I1","LST1","KIT","IL7R","TPT1","FXYD5","NFKB1","TCF7","AHR","EEF1G","NAP1L1","NCOA7","LINC00299","RPLP1","IL1R1","TAGLN2","AFF3","TNFRSF18","H2AFY","PABPC1","CDKN1A","TMEM123","TNFRSF25","EEF1A1","CD44","IRF2BP2","AQP3","DLL1","FURIN","ELOVL5","PCDH9","SVIL","RPL8","RPS9","ZFP36L1","LTB","RPL13","JMY","GOLGA8A","RPS6KA3","SKIL","PRNP","EML4","CHKA","MBOAT7","GNA15","MAP3K8","AL136456.1","PRMT9","EEF2",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_ILC",
                  save="Lymphoid_ILC.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_ILC.markergene.0.xxx.pdf")

marker_genes_dict=["NINJ1","CAT","KMT2E","NRIP1","CERK","KLF4","RPL9","SPINK2","CXCR4","TOX2","ELL2","GOLGA8B","DENND4A","MGAT5","SSBP2","LTA4H","RORC","PIM3","LIF","RPL11","IL23R","SCN1B","RPS24","SPOCK2","IL4R","FES","TNFSF13B","ESYT2","BACH2","KRT86","IFFO2","GPR183","PER1","GNAS","DST","SNRNP70","GSN","GDE1","CMTM6","AHI1","TNFSF11","HNRNPA0","TNFAIP3","SRGN","VIM","IGFBP4","EIF4G2","TTN","PLEKHO1","PAK1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_ILC",
                  save="Lymphoid_ILC.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_ILC.markergene.1.xxx.pdf")
