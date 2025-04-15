import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Myeloid_Marco_LYVE1");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Myeloid.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["STAB1","RNASE1","FOLR2","ITSN1","ZSWIM6","F13A1","SELENOP","RBM47","DOCK4","AC016831.7","NEAT1","ITPR2","PEAK1","FRMD4B","SLC9A9","MAF","DYNC1H1","LGMN","ARHGAP26","ABCA1","NRP1","DAB2","FRMD4A","STARD13","MCTP1","SIK3","SIPA1L1","LILRB5","NRP2","QKI","FAM20A","AC009690.1","AP2A2","MALAT1","CD163L1","LYVE1","SEMA4A","CEMIP2","RGL1","SLC16A10","SRGAP1","CCDC152","SNX9","USP36","ADAP2","MTSS1","FMN1","ME1","TCF12","WWP1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_Marco_LYVE1",
                  save="Myeloid_Marco_LYVE1.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_Marco_LYVE1.markergene.0.xxx.pdf")

marker_genes_dict=["HIF1A","A2M","ATG7","TTN","FMNL2","MYO5A","SLC4A7","KIF1B","CD209","FCHSD2","SIGLEC1","SNX29","SAMD4A","PAPSS2","MGAT5","SPAG9","RUNX1","WLS","IGF1","CMAHP","AC017083.3","RAPH1","MERTK","PELI1","NPHS1","ABCC5","NPL","AC122718.1","GAS6","FOXO3","LPAR6","ZHX2","PLXDC1","SLC40A1","MEF2C","AFF4","AC069368.1","AFF1","GGTA1P","AC027644.4","LRMDA","FCHO2","PDGFC","SFMBT2","FHIT","ANKH","ELMO1","USP37","MACF1","MAML2",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_Marco_LYVE1",
                  save="Myeloid_Marco_LYVE1.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_Marco_LYVE1.markergene.1.xxx.pdf")
