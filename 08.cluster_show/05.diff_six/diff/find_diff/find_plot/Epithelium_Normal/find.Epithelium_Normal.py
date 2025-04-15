import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Epithelium_Normal");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Epithelium.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["MYBPC1","PIP","AZGP1","SLC39A6","ELOVL5","TAT","GFRA1","SFTPC","AGR3","TAT-AS1","PIEZO2","GALNT5","TBC1D9","SFTPA1","SCUBE2","SYCE3","SNHG11","SERPINA3","NINJ1","PHGR1","SFTPA2","PPP4R4","NEK10","B4GALT1","FAM102A","LPP","DEGS2","COX7A1","PI15","SFTPD","EVL","PCAT18","ACAT1","TGFBR3","AC103957.2","SCUBE1","TMC3","IL6ST","CLIC6","CYP4F22","ENTPD5","AL121790.1","LAMP3","RHOH","TTN","ERO1B","LGALS4","CHAD","SH3BGRL","CRIP1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Normal",
                  save="Epithelium_Normal.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Normal.markergene.0.xxx.pdf")

marker_genes_dict=["AC106045.1","PLAT","MAGED2","PREX1","ITPR1","HAMP","CLGN","GHRH","CA12","IGHA2","AC110619.1","ALDOB","CLDN18","UCP1","LINC00504","FAM3D","ESR1","WIF1","KIAA1324","ABCA3","BATF","MRPS30-DT","CDX1","H2AFJ","CDON","SLC1A4","PEBP4","SORD","TGFB2","GALNT7","ADH1C","CGNL1","MGP","TPRG1","SLC1A1","DLC1","PLIN5","KIF5C","C15orf48","SERPINA11","ARHGEF6","TPPP3","SULF2","PGLYRP2","HIVEP3","SLC6A4","SLC16A6","LPL","CACNA2D2","C12orf60",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Normal",
                  save="Epithelium_Normal.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Normal.markergene.1.xxx.pdf")
