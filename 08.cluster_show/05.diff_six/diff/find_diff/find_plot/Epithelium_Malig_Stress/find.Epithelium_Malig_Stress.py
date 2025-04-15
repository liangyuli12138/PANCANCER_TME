import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Epithelium_Malig_Stress");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Epithelium.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["FOSB","HSP90AB1","JUN","ERCC1","IER2","HSPA1B","FOS","DNAJA1","EGR1","PPP1R15A","CCNL1","CLK1","EIF1","DNAJB1","EIF5","UBC","MCL1","SERTAD1","SRSF2","JUNB","BRD2","SRSF3","H3F3B","EIF4A2","HSPA8","JUND","HNRNPDL","MT-CO1","PNRC1","HSP90AA1","ZFAND5","CSRNP1","IFRD1","FAU","KLF4","SFPQ","DDX5","TRIB1","HNRNPA2B1","LINC01578","HSPE1","RHOB","MT-CO3","HSPD1","CALM1","DUSP1","SRSF5","CEBPB","KLF10","UBB",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Malig_Stress",
                  save="Epithelium_Malig_Stress.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Malig_Stress.markergene.0.xxx.pdf")

marker_genes_dict=["CDKN1A","ATF4","HNRNPA1","DYNLL1","PTMA","SRSF7","MT-ND4","SENP3-EIF4A1","ATF3","PTGES3","RPS20","FUS","RPL37","NFKBIZ","MT-ATP6","DDX3X","TUBB4B","NR4A1","MAP1LC3B","TOB1","RPL24","MT-CYB","NPM1","HNRNPH3","RSRP1","SON","HSPA1A","RSRC2","MT-ND5","HNRNPH1","RPL15","YBX3","SBDS","ZFP36","GAS5","TCP1","TRA2B","MYLIP","GPX1","HNRNPA0","TRA2A","PPIA","HNRNPC","SNHG8","PCBP1","MORF4L2","MT-ND4L","ZFAS1","MIDN","IVNS1ABP",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Malig_Stress",
                  save="Epithelium_Malig_Stress.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Malig_Stress.markergene.1.xxx.pdf")
