import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Fibroblast_iCAF_IGFBP6");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Fibroblast.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["CFD","DCN","SCARA5","IGFBP6","GSN","S100A4","FBLN2","CHRDL1","PLAC9","GPX3","MFAP5","FBLN1","C1R","ACKR3","CD34","GPNMB","MFAP4","CST3","CYBRD1","SLIT3","CCDC80","EXOSC7","S100A6","ADH1B","S100A10","CD248","PTGIS","TNXB","PI16","ANXA2","SERPING1","EFEMP1","OMD","C3","MGST1","GALNT15","CADM3","AHNAK","VAT1","PCOLCE2","CCN5","PSAP","PODN","CILP","OGN","PMP22","ABI3BP","ADAMTSL4","LTBP4","ADD3",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Fibroblast_iCAF_IGFBP6",
                  save="Fibroblast_iCAF_IGFBP6.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Fibroblast_iCAF_IGFBP6.markergene.0.xxx.pdf")

marker_genes_dict=["SEMA3C","F10","CLEC3B","SFRP2","ADAM33","TXNIP","CD55","FBN1","SH3BGRL3","TIMP2","PRELP","C1S","EBF1","SOD3","SERPINF1","ABLIM1","PIGT","FSTL1","VIT","C16orf89","UAP1","SLPI","CLU","CPVL","EMP3","MEDAG","PXN","TIMP3","APP","FCGRT","CRIP1","ITGBL1","TFPI","METTL7A","S100A13","LINC01133","DBN1","ANXA1","LRRN4CL","PCOLCE","FHL1","HSPB8","VEGFD","DDAH2","CTSK","GALNT12","NDRG1","NUPR1","SMIM14","CELF2",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Fibroblast_iCAF_IGFBP6",
                  save="Fibroblast_iCAF_IGFBP6.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Fibroblast_iCAF_IGFBP6.markergene.1.xxx.pdf")
