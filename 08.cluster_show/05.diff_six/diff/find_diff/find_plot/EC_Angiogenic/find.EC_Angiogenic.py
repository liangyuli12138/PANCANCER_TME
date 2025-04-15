import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/EC_Angiogenic");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.EC.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["SPARC","INSR","HSPG2","COL4A1","LAMB1","COL4A2","VWA1","PLVAP","IGFBP7","RGCC","LAMA4","MIR4435-2HG","IVNS1ABP","GRB10","FLT1","NOTCH4","MCAM","HTRA1","ESAM","ITGB1","CD93","RFLNB","TCF4","ETS1","ARHGAP18","SPRY4","EFNA1","PLXND1","JUP","MLEC","EDNRB","GNAS","ESM1","ARHGAP29","KDR","APLP2","GSN","A2M","VWF","RGS3","ITGA1","MYO1B","ARHGDIB","CXCR4","FYN","CD34","CYTOR","COL18A1","TP53I11","AFAP1L1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Angiogenic",
                  save="EC_Angiogenic.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Angiogenic.markergene.0.xxx.pdf")

marker_genes_dict=["F2R","TNFAIP8L1","PODXL","DLL4","PLXNA2","COL15A1","DGKD","SEMA3F","CTNNB1","FSCN1","UNC5B","CDK2AP1","C1QTNF5","NID1","YWHAH","PXDN","CLEC14A","MYO10","GNAI2","CA2","PDGFD","LBH","ADGRL4","OAZ2","NDST1","SERPINH1","DYSF","MYO6","LAMC1","KCNE3","PDGFB","VASH1","NRP1","SEC14L1","CALM1","THY1","CALCRL","RASGRP3","INHBB","ENG","EHD4","TM4SF18","TMEM204","GJA1","CDH5","SIPA1L2","GPR4","ANGPT2","CFL1","APP",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Angiogenic",
                  save="EC_Angiogenic.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Angiogenic.markergene.1.xxx.pdf")
