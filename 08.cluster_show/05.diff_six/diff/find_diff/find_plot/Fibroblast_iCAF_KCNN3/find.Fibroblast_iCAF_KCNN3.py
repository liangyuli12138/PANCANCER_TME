import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Fibroblast_iCAF_KCNN3");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Fibroblast.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["C7","KCNN3","DPT","FXYD6","SCN7A","NDRG2","EGR1","MATN2","SPARCL1","FOSB","FOS","ERCC1","SRSF5","HEXIM1","JUN","P2RY1","FXYD1","ALDH1A1","LAMB1","PPP1R10","DDX5","RGS2","MBP","PID1","ABCA9","ADAMTSL3","NFIA","JUNB","DDX3X","CIRBP","RHOJ","DDIT3","JAM2","PDGFRA","PTN","RERG","SRSF3","ZNF106","CXCL12","SPTBN1","PTGER3","CD302","SLC8A1","IER2","PCF11","SLIT2","BRD2","HNRNPH1","SELENOP","MCL1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Fibroblast_iCAF_KCNN3",
                  save="Fibroblast_iCAF_KCNN3.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Fibroblast_iCAF_KCNN3.markergene.0.xxx.pdf")

marker_genes_dict=["PREX2","PBX1","SPRY1","RORA","DNAJB1","MGP","NFIB","ADH1B","H3F3B","IL6ST","ADCYAP1R1","RBBP6","C12orf60","MEIS1","SAP18","SRSF7","CCNL1","COLEC12","NR2F1","GUCY1A1","BOC","ARL6IP1","EIF5","HNRNPA2B1","ECHDC2","DNAJA1","GSTM5","TCEAL7","HOOK2","AC119674.2","SMARCA2","PLPP3","UBB","FAM110B","PPP1R15A","VAMP2","COL14A1","SON","LY6H","SVEP1","C6orf62","SFPQ","H2AFZ","WSB1","PHIP","CCL11","CALM1","AR","TENT5C","SGCD",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Fibroblast_iCAF_KCNN3",
                  save="Fibroblast_iCAF_KCNN3.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Fibroblast_iCAF_KCNN3.markergene.1.xxx.pdf")
