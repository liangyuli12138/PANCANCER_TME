import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/EC_Lymph");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.EC.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["TFPI","MMRN1","TFF3","COLEC12","LAPTM5","CCL21","PKHD1L1","PROX1","NR2F1","CHRDL1","LYVE1","PPFIBP1","ADD3","MRC1","NRP2","PTPRE","LAMA4","SHC1","PARD6G","FLT4","PDPN","S100A10","MPP7","TSPAN5","MARCKS","PGM5","GYPC","ATP5F1E","SEMA3D","TBX1","AKAP12","DKK3","EPB41L2","SEMA3A","RELN","ANXA2","AL162231.2","MYH10","FABP5","OAF","CAVIN2","ECSCR","SCN3B","FABP4","RHOJ","STAB2","WFS1","STON2","SLC38A1","SLC24A1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Lymph",
                  save="EC_Lymph.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Lymph.markergene.0.xxx.pdf")

marker_genes_dict=["KBTBD11","SMYD2","CD9","LOX","DSP","EFNA5","AC007319.1","AC007998.3","MAF","FAM189A2","EFEMP1","FAM107B","STMN1","IGF1","LIMS1","RAB11FIP1","NR2F2","APOD","MYZAP","FAM102A","MRVI1-AS1","PTX3","ACTN1","TC2N","DOCK5","ABI3BP","POLR2L","AC142391.1","PLIN5","CTSZ","PLSCR4","CCDC80","PIEZO1","GGTA1P","HOMER3","PRXL2A","CR381670.1","SYNM","REEP1","RARRES2","TUBA1A","FXYD6","SCN3A","DST","GRAP","PDE1A","TNFAIP8L3","GUSB","SGCE","CD151",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="EC_Lymph",
                  save="EC_Lymph.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="EC_Lymph.markergene.1.xxx.pdf")
