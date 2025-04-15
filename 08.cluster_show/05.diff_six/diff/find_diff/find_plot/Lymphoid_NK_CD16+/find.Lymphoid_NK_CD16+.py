import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_NK_CD16+");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_NK.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["FCGR3A","FGFBP2","NKG7","SPON2","ADGRG1","PRF1","GZMH","GZMB","CX3CR1","MYOM2","CST7","PLEK","GNLY","ITGB2","S1PR5","EFHD2","PRSS23","PLAC8","CD247","ARL4C","LAIR2","CYBA","PTGDS","CXCR2","PXN","SYNE1","IGFBP7","AKR1C3","KIR2DL1","C1orf21","CHST2","TTC38","CEP78","LITAF","TBX21","CD300A","RASA3","AL590385.2","KLRF1","HOPX","GK5","SH3BP5","GPR141","ADD3","RASSF4","C12orf75","TGFBR3","MTSS1","OSBPL5","ZEB2",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_NK_CD16+",
                  save="Lymphoid_NK_CD16+.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_NK_CD16+.markergene.0.xxx.pdf")

marker_genes_dict=["TYROBP","BIN2","ITGAM","MLC1","KLRD1","SLCO4C1","CTBP2","FCRL6","PLEKHG3","SYNE2","SELPLG","KIR2DS4","RIPOR2","CTSW","LAIR1","SH2D1B","FGL2","IGF2R","PTPN12","TFDP2","GNG2","EMP3","LGR6","RAP1GAP2","ANKRD20A11P","GNPTAB","PYHIN1","UBE2F","DTHD1","SSBP3","LINC02384","RAB29","KLRC3","ARPC2","FGR","LGALS1","B3GAT1","UCP2","KLRC2","SLC15A4","CMKLR1","KIR3DL1","GOLM1","DAB2","GSAP","MBP","KIR2DL3","LGALS9B","F2R","SBK1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_NK_CD16+",
                  save="Lymphoid_NK_CD16+.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_NK_CD16+.markergene.1.xxx.pdf")
