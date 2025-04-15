import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Mural_cell_Pericyte2");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Mural_cell.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["CD44","IGFBP5","GGT5","CCL2","RARRES2","SYNPO2","CLSTN2","STEAP4","SOCS3","C1R","CXCL12","IFITM2","C1S","TIMP3","AC003991.1","PRRX1","SERPING1","CFD","MYC","CEBPD","COL14A1","TM4SF1","TMEM176B","ECRG4","TMEM176A","FGF7","EXOSC7","LGI4","GPX3","CTSC","SYN3","TGFBI","NR2F2","PROCR","NAMPT","CDKN1A","THBS1","SSTR2","NNMT","EMP1","ZFP36","EPS8","ZEB2","RPS2","CHRDL1","C11orf96","EGR1","AC007563.2","LINC01197","ADH1B",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Mural_cell_Pericyte2",
                  save="Mural_cell_Pericyte2.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Mural_cell_Pericyte2.markergene.0.xxx.pdf")

marker_genes_dict=["FTL","ADAMTS4","MT2A","CFH","SYTL2","EPB41L2","ARHGAP15","APOE","IFITM1","ID4","ARID5B","CSF1","SENP3-EIF4A1","SOD2","GAS5","FMO2","S100A4","NFIL3","CRISPLD2","RGS16","SLC2A3","RPL13","NFIB","LHFPL6","PROS1","PMEPA1","EPHX1","ABCC9","JUNB","PPP1R15A","HSPA1A","GEM","LMCD1","FSIP1","OSMR","RPL12","ZFP36L2","RPS18","SLC38A2","SLC1A5","H3F3B","IFI16","MYH9","CP","CCL19","BRD2","FGF10","NOP53","BTG1","RPL11",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Mural_cell_Pericyte2",
                  save="Mural_cell_Pericyte2.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Mural_cell_Pericyte2.markergene.1.xxx.pdf")
