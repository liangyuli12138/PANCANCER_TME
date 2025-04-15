import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Mural_cell_Pericyte1");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Mural_cell.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["PDGFRB","MYO1B","ARHGDIB","HIGD1B","NOTCH3","ADGRF5","RGS5","CYGB","NDUFA4L2","SPARC","PLXDC1","DLC1","ITGA1","GJC1","COL4A2","CCDC102B","LHFPL6","NID1","ETS1","COL18A1","COL4A1","SEPTIN11","COL5A3","LPL","UACA","COX4I2","EGFLAM","THY1","FRMD3","EPS8","PHLDA1","MIR4435-2HG","SPTBN1","PAG1","GMFG","ARHGEF17","CAMK2N1","CD248","HLA-B","GUCY1A2","IFITM2","EPAS1","SEPTIN4","CD36","IFI27","COL5A2","CYTOR","TFPI","MYH9","SEMA5A",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Mural_cell_Pericyte1",
                  save="Mural_cell_Pericyte1.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Mural_cell_Pericyte1.markergene.0.xxx.pdf")

marker_genes_dict=["CHN1","SPRY4","HLA-A","RFTN1","KCNJ8","NRP1","TNFRSF21","ENPEP","EBF1","EPB41L1","IGFBP7","CHST2","TRIB2","NEURL1B","BGN","FJX1","ANO1","NCK2","TRPC6","FAM162B","ECE1","F2R","TMEM204","ANGPT2","COL3A1","ABCC9","PMEPA1","HLA-C","STEAP4","C1QTNF5","ADCY3","TCF4","TBX2","LZTS1","GUCY1B1","PELO","SH2D3C","ADAP2","SPRY1","ITGA4","PEAK1","SH2B3","CSPG4","CDH6","ASAP1","DPYSL2","ECM1","CPE","PCOLCE","CTSC",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Mural_cell_Pericyte1",
                  save="Mural_cell_Pericyte1.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Mural_cell_Pericyte1.markergene.1.xxx.pdf")
