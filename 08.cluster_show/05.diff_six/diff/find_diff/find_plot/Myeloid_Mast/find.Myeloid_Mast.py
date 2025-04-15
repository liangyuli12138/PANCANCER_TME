import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Myeloid_Mast");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Myeloid.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["TPSAB1","TPSB2","MS4A2","CPA3","VWA5A","KIT","GATA2","AC092979.1","IL1RL1","HDC","CD82","HPGDS","CLU","LMNA","SLC18A2","TPSD1","MAOB","MYADM","SIGLEC6","PALM2-AKAP2","BACE2","AC109460.3","RHEX","ACSL4","LAPTM4A","RGS13","SIGLEC17P","ANXA1","IL18R1","PTGS1","FOSB","SLC2A3","CTSG","GCSAML","LAX1","SLC24A3","IDS","SAMSN1","RAC2","TIMP3","ERCC1","LMO4","MLPH","AL157895.1","LTC4S","TESPA1","YWHAZ","STX3","TNIK","RHOH",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_Mast",
                  save="Myeloid_Mast.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_Mast.markergene.0.xxx.pdf")

marker_genes_dict=["LEO1","FERMT2","VIM","GRAP2","RAB27B","CALB2","CAVIN1","DNAJA1","DLC1","AC008440.2","HS6ST1","CSF1","HSP90AA1","ANKRD28","ARHGEF6","CTNNBL1","SLC45A3","FOXP1","UBB","TUBA1A","CMA1","B4GALT5","SYTL2","MBOAT7","STXBP5","ENPP3","P2RX1","JUND","CD69","HPGD","SMYD3","ITM2C","GALC","SELENOK","SLC44A1","MORF4L2","RAB37","LIF","PAG1","SRGN","ATP6V0A2","CLEC4O","ITM2A","BATF","DUSP10","NDFIP2","ADRB2","G3BP2","PDZD8","SLC26A3",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_Mast",
                  save="Myeloid_Mast.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_Mast.markergene.1.xxx.pdf")
