import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Epithelium_Malig_Migration");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Epithelium.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["MALAT1","NFAT5","NEAT1","RNF213","ZSWIM6","ZFAND3","PARD3","SRRM2","FTX","FAM160A1","PATJ","TJP1","AC027644.4","SOX4","PTK2","ADAMTSL4-AS1","STAG1","ARIH1","STK3","LINC00511","MACF1","EXT1","KMT2C","FNDC3B","DLG1","BTBD9","FAF1","SSH2","CUX1","EXOC4","MECOM","SIK3","MYOF","EXOC6B","CHD9","LINC01876","MYH14","ANKRD11","SYNE2","DLEU1","SCAPER","MYO1E","BIRC6","RBM6","LAMA5","ERC1","LRBA","ASH1L","NAALADL2","ANKHD1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Malig_Migration",
                  save="Epithelium_Malig_Migration.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Malig_Migration.markergene.0.xxx.pdf")

marker_genes_dict=["ADK","KIF13A","LUC7L3","ANKRD36","GPHN","KIAA1217","MUC4","AKAP13","ATRX","RABGAP1L","MAPK8IP3","TRIO","GNGT1","WWOX","XPR1","SMYD3","ACAP2","SPIDR","USP34","ARHGAP26","VPS13B","SMURF1","JMJD1C","MYO1D","ARHGEF3","ZBTB20","HIST1H3J","BDP1","NF1","MAML2","WDR60","AOPEP","CTNNB1","SIPA1L3","RERE","NIPBL","PLEKHA5","ANKRD36B","UBR5","AC020916.1","PRRC2C","CAMK1D","RIPK1","ATP11A","SVIL","IMMP2L","NPAS2","LLGL2","ANGPTL4","BAIAP2L1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Malig_Migration",
                  save="Epithelium_Malig_Migration.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Malig_Migration.markergene.1.xxx.pdf")
