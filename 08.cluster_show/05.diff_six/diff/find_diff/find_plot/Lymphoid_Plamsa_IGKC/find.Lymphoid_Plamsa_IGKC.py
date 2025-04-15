import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_Plamsa_IGKC");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_B.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["IGKC","ANKRD36BP2","SSR4","RRBP1","FKBP2","AC073610.2","MZB1","FNDC3B","FKBP11","HDLBP","P4HB","JUN","SPATS2","MALAT1","SEC61B","XBP1","CYBA","DNAJC1","HM13","SEC11C","SPCS2","PRDX4","SELENOS","SLAMF7","FBXW7","SDF2L1","SSR3","TMED9","SND1","FNDC3A","SEC61A1","TENT5C","TXNDC5","HSPA1B","HERPUD1","RABAC1","SEL1L","TXNDC11","SDC1","TMEM258","UBE2J1","SPCS3","FER1L4","PSAP","PPIB","SPCS1","NEAT1","HSP90B1","ITM2C","ACADVL",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_Plamsa_IGKC",
                  save="Lymphoid_Plamsa_IGKC.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_Plamsa_IGKC.markergene.0.xxx.pdf")

marker_genes_dict=["RPN2","PDIA4","MAN1A2","ANKRD28","CREB3L2","MYDGF","LINGO1","PDK1","MYO1D","CD38","CLIC4","PDIA6","CLPTM1L","LINC02362","LARP1B","DPP7","LMAN2","LMAN1","DERL3","JSRP1","IFNAR2","SEC61G","GSTP1","SLC44A1","PRDM1","NDUFA1","EDF1","NUCB2","CPEB4","IGHG1","TRAM1","FCRL5","KRTCAP2","ERC1","TMEM59","SEC14L1","SIL1","MAN1A1","AC244205.1","CRELD2","POU2AF1","OS9","SEC24D","ERGIC3","USO1","HSH2D","H1FX","CYTOR","MANF","IGHG3",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_Plamsa_IGKC",
                  save="Lymphoid_Plamsa_IGKC.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_Plamsa_IGKC.markergene.1.xxx.pdf")
