import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_Plamsa_IGLC");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_B.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["SSR4","IGLC2","FKBP2","AC073610.2","FKBP11","MZB1","SEC11C","SLC26A3","RRBP1","SPCS2","CYBA","SEC61B","HM13","SELENOS","RABAC1","SSR3","TMED9","P4HB","PRDX4","TMEM258","PSAP","RPN2","DERL3","TMEM59","SPCS1","HDLBP","MYDGF","SLAMF7","SEC61A1","HSPA1B","SDF2L1","PPIB","KRTCAP2","XBP1","MANF","IGLC3","TXNDC5","HERPUD1","JUN","SSR2","TRAM1","HSP90B1","NDUFA1","LMAN2","RPL36AL","NUCB2","SPCS3","LMAN1","JCHAIN","PDIA6",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_Plamsa_IGLC",
                  save="Lymphoid_Plamsa_IGLC.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_Plamsa_IGLC.markergene.0.xxx.pdf")

marker_genes_dict=["CD63","MAN1A2","SEL1L","UBE2J1","ERLEC1","ACADVL","UBC","SEC61G","HLA-C","ITM2C","SPATS2","SDC1","EDF1","RPN1","PDIA4","TENT5C","GSTP1","ERGIC3","SELENOK","SERP1","HSPA1A","CANX","OST4","FNDC3B","NPC2","TMED10","OS9","SERF2","CRELD2","IGLL5","FBXW7","SRPRB","CUTA","DNAJC1","AC009570.2","FTL","SRPRA","CLPTM1L","DAD1","UQCRQ","TXNDC15","TMEM208","IGHA1","LINGO1","KDELR2","TXNDC11","NME1-NME2","SEC62","AL121845.3","AP000350.4",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_Plamsa_IGLC",
                  save="Lymphoid_Plamsa_IGLC.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_Plamsa_IGLC.markergene.1.xxx.pdf")
