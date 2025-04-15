import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Epithelium_Malig_Cycle");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Epithelium.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["TUBA1B","SNRPG","TUBB","RAN","PPIA","STMN1","H2AFZ","SNRPD1","HMGN2","RANBP1","PTMA","SEM1","NME1","SNRPF","SNRPD2","H2AFV","HMGB1","DEK","GAPDH","HMGB2","MRPL51","HNRNPA3","PA2G4","AP000350.4","NME1-NME2","UQCRH","CFL1","RPS2","TMA7","MCM7","LSM3","ANP32B","LSM4","ENO1","UQCRQ","YBX1","ATP5MF","SNHG6","TPI1","PCNA","PFN1","SNRPE","NUCKS1","TYMS","PSME2","SNRPB","PDCD5","CCT5","SLC25A5","PSMA7",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Malig_Cycle",
                  save="Epithelium_Malig_Cycle.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Malig_Cycle.markergene.0.xxx.pdf")

marker_genes_dict=["RPS19","RPSA","TRMT112","PCLAF","HINT1","CBX3","NDUFS6","ERH","TK1","ATP5F1E","NDUFS5","CHCHD2","TAGLN2","SLIRP","RPS20","H3F3A","MCM4","SET","UBE2T","ANAPC11","TXNDC17","SSBP1","RPLP0","CKS1B","H2AFY","SUMO2","TMSB10","DUT","RPL35","DNAJC9","DHFR","TOMM5","ZWINT","RPS7","RPS16","PSMA4","LSM5","GGCT","MAD2L1","RPA3","DYNLL1","ELOB","UBA52","COX6B1","AL365205.1","HNRNPA2B1","ATP5PF","SMC4","HMGA1","POLR2G",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Epithelium_Malig_Cycle",
                  save="Epithelium_Malig_Cycle.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Epithelium_Malig_Cycle.markergene.1.xxx.pdf")
