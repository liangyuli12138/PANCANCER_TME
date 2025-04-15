import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_CD8_Tn");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_CD8.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["IL7R","CCR7","PABPC1","RPL9","RPL32","RPL34","RPS12","RPL13","SARAF","EEF1A1","RPS6","RPL3","EEF1G","RPS3A","CD55","RPS14","RPS13","RPS16","RPS4X","EZR","TPT1","RPS18","RPL30","RPL11","RPL5","RPS2","RPL10","RPL31","RPL4","TAGLN2","RPS23","RPL18","RPS3","RPL7","RPL21","RPL18A","RPL13A","RPS28","RPL27","RPL19","RPLP0","RPLP2","RPS27A","RPL14","RPS8","RPL35A","RPL36","RPS20","RPS21","RPL38",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD8_Tn",
                  save="Lymphoid_CD8_Tn.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD8_Tn.markergene.0.xxx.pdf")

marker_genes_dict=["RPS5","EEF2","AC024293.1","TCF7","RPS25","RPL8","NFKB1","RPS9","RACK1","RPS10-NUDT3","RPS10","RPL10A","RPL37","CHD2","RPSA","FAM107B","AC116533.1","CD44","RPL36A","RPS17","RPL23A","FTH1","GPR183","EEF1B2","RPLP1","SPOCK2","FAM177A1","RPS15A","RPS27","RPL37A","RPL26","RPL27A","NAP1L1","RPL28","RPL29","AC245033.4","AP001324.1","LMNA","BACH2","VIM","RPL23","RPS7","SMAP2","RPS29","SATB1","RPL17-C18orf32","NOP53","RPL7A","YPEL5","RPL12",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD8_Tn",
                  save="Lymphoid_CD8_Tn.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD8_Tn.markergene.1.xxx.pdf")
