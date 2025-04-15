import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Myeloid_cDC1");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Myeloid.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["RPS11","RPS23","RPS16","RPL18A","NAPSB","RPS2","RPS21","RPL8","RPS7","RPS6","RPL4","RPL27A","RPS18","HLA-DPA1","HLA-DPB1","CD74","CPVL","RPL13A","RPS14","RPLP1","RPS3A","RPS13","EEF1B2","RPS12","RPL30","WDFY4","RPS29","ACTG1","RPL19","RPL32","LSP1","RPL18","RPS15A","RPL23A","RPSA","RPL31","C1orf54","RPL5","RPL14","RPS24","RPS27","RPL10","RPL6","RPL10A","RPS8","RPL23","NAP1L1","RPLP0","RPL21","RPL28",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_cDC1",
                  save="Myeloid_cDC1.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_cDC1.markergene.0.xxx.pdf")

marker_genes_dict=["RPL9","RPL27","RPS20","SNX3","RPL37","RPS3","FAU","RPL35","RPS25","IRF8","RPS15","LCP1","RPS17","RPS5","RPL7A","SNHG5","RPL35A","RPL41","RACK1","RPL11","PABPC1","RPL7","HLA-DQB1","RPL13","EEF1A1","RPS4X","HLA-DRA","AC116533.1","RPL26","EEF1G","RPL24","SUB1","NME1-NME2","HLA-DQA1","RPL36","RPS28","RPLP2","RPL29","RPL3","RPS27A","BASP1","RPL12","UBA52","RPL37A","SLAMF7","RPL34","NACA","AC245033.4","BTF3","CPNE3",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_cDC1",
                  save="Myeloid_cDC1.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_cDC1.markergene.1.xxx.pdf")
