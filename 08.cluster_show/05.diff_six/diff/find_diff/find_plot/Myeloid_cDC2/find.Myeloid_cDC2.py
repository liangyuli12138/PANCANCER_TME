import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Myeloid_cDC2");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Myeloid.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["RPS11","HLA-DPB1","RPS16","RPS2","RPS7","RPS14","RPS6","HLA-DQB1","RPS18","RPS3A","RPS19","RPS23","RPS13","RPS24","RPL9","HLA-DQA1","EEF1G","HLA-DRA","RPS15A","RPL30","RPL18","RPL4","HLA-DPA1","RPL18A","RPL27","RPL34","RPS27","RPL32","RPL37","RPS3","RPS29","RPS21","AC116533.1","RPL11","RPL6","RPS9","RPL31","RPL13A","RPL13","RPL10A","RPL10","RPL35A","RPSA","RPS12","HLA-DRB1","RPL7A","RPL28","RPL23A","RPL41","RPS17",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_cDC2",
                  save="Myeloid_cDC2.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_cDC2.markergene.0.xxx.pdf")

marker_genes_dict=["RPS4X","RPL26","RPL19","RPL36A","RPL15","CD1C","RPL5","RPL3","RACK1","RPS8","RPL8","RPS15","RPL27A","RPL37A","EEF1A1","RPL7","RPL23","GPR183","FAU","RPLP0","RPL17-C18orf32","PEA15","RPS10-NUDT3","CD74","AC024293.1","RPS27A","RPL35","UBA52","RPL38","RPLP1","PABPC1","EEF1B2","AC245033.4","RPS5","RPLP2","JAML","RPL29","RPL24","COTL1","RPL26P19","TPT1","RPL21","RPS25","LST1","RPL36","RPL14","NAPSB","NACA","RPS10","FCER1A",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_cDC2",
                  save="Myeloid_cDC2.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_cDC2.markergene.1.xxx.pdf")
