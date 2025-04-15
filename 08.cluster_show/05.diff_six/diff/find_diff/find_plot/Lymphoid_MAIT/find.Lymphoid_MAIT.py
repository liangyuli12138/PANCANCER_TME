import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_MAIT");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_NK.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["IL7R","SLC4A10","LTB","RPS12","CD3D","RPLP1","JAML","TRAC","JUNB","SPOCK2","RPL13","RPS20","TPT1","RPS29","RPL32","RPLP0","FOSB","RPL28","JUN","EGR1","RPLP2","EEF1A1","DPP4","RPL30","RPS25","KLRB1","RPL34","RPS2","EEF1G","S100A4","RPL9","Z94721.2","RPS13","CD6","SLAMF1","RPL11","FOS","RPS8","ERCC1","RPL36","RPL4","RPS21","AC116533.1","RPL19","RPS16","RPS18","RPL10","RPL41","TRBC2","RPL31",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_MAIT",
                  save="Lymphoid_MAIT.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_MAIT.markergene.0.xxx.pdf")

marker_genes_dict=["RPL18A","DUSP1","RPS28","SRSF7","RPL36A","RPL35A","RPL17-C18orf32","RPS6","LAPTM5","CD28","CD8A","RPL3","RPL21","AQP3","RPS27A","AP001324.1","RPL38","RPS14","EIF1","CITED2","RPSA","EEF1B2","RPL13A","RPS17","TNFRSF25","RPL23A","CAMK4","TNFAIP3","RPL8","RPS3A","RPS23","H3F3B","S100A6","PAG1","UBE2S","RPL12","RPL10A","RPS19","TNF","PBXIP1","RPS27","RGS2","RPL5","HOOK2","IL32","RPS3","AC245033.4","BTG2","HSPA8","RPL27A",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_MAIT",
                  save="Lymphoid_MAIT.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_MAIT.markergene.1.xxx.pdf")
