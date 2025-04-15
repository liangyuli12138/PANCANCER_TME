import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Lymphoid_CD4_Tstr");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Lymphoid_CD4.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["FOS","GADD45B","IER2","JUN","FOSB","JUNB","BRD2","DUSP1","ERCC1","ID2","TNF","AL499604.1","AL691403.1","AC022217.3","EGR1","HEXIM1","HOOK2","HSPA1B","GADD45G","CITED2","ILF3-DT","PPP1R10","SNHG5","RASD1","HIST1H2BG","BTG2","DDIT3","SNAI1","IER3","NR4A1","PLK2","PRR3","NFKBIA","AL671762.1","RPS12","TRMO","DDIT4","HIST1H3J","H2AFX","DNAJB1","SRSF3","SERTAD1","CSKMT","ZFP36","RGS16","RPL13P5","SGK1","RPLP2","AP001160.1","AL592429.2",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD4_Tstr",
                  save="Lymphoid_CD4_Tstr.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD4_Tstr.markergene.0.xxx.pdf")

marker_genes_dict=["TRA2B","ZFAS1","AL355075.4","ARL4D","UBE2S","MYL12A","AL118516.1","ARRDC3","AC023157.3","RPS29","PEBP1P3","EEF1B2","SCML4","RPL30","TEX14","EIF4A2","RPL34","HIST1H2BD","NXF1","CD48","AP001324.1","RPL27A","SAP18","KLF2","RPS27A","RPL31","RPS21","GIMAP1","AC109326.1","AC011446.2","SNHG8","RPL38","RPS26P43","TOB1","AC116533.1","GIMAP4","POLG2","CD52","EIF4A1","RN7SK","RPL32","RPL36A","IL4I1","RPL23A","RPS27","SAE1","RPS18","RHOB","RPS20","AC011452.1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Lymphoid_CD4_Tstr",
                  save="Lymphoid_CD4_Tstr.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Lymphoid_CD4_Tstr.markergene.1.xxx.pdf")
