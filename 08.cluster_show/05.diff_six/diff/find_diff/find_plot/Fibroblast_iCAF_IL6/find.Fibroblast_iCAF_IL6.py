import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Fibroblast_iCAF_IL6");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Fibroblast.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["CDKN1A","GEM","NAMPT","IL6","SOD2","SOCS3","SRPX","NNMT","GFPT2","C3","EIF4A3","ZFP36","CXCL12","FGF7","MAT2A","SLC2A3","CXCL2","GPRC5A","PNRC1","CTSL","KLF4","EMP1","MYC","MAFB","TNFAIP6","SENP3-EIF4A1","NFIL3","CEBPD","CEBPB","KLF9","ZFAS1","SERPING1","DDX5","JUNB","C1S","IFI16","RASD1","NFKBIA","MAP1LC3B","ELL2","EIF1","FTL","RND3","ANXA1","ZC3H12A","ARID5B","EFEMP1","H2AFZ","MCL1","APOE",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Fibroblast_iCAF_IL6",
                  save="Fibroblast_iCAF_IL6.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Fibroblast_iCAF_IL6.markergene.0.xxx.pdf")

marker_genes_dict=["CHD1","PLIN2","H3F3B","ATF3","IFRD1","MT2A","CCL2","DNAJB1","UAP1","CCNL1","DAB2","RSPO3","FOSB","PXDC1","SELENOP","CFD","MGP","GADD45B","IGFBP4","THBS1","STEAP1","IER3","ERRFI1","HLA-E","CSRNP1","VIM","NR4A1","DNAJA1","NABP1","SERTAD1","APOD","SLIT3","SFPQ","UBB","ID2","C11orf96","MAFF","PTGS2","CFH","MEDAG","YBX3","HSPB8","SLC19A2","PPP1R15A","ELOC","HAS2","MGST1","ADAMTS1","SERPINF1","HMOX1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Fibroblast_iCAF_IL6",
                  save="Fibroblast_iCAF_IL6.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Fibroblast_iCAF_IL6.markergene.1.xxx.pdf")
