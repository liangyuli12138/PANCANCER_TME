import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib as mpl

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/04.marker_plot")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Myeloid_Mono.0905.h5ad")
marker_genes_dict = {
"Myeloid_Mono":["S100A8","S100A9","SLC11A1","SERPINA1","FCN1"],
"Myeloid_cDC1":["CLEC9A","XCR1","WDFY4","SNX3","CPNE3"],
"Myeloid_cDC2":["CD1C","FCER1A","CLEC10A","JAML","LST1"],
"Myeloid_cDC3":["LAMP3","CCR7","FSCN1","CCL19","CCL17"],
"Myeloid_pDC":["TCF4","GZMB","LILRA4","PLAC8","IRF7"],
"Myeloid_Mast":["TPSAB1","TPSB2","MS4A2","CPA3","VWA5A","KIT"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.8, categories_order=["Myeloid_Mono","Myeloid_cDC1","Myeloid_cDC2","Myeloid_cDC3","Myeloid_pDC","Myeloid_Mast",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Myeloid_Mono.pdf")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Lymphoid_B.0905.h5ad")
marker_genes_dict = {
"Lymphoid_B_naive":["YBX3","TCL1A","IGHD","IL4R","FCER2","CD79A"],
"Lymphoid_B_memory":["GPR183","BCL2A1","CRIP1"],
"Lymphoid_Plamsa_IGKC":["IGHA1","IGHG1","MZB1","JCHAIN","IGKC"],
"Lymphoid_Plamsa_IGLC":["IGLC2","IGLC3","IGLC7"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.8, categories_order=["Lymphoid_B_naive","Lymphoid_B_memory","Lymphoid_Plamsa_IGKC","Lymphoid_Plamsa_IGLC",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Lymphoid_B.pdf")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.EC.0905.h5ad")
marker_genes_dict = {
"EC_Angiogenic":["INSR","COL4A1","LAMB1","COL4A2","VWA1","MCAM"],
"EC_Artery":["SEMA3G","IGF2","AIF1L","CXCL12","GLUL","SRGN","ASS1","KCTD12","FBLN2"],
"EC_Capillary":["CA4","SLCO2A1","RGCC"],
"EC_Vein":["ACKR1","CLU","SELE","SELP","OLFM1"],
"EC_Lymph":["MMRN1","TFF3","COLEC12","CCL21","PROX1","LYVE1"],
"EC_Alveolar":["FENDRR","VIPR1","SLC6A4","HPGD","TMEM100","TNFSF10"],
"EC_Glomerular":["IGFBP5","PLAT","EHD3","SOST","SYNE1","TGFBR2","ITGA8","SLC14A1"],
"EC_Sinusoidal":["CTSD","CLEC4G","CTSL","FCN2","LGMN","OIT3","ACP5"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.8, categories_order=["EC_Angiogenic","EC_Artery","EC_Capillary","EC_Vein","EC_Lymph","EC_Alveolar","EC_Glomerular","EC_Sinusoidal",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.EC.pdf")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Mural_cell.0905.h5ad")
marker_genes_dict = {
"Mural_cell_Pericyte1":["PDGFRB","MYO1B","ARHGDIB","HIGD1B","NOTCH3","RGS5","CYGB","SPARC","NDUFA4L2","COL4A1"],
"Mural_cell_Pericyte2":["CD44","IGFBP5","GGT5","CCL2","CLSTN2","CXCL12"],
"Mural_cell_SMC1":["ADIRF","MUSTN1","RERGL","PLN","SORBS2","CRIP1","NET1"],
"Mural_cell_SMC2":["RAMP1","PALLD","MYLK","TCEAL4","NCAM1","COLEC12","PSAP","SFRP4"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.8, categories_order=["Mural_cell_Pericyte1","Mural_cell_Pericyte2","Mural_cell_SMC1","Mural_cell_SMC2",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Mural_cell.pdf")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Myeloid_Marco.0905.h5ad")
marker_genes_dict = {
"Myeloid_Marco_C1QC":["C1QC","C1QB","APOC1"],
"Myeloid_Marco_SPP1":["SPP1","CCL4","CCL3L1","CXCR4"],
"Myeloid_Marco_LYVE1":["LYVE1","FOLR2","STAB1","F13A1"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.8, categories_order=["Myeloid_Marco_C1QC","Myeloid_Marco_SPP1","Myeloid_Marco_LYVE1",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Myeloid_Marco.pdf")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Lymphoid_NK.0905.h5ad")
marker_genes_dict = {
"Lymphoid_NK_CD56+":["IL2RB","FCER1G","CD7","IRF8","LAT2","PIK3R1","GSTP1","CD160"],
"Lymphoid_NK_CD16+":["FCGR3A","FGFBP2","NKG7","SPON2","ADGRG1","PRF1","GZMB","GNLY","CX3CR1","CST7"],
"Lymphoid_NKT":["GZMH","CD3D","CD8A","CD8B","CD52","S100A10","CD3G","PATL2","ITGB1"],
"Lymphoid_MAIT":["SLC4A10","TRAC","DPP4","SLAMF1","DUSP1","CD28","CITED2","TNF","EGR1"],
"Lymphoid_ILC":["IL4I1","LST1","KIT","FXYD5","NFKB1","AHR","NCOA7","AFF3","IL1R1"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.8, categories_order=["Lymphoid_NK_CD56+","Lymphoid_NK_CD16+","Lymphoid_NKT","Lymphoid_MAIT","Lymphoid_ILC",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Lymphoid_NK.pdf")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Fibroblast.0905.h5ad")
marker_genes_dict = {
"Fibroblast_mCAF_POSTN":["BGN","SPARC","POSTN","CTHRC1","COL1A1","COL3A1","COL8A1","VCAN"],
"Fibroblast_mCAF_KRT19":["MYLK","KRT19","TCEAL2","TCEAL4"],
"Fibroblast_mCAF_WNT5A":["TBX3","WNT5A","FHL2","TGFBI"],
"Fibroblast_iCAF_IGFBP6":["IGFBP6","SCARA5","FBLN2","MFAP5"],
"Fibroblast_iCAF_IL6":["IL6","GEM","NAMPT","SOD2","SOCS3","CXCL2","SLC2A3"],
"Fibroblast_iCAF_KCNN3":["KCNN3","FXYD6","SCN7A","NDRG2","MATN2","MBP","SPARCL1"],
"Fibroblast_apCAF_CD74":["CD74","HLA-DRB1","HLA-DRB5","HLA-DRA"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.8, categories_order=["Fibroblast_mCAF_POSTN","Fibroblast_mCAF_KRT19","Fibroblast_mCAF_WNT5A","Fibroblast_iCAF_IGFBP6","Fibroblast_iCAF_IL6","Fibroblast_iCAF_KCNN3","Fibroblast_apCAF_CD74",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Fibroblast.pdf")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Epithelium.0905.h5ad")
marker_genes_dict = {
"Epithelium_Malig_Migration":["MALAT1","NEAT1","NFAT5","ZSWIM6","FTX","SRRM2","SOX4","STK3"],
"Epithelium_Malig_Cycle":["TUBA1B","SNRPG","TUBB","STMN1","H2AFZ","SNRPD1","HMGN2","HMGB1"],
"Epithelium_Malig_cEMT":["TPM2","CAVIN1","COL17A1","CALD1","MYL9","ACTA2","MYLK","ACTN1","SPARCL1","TAGLN"],
"Epithelium_Malig_Interferon":["B2M","STAT1","UBE2L6","PSMB8","HLA-E","XAF1","HLA-A","IFI6","CD74"],
"Epithelium_Malig_Stress":["FOSB","JUN","ERCC1","IER2","HSPA1B","EGR1","HSP90AB1","DNAJA1","CLK1","MCL1"],
"Epithelium_Malig_Basal":["ANXA2","YWHAZ","KRT19","S100A6","KRT18"],
"Epithelium_Malig_Glandular":["TMED3","RABAC1","TTC3","TFF3","AGR2","WFDC2","CKAP4","MUC2","CKAP4"],
"Epithelium_Normal":["MYBPC1","MGP","AZGP1","TAT","SFTPC","AGR3","PHGR1","ERO1B","PIP","GFRA1"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.8, categories_order=["Epithelium_Malig_Migration","Epithelium_Malig_Cycle","Epithelium_Malig_cEMT","Epithelium_Malig_Interferon","Epithelium_Malig_Stress","Epithelium_Malig_Basal","Epithelium_Malig_Glandular","Epithelium_Normal",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Epithelium.pdf")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Lymphoid_CD4.0905.h5ad")
marker_genes_dict = {
"Lymphoid_CD4_Tn":["CCR7","LEF1","SELL","FTH1"],
"Lymphoid_CD4_Tcm":["GPR183","IL7R","ANXA1","ANKRD28"],
"Lymphoid_CD4_Tfh":["CXCL13","PDCD1","TOX2","IL21","GNG4","CD200"],
"Lymphoid_CD4_Th17":["IL17A","KLRB1","CCL20","RORA","AHR"],
"Lymphoid_CD4_Treg":["FOXP3","IL2RA","CTLA4","TNFRSF4","IKZF2","TIGIT"],
"Lymphoid_CD4_CTL":["GZMA","GZMK","CCL5"],
"Lymphoid_CD4_Tstr":["FOS","JUN","TNF","FOSB","JUNB","DUSP1"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.8, categories_order=["Lymphoid_CD4_Tn","Lymphoid_CD4_Tcm","Lymphoid_CD4_Tfh","Lymphoid_CD4_Th17","Lymphoid_CD4_Treg","Lymphoid_CD4_CTL","Lymphoid_CD4_Tstr",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Lymphoid_CD4.pdf")


adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.Lymphoid_CD8.0905.h5ad")
marker_genes_dict = {
"Lymphoid_CD8_Tn":["CCR7","TCF7","SELL","LEF1"],
"Lymphoid_CD8_Tm":["GPR183","IL7R","ANXA1"],
"Lymphoid_CD8_Tex":["DUSP4","CTLA4","TNFRSF9","RBPJ","RGS1","PDCD1","TIGIT","CD82","CXCL13","LAG3"],
"Lymphoid_CD8_Teff":["GZMA","GZMH","NKG7"],
"Lymphoid_CD8_Tstr":["FOSB","ERCC1","JUN","EGR1","BRD2","IER2","JUNB","NR4A1"],
"Lymphoid_CD8_Tisg":["IFIT3","IFIT2","LIPA","HERC5","IFI44","DDX58","IFIH1"],
}
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri",dot_max=0.8, categories_order=["Lymphoid_CD8_Tn","Lymphoid_CD8_Tm","Lymphoid_CD8_Tex","Lymphoid_CD8_Teff","Lymphoid_CD8_Tstr","Lymphoid_CD8_Tisg",],use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".sub.cluster.markergene.Lymphoid_CD8.pdf")


