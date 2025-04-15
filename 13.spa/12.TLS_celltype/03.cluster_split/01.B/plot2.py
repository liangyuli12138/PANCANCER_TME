import scanpy as sc
import pandas as pd
import os
import anndata


os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/03.cluster_split/01.B")
adata_concat = sc.read_h5ad("leiden_resolution_1.0.h5ad")

marker_genes_dict = {
"Lymphoid_B_naive":["YBX3","TCL1A","IGHD","IL4R","FCER2","CD79A"],
"Lymphoid_B_memory":["GPR183","BCL2A1","CRIP1"],
"Lymphoid_Plamsa_IGKC":["IGHA1","IGHG1","MZB1","JCHAIN","IGKC"],
"Lymphoid_Plamsa_IGLC":["IGLC2","IGLC3","IGLC7"],
"Myeloid_Mono":["S100A8","S100A9","SLC11A1","SERPINA1","FCN1"],
"Myeloid_cDC1":["CLEC9A","XCR1","WDFY4","SNX3","CPNE3"],
"Myeloid_cDC2":["CD1C","FCER1A","CLEC10A","JAML","LST1"],
"Myeloid_cDC3":["LAMP3","CCR7","FSCN1","CCL19","CCL17"],
"Myeloid_pDC":["TCF4","GZMB","LILRA4","PLAC8","IRF7"],
"Myeloid_Mast":["TPSAB1","TPSB2","MS4A2","CPA3","VWA5A","KIT"],
"Lymphoid_NK_CD56+":["IL2RB","FCER1G","CD7","IRF8","LAT2","PIK3R1","GSTP1","CD160"],
"Lymphoid_NK_CD16+":["FCGR3A","FGFBP2","NKG7","SPON2","ADGRG1","PRF1","GZMB","GNLY","CX3CR1","CST7"],
"Lymphoid_NKT":["GZMH","CD3D","CD8A","CD8B","CD52","S100A10","CD3G","PATL2","ITGB1"],
"Lymphoid_MAIT":["SLC4A10","TRAC","DPP4","SLAMF1","DUSP1","CD28","CITED2","TNF","EGR1"],
"Lymphoid_ILC":["IL4I1","LST1","KIT","FXYD5","NFKB1","AHR","NCOA7","AFF3","IL1R1"],
"Lymphoid_CD4_Tn":["CCR7","LEF1","SELL","FTH1"],
"Lymphoid_CD4_Tcm":["GPR183","IL7R","ANXA1","ANKRD28"],
"Lymphoid_CD4_Tfh":["CXCL13","PDCD1","TOX2","IL21","GNG4","CD200"],
"Lymphoid_CD4_Th17":["KLRB1","CCL20","RORA","AHR"],
"Lymphoid_CD4_Treg":["FOXP3","IL2RA","CTLA4","TNFRSF4","IKZF2","TIGIT"],
"Lymphoid_CD4_CTL":["GZMA","GZMK","CCL5"],
"Lymphoid_CD4_Tstr":["FOS","JUN","TNF","FOSB","JUNB","DUSP1"],
"Lymphoid_CD8_Tn":["CCR7","TCF7","SELL","LEF1"],
"Lymphoid_CD8_Tm":["GPR183","IL7R","ANXA1"],
"Lymphoid_CD8_Tex":["DUSP4","CTLA4","TNFRSF9","RBPJ","RGS1","PDCD1","TIGIT","CD82","CXCL13","LAG3"],
"Lymphoid_CD8_Teff":["GZMA","GZMH","NKG7"],
"Lymphoid_CD8_Tstr":["FOSB","ERCC1","JUN","EGR1","BRD2","IER2","JUNB","NR4A1"],
"Lymphoid_CD8_Tisg":["IFIT3","IFIT2","LIPA","HERC5","IFI44","DDX58","IFIH1"],
}
import matplotlib.pyplot as plt

fig, axs = plt.subplots(ncols=4)

fig, axs = plt.subplots(ncols=4)

for ax, (cell_type, genes) in zip(axs.flat, marker_genes_dict.items()):
    sc.pl.umap(adata_concat, color=genes, s=2, frameon=False, vmax='p99.9', na_color='grey', 
               ax=ax, show=False)
    ax.set_title(f"{cell_type}: {genes[0]}")

plt.tight_layout()
plt.savefig("marker_genes_umap.png")
