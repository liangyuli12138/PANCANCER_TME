import scanpy as sc
import pandas as pd
import os
import anndata


os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/02.cluster_tls")
adata_concat = sc.read_h5ad("leiden_resolution_1.0.h5ad")

genes=[
"S100A8","S100A9","SLC11A1","SERPINA1","FCN1",
"CLEC9A","XCR1","WDFY4","SNX3","CPNE3",
"CD1C","FCER1A","CLEC10A","JAML","LST1",
"LAMP3","CCR7","FSCN1","CCL19","CCL17",
"TCF4","GZMB","LILRA4","PLAC8","IRF7",
"TPSAB1","TPSB2","MS4A2","CPA3","VWA5A","KIT",
"YBX3","TCL1A","IGHD","IL4R","FCER2","CD79A",
"GPR183","BCL2A1","CRIP1",
"IGHA1","IGHG1","MZB1","JCHAIN","IGKC",
"IGLC2","IGLC3","IGLC7",
"C1QC","C1QB","APOC1",
"SPP1","CCL4","CCL3L1","CXCR4",
"LYVE1","FOLR2","STAB1","F13A1",
"IL2RB","FCER1G","CD7","IRF8","LAT2","PIK3R1","GSTP1","CD160",
"FCGR3A","FGFBP2","NKG7","SPON2","ADGRG1","PRF1","GZMB","GNLY","CX3CR1","CST7",
"GZMH","CD3D","CD8A","CD8B","CD52","S100A10","CD3G","PATL2","ITGB1",
"SLC4A10","TRAC","DPP4","SLAMF1","DUSP1","CD28","CITED2","TNF","EGR1",
"IL4I1","LST1","KIT","FXYD5","NFKB1","AHR","NCOA7","AFF3","IL1R1",
"CCR7","LEF1","SELL","FTH1",
"GPR183","IL7R","ANXA1","ANKRD28",
"CXCL13","PDCD1","TOX2","IL21","GNG4","CD200",
"KLRB1","CCL20","RORA","AHR",
"FOXP3","IL2RA","CTLA4","TNFRSF4","IKZF2","TIGIT",
"GZMA","GZMK","CCL5",
"FOS","JUN","TNF","FOSB","JUNB","DUSP1",
"CCR7","TCF7","SELL","LEF1",
"GPR183","IL7R","ANXA1",
"DUSP4","CTLA4","TNFRSF9","RBPJ","RGS1","PDCD1","TIGIT","CD82","CXCL13","LAG3",
"GZMA","GZMH","NKG7",
"FOSB","ERCC1","JUN","EGR1","BRD2","IER2","JUNB","NR4A1",
"IFIT3","IFIT2","LIPA","HERC5","IFI44","DDX58","IFIH1",
]

sc.pl.umap(adata_concat, color=genes,
 s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save="st.tls.umap.markergene.all.png")

