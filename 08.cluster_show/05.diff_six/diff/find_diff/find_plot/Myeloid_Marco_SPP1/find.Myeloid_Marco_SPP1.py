import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/05.diff_six/diff/find_diff/find_plot/Myeloid_Marco_SPP1");
adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/08.cluster_show/01.all.cluster/split_diff/pancancer.split.Myeloid.0831.h5ad")
sc.set_figure_params(figsize=(15, 15))


marker_genes_dict=["CCL3","CCL4","CCL3L1","TREM2","NPC2","SGK1","AC243829.4","IER3","CCL4L2","A2M","APOE","FCGR2A","AC243829.1","GPX1","FOS","GADD45B","CD14","C3AR1","MS4A6A","MFSD1","RGS1","RHOB","PLAU","MAFB","MARCKS","ATF3","ZFP36L1","EGR2","NR4A2","CTSB","GPNMB","LILRB4","HLA-DMB","TGFBI","KCTD12","IER2","HLA-DMA","C1QC","C3","LGMN","PSAP","CSF1R","MGAT4A","SDS","LPAR6","LHFPL2","ITM2B","RNF130","RNASET2","DUSP1",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_Marco_SPP1",
                  save="Myeloid_Marco_SPP1.markergene.0.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_Marco_SPP1.markergene.0.xxx.pdf")

marker_genes_dict=["FYB1","PRDM1","FPR3","SLCO2B1","ATP6AP2","CREG1","FCGR3A","TAGAP","TIMP2","OTUD1","TYROBP","RASSF4","FCGR2B","EGR1","JUN","ADAP2","PLD3","TMEM51","GPR34","PLXDC2","KLF6","CLEC7A","LAIR1","TMEM176B","ZFP36L2","RB1","FCGR1A","RNASE6","CTSA","HEXA","CTSZ","HLA-DOA","AL590385.2","HERPUD1","CEBPD","LIMS1","IER5","CD74","CXCR4","TLR2","CYFIP1","KLF4","YWHAH","AKR1B1","SCARB2","C1QB","ITGB2","GADD45G","GBP2","GAL3ST4",]
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", var_group_labels="Myeloid_Marco_SPP1",
                  save="Myeloid_Marco_SPP1.markergene.1.rawpdf")
sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="groups_pri", colorbar_title="mean z-score",use_raw=False, vmin=-1, vmax=1, cmap="RdBu_r", save="Myeloid_Marco_SPP1.markergene.1.xxx.pdf")
