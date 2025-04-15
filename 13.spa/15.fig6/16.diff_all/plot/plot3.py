import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib as mpl

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/16.diff_all/plot")

adata_concat = sc.read_h5ad("pancancer.icar.all.cell.filter.h5ad")

marker_genes_dict = {
"Lymphoid1":["CXCL13","AICDA","IL21R","IL21","CXCR5","CD40","FCGR2B","FCER2","CR1","TNFRSF13C","CD79A","CD79B","POU2AF1","MS4A1","HLA-DRA","CCL21",],
"Lymphoid2":["BCL2","BCL6","TXNIP","FXYD2","A2M","C7","CD40LG","IL7","IL7R",],
"Lymphoid3":["CXCL14","CCR7","CCL19","C3","GZMA","GZMB","GZMK","TIGIT","CTLA4","FOXP3","TNFSF13B","STAT1","ISG15","CCL5","S100A8","S100A9","MDK",],
}

adata_concat = adata_concat[adata_concat.obs['new_groups'].isin(['Lymphoid1', 'Lymphoid2','Lymphoid3']), :]
sc.pp.scale(adata_concat, max_value=10)

sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="new_groups",dot_max=0.8,use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".bluk.lym123.diff.zscor.dotplot.pdf")

sc.pl.matrixplot(adata_concat, marker_genes_dict, groupby="new_groups",dendrogram=False, colorbar_title='mean z-score', use_raw=False, vmin=-0.1, vmax=0.1, cmap='RdBu_r',save=".bluk.lym123.diff.zscor.matrixplot.0620.pdf")

sc.pl.violin(adata_concat, ["FDCSP","FCER2","C7","CXCL13","CXCL9","CXCL5","CXCL1","CXCL2","CXCL3","CXCL8","CXCL11","CXCL17","CXCL16","CXCL10","CXCL12","CXCL14","S100A8","CD2","S100A9","C3","IL10","HLA-A","HLA-B","HLA-C"], groupby="new_groups",stripplot=False)

import matplotlib.pyplot as plt
plt.savefig("violin_plots.pdf")



