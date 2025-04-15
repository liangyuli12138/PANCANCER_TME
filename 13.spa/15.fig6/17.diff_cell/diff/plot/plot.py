import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib as mpl

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/17.diff_cell/diff/plot")

adata_concat = sc.read_h5ad("pancancer.icar.all.cell.mergetype.h5ad")

marker_genes_dict = {
"Lymphoid_B_Lymphoid1":["AICDA","CD83","BCL6","CD79A","CXCR5","MS4A1","HLA-DRA","CD40","CXCR4","IL21R","TNFRSF13C","MZB1","IL10"],
"Lymphoid_B_Lymphoid2":["AICDA","CD83","BCL6","CD79A","CXCR5","MS4A1","HLA-DRA","CD40","CXCR4","IL21R","TNFRSF13C","MZB1","IL10"],
"Lymphoid_B_Lymphoid3":["AICDA","CD83","BCL6","CD79A","CXCR5","MS4A1","HLA-DRA","CD40","CXCR4","IL21R","TNFRSF13C","MZB1","IL10"],
}

adata_concat = adata_concat[adata_concat.obs['difftype'].isin(['Lymphoid_B_Lymphoid1', 'Lymphoid_B_Lymphoid2','Lymphoid_B_Lymphoid3']), :]
sc.pp.scale(adata_concat, max_value=10)

sc.pl.dotplot(adata_concat, marker_genes_dict, groupby="new_groups",dot_max=0.8,use_raw=False, colorbar_title="mean z-score", vmin=-1, vmax=1, cmap="RdBu_r",save=".bluk.lym123.diff.zscor.dotplot.pdf")

sc.pl.matrixplot(adata_concat, marker_genes_dict, groupby="new_groups",dendrogram=False, colorbar_title='mean z-score', use_raw=False, vmin=-0.1, vmax=0.1, cmap='RdBu_r',save=".bluk.b.diff.zscor.matrixplot.png")

sc.pl.violin(adata_concat, ["FDCSP","FCER2","C7","CXCL13","CXCL9","CXCL5","CXCL1","CXCL2","CXCL3","CXCL8","CXCL11","CXCL17","CXCL16","CXCL10","CXCL12","CXCL14","S100A8","CD2","S100A9","C3","IL10","HLA-A","HLA-B","HLA-C"], groupby="new_groups",stripplot=False)

import matplotlib.pyplot as plt
plt.savefig("violin_plots.pdf")



