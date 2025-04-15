import scanpy as sc
import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.sparse import csr_matrix
import glob
from anndata import AnnData

df = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/05.stat/leiden2/cluster.immune.cell.csv", header=0, index_col=0)

adata = AnnData(df)
sc.pp.neighbors(adata,n_neighbors = 15)
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")


sc.tl.leiden(adata, resolution=0.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.r0.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.r0.5.obs")
adata.write_h5ad("immune.cluster.r0.5.h5ad")

sc.tl.leiden(adata, resolution=0.6)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.r0.6.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.r0.6.obs")
adata.write_h5ad("immune.cluster.r0.6.h5ad")

