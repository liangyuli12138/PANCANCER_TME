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

os.chdir=("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/11.icm/02.border/03.dist/cluster/leiden")
df = pd.read_csv("all.immune.dist.stat.csv", header=0, index_col=0)

adata = AnnData(df)
sc.pp.neighbors(adata,n_neighbors = 10)
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=0.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.r0.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.r0.2.obs")
adata.write_h5ad("immune.cluster.r0.2.h5ad")

sc.tl.leiden(adata, resolution=0.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.r0.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.r0.5.obs")
adata.write_h5ad("immune.cluster.r0.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=300, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.r2.obs")
adata.write_h5ad("immune.cluster.r2.h5ad")

sc.tl.leiden(adata, resolution=1)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.01)
sc.pl.umap(adata, color=['region_cluster'], size=300, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.r1.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.r1.obs")
adata.write_h5ad("immune.cluster.r1.h5ad")





