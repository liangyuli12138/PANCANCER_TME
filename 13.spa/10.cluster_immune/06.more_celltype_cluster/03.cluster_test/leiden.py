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

os.chdir=("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/02.cluster")
df = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/02.cluster/cluster.immune.cell.csv", header=0, index_col=0)

adata = AnnData(df)
sc.pp.pca(adata,n_comps=9)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.r1.2.obs")
adata.write_h5ad("immune.cluster.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.r1.5.obs")
adata.write_h5ad("immune.cluster.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.r2.obs")
adata.write_h5ad("immune.cluster.r2.h5ad")

sc.tl.leiden(adata, resolution=3)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.r3.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.r3.obs")
adata.write_h5ad("immune.cluster.r3.h5ad")





