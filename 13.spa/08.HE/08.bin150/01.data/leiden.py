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

df = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/08.bin150/01.data/bin150.merge.sub.percent.csv.filter", header=0, index_col=0)

adata = AnnData(df)
sc.pp.neighbors(adata,n_neighbors = 15)
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=0.1)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.5)
sc.pl.umap(adata, color=['region_cluster'], size=4, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("bin150.r0.1.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("bin150.r0.1.obs")
adata.write_h5ad("bin150.r0.1.h5ad")


sc.tl.leiden(adata, resolution=0.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.5)
sc.pl.umap(adata, color=['region_cluster'], size=4, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("bin150.r0.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("bin150.r0.2.obs")
adata.write_h5ad("bin150.r0.2.h5ad")

sc.tl.leiden(adata, resolution=0.3)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.5)
sc.pl.umap(adata, color=['region_cluster'], size=4, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("bin150.r0.3.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("bin150.r0.3.obs")
adata.write_h5ad("bin150.r0.3.h5ad")

sc.tl.leiden(adata, resolution=0.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.5)
sc.pl.umap(adata, color=['region_cluster'], size=4, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("bin150.r0.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("bin150.r0.5.obs")
adata.write_h5ad("bin150.r0.5.h5ad")

sc.tl.leiden(adata, resolution=0.8)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.5)
sc.pl.umap(adata, color=['region_cluster'], size=4, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("bin150.r0.8.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("bin150.r0.8.obs")
adata.write_h5ad("bin150.r0.8.h5ad")

sc.tl.leiden(adata, resolution=1)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.5)
sc.pl.umap(adata, color=['region_cluster'], size=4, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("bin150.r1.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("bin150.r1.obs")
adata.write_h5ad("bin150.r1.h5ad")

