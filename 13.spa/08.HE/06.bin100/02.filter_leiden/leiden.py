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

df = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/06.bin100/02.filter_leiden/bin100.merge.sub.percent.csv.filter", header=0, index_col=0)

adata = AnnData(df)
sc.pp.neighbors(adata,n_neighbors = 30)
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=0.1)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.5)
sc.pl.umap(adata, color=['region_cluster'], size=4, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("bin200.r0.1.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("bin200.r0.1.obs")

sc.tl.leiden(adata, resolution=0.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.5)
sc.pl.umap(adata, color=['region_cluster'], size=4, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("bin200.r0.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("bin200.r0.2.obs")

sc.tl.leiden(adata, resolution=0.3)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.5)
sc.pl.umap(adata, color=['region_cluster'], size=4, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("bin200.r0.3.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("bin200.r0.3.obs")

sc.tl.leiden(adata, resolution=0.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.5)
sc.pl.umap(adata, color=['region_cluster'], size=4, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("bin200.r0.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("bin200.r0.5.obs")

sc.tl.leiden(adata, resolution=0.8)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.5)
sc.pl.umap(adata, color=['region_cluster'], size=4, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("bin200.r0.8.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("bin200.r0.8.obs")

sc.tl.leiden(adata, resolution=1)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.5)
sc.pl.umap(adata, color=['region_cluster'], size=4, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("bin200.r1.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("bin200.r1.obs")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.5)
sc.pl.umap(adata, color=['region_cluster'], size=4, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("bin200.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("bin200.r1.5.obs")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.5)
sc.pl.umap(adata, color=['region_cluster'], size=4, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("bin200.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("bin200.r2.obs")

sc.tl.leiden(adata, resolution=3)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.8)
sc.pl.umap(adata, color=['region_cluster'], size=4, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("bin200.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("bin200.r2.obs")

