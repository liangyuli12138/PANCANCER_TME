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

os.chdir=("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/03.cluster_test/out")
df = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/02.cluster/cluster.immune.cell.csv", header=0, index_col=0)


#########pca=3##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=3)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.3.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.3.r1.2.obs")
adata.write_h5ad("immune.cluster.3.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.3.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.3.r1.5.obs")
adata.write_h5ad("immune.cluster.3.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.3.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.3.r2.obs")
adata.write_h5ad("immune.cluster.3.r2.h5ad")
#############################

#########pca=4##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=4)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.4.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.4.r1.2.obs")
adata.write_h5ad("immune.cluster.4.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.4.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.4.r1.5.obs")
adata.write_h5ad("immune.cluster.4.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.4.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.4.r2.obs")
adata.write_h5ad("immune.cluster.4.r2.h5ad")
#############################

#########pca=5##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=5)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.5.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.5.r1.2.obs")
adata.write_h5ad("immune.cluster.5.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.5.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.5.r1.5.obs")
adata.write_h5ad("immune.cluster.5.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.5.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.5.r2.obs")
adata.write_h5ad("immune.cluster.5.r2.h5ad")
#############################

#########pca=6##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=6)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.6.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.6.r1.2.obs")
adata.write_h5ad("immune.cluster.6.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.6.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.6.r1.5.obs")
adata.write_h5ad("immune.cluster.6.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.6.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.6.r2.obs")
adata.write_h5ad("immune.cluster.6.r2.h5ad")
#############################

#########pca=7##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=7)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.7.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.7.r1.2.obs")
adata.write_h5ad("immune.cluster.7.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.7.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.7.r1.5.obs")
adata.write_h5ad("immune.cluster.7.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.7.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.7.r2.obs")
adata.write_h5ad("immune.cluster.7.r2.h5ad")
#############################

#########pca=8##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=8)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.8.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.8.r1.2.obs")
adata.write_h5ad("immune.cluster.8.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.8.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.8.r1.5.obs")
adata.write_h5ad("immune.cluster.8.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.8.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.8.r2.obs")
adata.write_h5ad("immune.cluster.8.r2.h5ad")
#############################

#########pca=9##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=9)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.9.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.9.r1.2.obs")
adata.write_h5ad("immune.cluster.9.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.9.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.9.r1.5.obs")
adata.write_h5ad("immune.cluster.9.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.9.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.9.r2.obs")
adata.write_h5ad("immune.cluster.9.r2.h5ad")
#############################

#########pca=10##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=10)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.10.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.10.r1.2.obs")
adata.write_h5ad("immune.cluster.10.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.10.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.10.r1.5.obs")
adata.write_h5ad("immune.cluster.10.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.10.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.10.r2.obs")
adata.write_h5ad("immune.cluster.10.r2.h5ad")
#############################

#########pca=11##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=11)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.11.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.11.r1.2.obs")
adata.write_h5ad("immune.cluster.11.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.11.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.11.r1.5.obs")
adata.write_h5ad("immune.cluster.11.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.11.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.11.r2.obs")
adata.write_h5ad("immune.cluster.11.r2.h5ad")
#############################

#########pca=12##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=12)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.12.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.12.r1.2.obs")
adata.write_h5ad("immune.cluster.12.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.12.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.12.r1.5.obs")
adata.write_h5ad("immune.cluster.12.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.12.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.12.r2.obs")
adata.write_h5ad("immune.cluster.12.r2.h5ad")
#############################

#########pca=13##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=13)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.13.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.13.r1.2.obs")
adata.write_h5ad("immune.cluster.13.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.13.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.13.r1.5.obs")
adata.write_h5ad("immune.cluster.13.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.13.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.13.r2.obs")
adata.write_h5ad("immune.cluster.13.r2.h5ad")
#############################

#########pca=14##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=14)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.14.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.14.r1.2.obs")
adata.write_h5ad("immune.cluster.14.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.14.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.14.r1.5.obs")
adata.write_h5ad("immune.cluster.14.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.14.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.14.r2.obs")
adata.write_h5ad("immune.cluster.14.r2.h5ad")
#############################

#########pca=15##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=15)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.15.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.15.r1.2.obs")
adata.write_h5ad("immune.cluster.15.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.15.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.15.r1.5.obs")
adata.write_h5ad("immune.cluster.15.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.15.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.15.r2.obs")
adata.write_h5ad("immune.cluster.15.r2.h5ad")
#############################

#########pca=16##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=16)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.16.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.16.r1.2.obs")
adata.write_h5ad("immune.cluster.16.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.16.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.16.r1.5.obs")
adata.write_h5ad("immune.cluster.16.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.16.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.16.r2.obs")
adata.write_h5ad("immune.cluster.16.r2.h5ad")
#############################

#########pca=17##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=17)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.17.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.17.r1.2.obs")
adata.write_h5ad("immune.cluster.17.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.17.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.17.r1.5.obs")
adata.write_h5ad("immune.cluster.17.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.17.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.17.r2.obs")
adata.write_h5ad("immune.cluster.17.r2.h5ad")
#############################

#########pca=18##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=18)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.18.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.18.r1.2.obs")
adata.write_h5ad("immune.cluster.18.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.18.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.18.r1.5.obs")
adata.write_h5ad("immune.cluster.18.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.18.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.18.r2.obs")
adata.write_h5ad("immune.cluster.18.r2.h5ad")
#############################

#########pca=19##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=19)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.19.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.19.r1.2.obs")
adata.write_h5ad("immune.cluster.19.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.19.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.19.r1.5.obs")
adata.write_h5ad("immune.cluster.19.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.19.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.19.r2.obs")
adata.write_h5ad("immune.cluster.19.r2.h5ad")
#############################

#########pca=20##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=20)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.20.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.20.r1.2.obs")
adata.write_h5ad("immune.cluster.20.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.20.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.20.r1.5.obs")
adata.write_h5ad("immune.cluster.20.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.20.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.20.r2.obs")
adata.write_h5ad("immune.cluster.20.r2.h5ad")
#############################

