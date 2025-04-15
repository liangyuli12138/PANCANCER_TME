#########pca=aaaa##############
adata = AnnData(df)
sc.pp.pca(adata,n_comps=aaaa)
sc.pp.neighbors(adata,n_neighbors = 50, use_rep = "X_pca")
#adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")

sc.tl.leiden(adata, resolution=1.2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.aaaa.r1.2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.aaaa.r1.2.obs")
adata.write_h5ad("immune.cluster.aaaa.r1.2.h5ad")

sc.tl.leiden(adata, resolution=1.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.aaaa.r1.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.aaaa.r1.5.obs")
adata.write_h5ad("immune.cluster.aaaa.r1.5.h5ad")

sc.tl.leiden(adata, resolution=2)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.aaaa.r2.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.aaaa.r2.obs")
adata.write_h5ad("immune.cluster.aaaa.r2.h5ad")
#############################

