adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/07.gene/merge_25chip_immune_area.nor.h5ad")
sc.pp.scale(adata, max_value=10)

df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
#df.to_csv("zscore.csv")
mean_values = df.groupby(adata.obs['new_groups']).mean()
mean_values.to_csv("zscore.csv")


df = pd.read_csv("zscore.csv", index_col=0)
df_transposed = df.transpose()
df_transposed.to_csv("zscore_transposed.csv")

