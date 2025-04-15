import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt


# 读取h5ad文件
adata = sc.read('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/07.gene/merge_25chip_immune_area.nor.h5ad')

# 根据cluster列对细胞进行聚合并求和
cluster_counts = adata.obs.groupby('cluster').sum()

# 提取基因名称
genes = adata.var_names

# 创建一个空的数据框来存储结果
cluster_expression = pd.DataFrame(index=cluster_counts.index, columns=genes)

# 遍历每个cluster，计算每个基因的表达量总和
for cluster in cluster_counts.index:
    cells = adata[adata.obs['cluster'] == cluster]
    gene_expression = cells.X.sum(axis=0)
    cluster_expression.loc[cluster] = gene_expression


# 创建新的h5ad文件
new_adata = sc.AnnData(X=cluster_expression.values, obs=pd.DataFrame(index=cluster_expression.index), var=pd.DataFrame(index=cluster_expression.columns))
sc.pp.normalize_total(new_adata)
sc.pp.log1p(new_adata)
sc.pp.scale(new_adata)
sc.pp.pca(new_adata, n_comps=50)
adata=new_adata

sc.pp.neighbors(adata,n_neighbors = 15, use_rep = "X_pca")
sc.tl.leiden(adata, resolution=0.5)
adata.obs["region_cluster"] = adata.obs["leiden"].astype("category")
sc.tl.umap(adata, min_dist = 0.1)
sc.pl.umap(adata, color=['region_cluster'], size=30, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)
plt.savefig("immune.cluster.r0.5.png",dpi=300, bbox_inches='tight')
adata.obs.to_csv("immune.cluster.r0.5.obs")
adata.write_h5ad("immune.cluster.r0.5.h5ad")

sc.pp.highly_variable_genes(adata, n_top_genes=10000)
variable_genes = adata.var.highly_variable
df = pd.DataFrame(adata.X[:, variable_genes], columns=adata.var_names[variable_genes], index=adata.obs_names)
df_transposed = df.transpose()
df_transposed.to_csv('immune.cluster.gene.X.transposed.xls')

gene_list = pd.read_csv('pick.gene.list', header=None, squeeze=True)
gene_indices = [adata.var_names.get_loc(gene_id) for gene_id in gene_list if gene_id in adata.var_names]
gene_expression = adata.X[:, gene_indices]
df = pd.DataFrame(gene_expression, columns=gene_list, index=adata.obs_names)
df.to_csv('gene_expression.pick.xls')

