import scanpy as sc
import pandas as pd

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

cluster2group = pd.read_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/07.gene/group.zscore/cluster2group.csv', index_col=0)
new_adata.obs=cluster2group

ori_expression = pd.DataFrame(new_adata.X, index=new_adata.obs.index, columns=new_adata.var.index).groupby(new_adata.obs['new_groups']).mean()

df_transposed = ori_expression.transpose()
df_transposed.to_csv("ori.normalize.csv")


# 提取名称包含"Lymphoid"的cluster
#lymphoid_clusters = new_adata.obs[new_adata.obs['new_groups'].str.contains('Lymphoid')].index
lymphoid_clusters = new_adata.obs[new_adata.obs['new_groups'].isin(['Lymphoid1', 'Lymphoid2', 'Lymphoid3'])].index

# 提取对应的表达量数据
#lymphoid_expression = new_adata.X[new_adata.obs.index.isin(lymphoid_clusters)]

lymphoid_adata = new_adata[new_adata.obs.index.isin(lymphoid_clusters)]

# 对gene表达量进行zscore处理
sc.pp.scale(lymphoid_adata)

# 计算new_groups的平均表达量
mean_expression = pd.DataFrame(lymphoid_adata.X, index=lymphoid_adata.obs.index, columns=lymphoid_adata.var.index).groupby(lymphoid_adata.obs['new_groups']).mean()

# 将结果保存为group.zscore.csv文件
mean_expression.to_csv('group.zscore.csv')

df_transposed = mean_expression.transpose()
df_transposed.to_csv("zscore_transposed.csv")
df_transposed.to_csv("zscore_transposed.lym.normalizecsv")


