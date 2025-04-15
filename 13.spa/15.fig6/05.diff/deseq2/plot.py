import scanpy as sc
import pandas as pd

# 读取pancancer.icar.all.cell.h5ad文件为adata
adata = sc.read_h5ad('pancancer.icar.all.cell.h5ad')

# 更新adata的obs 'new_cluster'列
adata.obs['new_cluster'] = adata.obs['merge_groups'].astype(str) + '_' + adata.obs['new_cluster'].astype(str)
adata.obs['new_cluster'] = adata.obs['new_cluster'].str.replace(r'Cluster\d+', '', regex=True)


# 读取diff.Lymphoid.filter文件为csv文件
diff_genes = pd.read_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/02.diff/deseq2/merge.deseq2.list')
diff_gene_list = diff_genes["gene"].tolist()
existing_genes = [gene for gene in diff_gene_list if gene in adata.var_names]


# 提取h5ad文件中在diff_genes中的基因
adata_genes = adata[:,existing_genes]
sc.pp.normalize_total(adata_genes)

# 创建一个空的DataFrame来存储聚类后的表达量
cluster_expr_df = pd.DataFrame(index=adata.obs['new_cluster'].unique(), columns=adata_genes.var_names)

# 对每个聚类进行遍历，计算表达量加和
for cluster in adata.obs['new_cluster'].unique():
    cluster_expr = adata_genes[adata.obs['new_cluster'] == cluster].X.mean(axis=0)
    cluster_expr_df.loc[cluster] = cluster_expr

# 对cluster_expr_df进行转置
cluster_expr_df = cluster_expr_df.T

# 保留含有"Lym"的列
cluster_expr_df = cluster_expr_df.filter(like='Lym', axis=1)

# 根据表头进行排序
cluster_expr_df = cluster_expr_df.sort_index(axis=1)

# 将cluster_expr_df保存为CSV文件
cluster_expr_df.to_csv('cluster_expr.csv')
