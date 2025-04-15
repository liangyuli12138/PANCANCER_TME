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

cluster_groups = pd.read_csv('cluster2group.csv', index_col=0)
new_adata.obs['new_groups'] = cluster_groups

# 筛选包含'Lymphoid'关键词的细胞
lymphoid_cells = new_adata.obs['new_groups'].str.contains('Lymphoid')
new_adata=new_adata[lymphoid_cells]

sc.pp.normalize_total(new_adata)
sc.pp.log1p(new_adata)
sc.pp.scale(new_adata)

gene_list = pd.read_csv('pick.gene.list', header=None, squeeze=True)
gene_list = [gene for gene in gene_list if gene in new_adata.var_names]
gene_indices = [new_adata.var_names.get_loc(gene_id) for gene_id in gene_list if gene_id in new_adata.var_names]
gene_expression = new_adata.X[:, gene_indices]
df = pd.DataFrame(gene_expression, columns=gene_list, index=new_adata.obs_names)
df.to_csv('gene_expression.pick.xls')

filtered_cluster_expr = df

from scipy.stats import pearsonr

num_cells = filtered_cluster_expr.shape[0]
num_genes = filtered_cluster_expr.shape[1]
correlation_matrix = np.zeros((num_cells, num_cells))

for i in range(num_cells):
    for j in range(num_cells):
        # 获取细胞i和j的基因表达
        gene_expr_i = filtered_cluster_expr.iloc[i]
        gene_expr_j = filtered_cluster_expr.iloc[j]
        # 计算基因表达的相关性
        correlation, _ = pearsonr(gene_expr_i, gene_expr_j)
        # 将相关性保存到矩阵中
        correlation_matrix[i, j] = correlation

gene_names = filtered_cluster_expr.columns.tolist()
cell_names = filtered_cluster_expr.index.tolist()
correlation_df = pd.DataFrame(correlation_matrix, index=cell_names, columns=cell_names)

correlation_df.to_csv("correlation.lym.xls")

