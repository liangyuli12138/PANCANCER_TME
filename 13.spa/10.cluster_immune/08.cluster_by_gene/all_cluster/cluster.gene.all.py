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

cluster_expr=cluster_expression

# 计算每个基因在每个细胞簇中的表达量和总表达量
cluster_expr_sum = np.sum(cluster_expr, axis=0)
cluster_expr_var = np.var(cluster_expr, axis=0)
cluster_expr_count = np.sum(cluster_expr > 0, axis=0)

import matplotlib.pyplot as plt
import os
os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/08.cluster_by_gene/all_cluster")

plt.figure(figsize=(10, 5))
plt.hist(cluster_expr_sum, bins=50, edgecolor='black')
plt.xlabel('cluster_expr_sum')
plt.ylabel('Frequency')
plt.title('Distribution of cluster_expr_sum')
plt.savefig('cluster_expr_sum_distribution.png')
plt.close()

# 绘制 cluster_expr_var 的分布曲线并保存为文件
plt.figure(figsize=(10, 5))
plt.hist(cluster_expr_var, bins=50, edgecolor='black')
plt.xlabel('cluster_expr_var')
plt.ylabel('Frequency')
plt.title('Distribution of cluster_expr_var')
plt.savefig('cluster_expr_var_distribution.png')
plt.close()

# 绘制 cluster_expr_count 的分布曲线并保存为文件
plt.figure(figsize=(10, 5))
plt.hist(cluster_expr_count, bins=50, edgecolor='black')
plt.xlabel('cluster_expr_count')
plt.ylabel('Frequency')
plt.title('Distribution of cluster_expr_count')
plt.savefig('cluster_expr_count_distribution.png')
plt.close()

# 过滤条件
gene_filter = (cluster_expr_count > 30) & (cluster_expr_sum > 1000) & (cluster_expr_var > 0.1)

# 获取满足过滤条件的基因索引
filtered_gene_indices = np.where(gene_filter)[0]

# 根据过滤条件过滤聚类表达矩阵
filtered_cluster_expr = cluster_expr.iloc[:, filtered_gene_indices]

# 更新adata对象的聚类表达矩阵
adata.obsm['X_cluster_expression'] = filtered_cluster_expr.values


import numpy as np
from scipy.stats import pearsonr

num_cells = filtered_cluster_expr.shape[0]
num_genes = filtered_cluster_expr.shape[1]

# 创建一个全零矩阵来保存相关性结果
correlation_matrix = np.zeros((num_cells, num_cells))

# 遍历每对细胞
for i in range(num_cells):
    for j in range(num_cells):
        # 获取细胞i和j的基因表达
        gene_expr_i = filtered_cluster_expr.iloc[i]
        gene_expr_j = filtered_cluster_expr.iloc[j]
        # 计算基因表达的相关性
        correlation, _ = pearsonr(gene_expr_i, gene_expr_j)
        # 将相关性保存到矩阵中
        correlation_matrix[i, j] = correlation

# 获取基因和细胞名称列表
gene_names = filtered_cluster_expr.columns.tolist()
cell_names = filtered_cluster_expr.index.tolist()
correlation_df = pd.DataFrame(correlation_matrix, index=cell_names, columns=cell_names)

# 使用 seaborn 创建热图
sns.clustermap(correlation_df, cmap='viridis')

# 调整热图大小
plt.gcf().set_size_inches(10, 10)

# 保存热图为 PNG 文件
plt.savefig('heatmap.png', dpi=300)

correlation_df.to_csv("correlation_df.xls")



