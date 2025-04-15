import os
import scanpy as sc
import pandas as pd

# 读取h5ad文件
adata = sc.read('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/07.gene/merge_25chip_immune_area.nor.h5ad')

cluster_counts = adata.obs.groupby('cluster').sum()

# 提取基因名称
genes = adata.var_names

# 创建一个空的数据框来存储结果
cluster_expression = pd.DataFrame(index=cluster_counts.index, columns=genes)

# 遍历每个cluster，计算每个基因的表达量总和
for cluster in cluster_counts.index:
    cells = adata[adata.obs['cluster'] == cluster]
    gene_expression = cells.X.meanll(axis=0)
    cluster_expression.loc[cluster] = gene_expression

# 创建新的h5ad文件
new_adata = sc.AnnData(X=cluster_expression.values, obs=pd.DataFrame(index=cluster_expression.index), var=pd.DataFrame(index=cluster_expression.columns))

# 对gene表达量进行zscore处理
sc.pp.scale(new_adata)

# 提取new_adata.X矩阵
matrix = new_adata.X

# 构建一个数据框
df = pd.DataFrame(matrix, index=new_adata.obs.index, columns=new_adata.var.index)

df_transposed = df.transpose()
df_transposed.to_csv("ori.groups.zscore_transposed.csv")

