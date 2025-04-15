import anndata
import scanpy as sc
import pandas as pd
import numpy as np
import os

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/05.diff/all2")
adata = sc.read_h5ad("pancancer.icar.all.cell.mergetype.h5ad")
adata = adata[adata.obs['difftype'] != 'other']
adata = adata[adata.obs['merge_celltype'].str.contains('Lymphoid_B')]

filtered_adata=adata

adata.obs['groups_clusters'] = adata.obs['new_groups'].astype(str) + '_' + adata.obs['new_cluster'].astype(str)

clusters = sorted(adata.obs['groups_clusters'].unique())

# 获取var中的基因id
genes = adata.var_names

# 创建空的结果矩阵
result_matrix = pd.DataFrame(index=genes, columns=clusters)

# 遍历每个cluster
for cluster in clusters:
    # 获取属于当前cluster的obs索引
    obs_index = adata.obs['groups_clusters'] == cluster
    # 获取属于当前cluster的obs对应的表达矩阵
    expression_matrix = adata[obs_index].X
    # 计算每个基因在当前cluster下的平均表达量
    avg_expression = np.mean(expression_matrix, axis=0)
    # 将平均表达量存入结果矩阵
    result_matrix.loc[genes, cluster] = avg_expression

result_matrix = result_matrix.T
# 将结果矩阵保存为csv文件
result_matrix.to_csv('gene.exp.cluster.xls')
