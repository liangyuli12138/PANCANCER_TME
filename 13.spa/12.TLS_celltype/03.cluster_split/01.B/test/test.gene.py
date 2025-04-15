import pandas as pd
import numpy as np

# 假设你的单细胞adata对象名为adata
# 假设你的基因列表已经保存在gene_list变量中

# 创建一个空的DataFrame来存储结果
result_df = pd.DataFrame(index=adata.obs['celltype'].unique(), columns=gene_list)

# 遍历每个细胞类型
for cell_type in adata.obs['celltype'].unique():
    # 获取该细胞类型的细胞索引
    cell_indices = adata.obs['celltype'] == cell_type
    # 计算该细胞类型下每个基因的平均表达量
    gene_expression = np.mean(adata[cell_indices, gene_list].X, axis=0)
    # 将结果存储到DataFrame中
    result_df.loc[cell_type] = gene_expression

# 将结果保存为CSV文件
result_df.to_csv('gene_expression_matrix.csv')

