import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import os

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/07.fig5.stage/04.gene_CD4_CD8/00.check_gene")

adata = sc.read_h5ad("merge_25chip_immune_area.h5ad")
csv_data = pd.read_csv('all_cell_meta.csv.at', index_col=0)
filtered_adata = adata[adata.obs.index.isin(csv_data.index)]
adata=filtered_adata
adata.obs = adata.obs.merge(csv_data, left_index=True, right_index=True)
#adata.obs.rename(columns={'clustr': 'cluster'}, inplace=True)

adata.obs['celltype'].unique().tolist()

adata_subset = adata[adata.obs['celltype'] == 'Lymphoid_CD8_Tex', :]

pdcd1_expression = adata_subset[:, 'PDCD1'].X

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

pdcd1_expression_dense = pdcd1_expression.A.flatten() 
sns.histplot(pdcd1_expression_dense, kde=True)
plt.xlabel('PDCD1 Expression')
plt.ylabel('Count')
plt.title('Distribution of PDCD1 Expression in Lymphoid_CD8_Tex Cells')
plt.savefig('check.pdcd1.png')

np.sum(pdcd1_expression_dense)

#函数
def calculate_expression_distribution(adata, obslist, celltype, gene_name):
    # 提取指定细胞类型的细胞
    adata_subset = adata[adata.obs[obslist] == celltype, :]
    # 提取指定基因的表达量
    gene_expression = adata_subset[:, gene_name].X
    gene_expression_dense = gene_expression.A.flatten()  # 将稀疏矩阵转换为稠密矩阵
    # 计算基因表达量的总和
    total_expression = np.sum(gene_expression_dense)
    # 计算对应cell的总数
    total_cells = len(gene_expression_dense)
    # 计算大于0的cell的数目
    positive_cells = np.count_nonzero(gene_expression_dense)
    # 计算总的平均表达量
    mean_expression_total = np.mean(gene_expression_dense)
    # 计算大于0的cell的平均表达量
    mean_expression_positive = np.mean(gene_expression_dense[gene_expression_dense > 0])
    # 绘制基因表达量的数量分布图
    fig, ax = plt.subplots()
    sns.histplot(gene_expression_dense, kde=True)
    plt.xlabel(gene_name + ' Expression')
    plt.ylabel('Count')
    plt.title('Distribution of ' + gene_name + ' Expression in ' + celltype + ' Cells')
    plt.yscale("log")  # 应用对数坐标轴
        # 添加数值标签
    for i, count in enumerate(np.bincount(gene_expression_dense.astype(int))):
        if count > 0:
            ax.text(i, count, str(count), ha='center', va='bottom')
    plt.savefig(obslist + celltype + gene_name + ".png")
    # 返回基因表达量的总和
    header = "total_expression, total_cells, positive_cells, mean_expression_total, mean_expression_positive\n"
    return header, total_expression, total_cells, positive_cells, mean_expression_total, mean_expression_positive
    plt.close()

calculate_expression_distribution(adata,'celltype','Lymphoid_CD4_Treg','FOXP3')
calculate_expression_distribution(adata,'new_groups', 'Lymphoid3','PDCD1')


