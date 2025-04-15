import scanpy as sc
import pandas as pd
import numpy as np

# 读取h5ad文件
adata = sc.read('pancancer.icar.background.cell.h5ad')

# 指定基因列表
genes_list = [
    'ACTA1',
    'CD14',
    'CD163',
    'CD38',
    'CD47',
    'CD68',
    'CEACAM8',
    'COL1A1',
    'FAP',
    'FUT4',
    'HAVCR2',
    'ITGAM',
    'ITGAX',
    'KLRD1',
    'KRT19',
    'LAG3',
    'MPO',
    'MRC1',
    'NCAM1',
    'PECAM1',
    'TNFRSF4'
]

# 创建结果文件并写入列名
result_file = 'result.csv'
pd.Series(adata.obs['celltype_merge'].unique()).to_frame().T.to_csv(result_file, index=True)

# 计算每个细胞类型中基因的平均表达量并追加到结果文件
for gene in genes_list:
    result = pd.DataFrame(index=[gene], columns=adata.obs['celltype_merge'].unique())
    
    for celltype in adata.obs['celltype_merge'].unique():
        mean_expression = adata[adata.obs['celltype_merge'] == celltype].X[:, adata.var_names == gene].mean()
        result.loc[:, celltype] = mean_expression
    
    result.to_csv(result_file, mode='a', header=False, index=True)
