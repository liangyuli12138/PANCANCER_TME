import os
import datetime
import json
import argparse
import pandas as pd
import anndata
import scanpy as sc

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/07.gene_20240701")
adata = sc.read_h5ad("merge_25chip_immune_area.h5ad")
cellist = pd.read_csv("merge_25chip_immune_area.obs.new.input")
atlist = pd.read_csv("merge_25chip_immune_area.obs.new.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)

sc.pp.filter_cells(adata, min_genes=100)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt < 20, :]
sc.pp.calculate_qc_metrics(adata)
sc.pp.normalize_total(adata)
adata.raw = adata.copy()

adata.write_h5ad("merge_25chip_immune_area.nor.h5ad")
adata.obs.to_csv("merge_25chip_immune_area.nor.obs")
adata.var.to_csv("merge_25chip_immune_area.nor.var")

cluster_counts = adata.obs.groupby('cluster').sum()
# 提取基因名称
genes = adata.var_names

# 创建一个空的数据框来存储结果
cluster_expression = pd.DataFrame(index=cluster_counts.index, columns=genes)

# 遍历每个cluster，计算每个基因的表达量总和
for cluster in cluster_counts.index:
    cells = adata[adata.obs['cluster'] == cluster]
    gene_expression = cells.X.mean(axis=0)
    cluster_expression.loc[cluster] = gene_expression

genes_of_interest = ["CXCL13","CXCR5","FCER2","MS4A1","BCL6","CCL19","CCL21","CCR7","SELL","LAMP3","CXCR4","CD86","CD3D"]
gene_expression_subset = cluster_expression[genes_of_interest]
gene_expression_subset.to_csv("group.gene.normalize.mean.csv")


# 创建新的h5ad文件
#new_adata = sc.AnnData(X=cluster_expression.values, obs=pd.DataFrame(index=cluster_expression.index), var=pd.DataFrame(index=cluster_expression.columns))

# 对gene表达量进行zscore处理

