import pandas as pd
import numpy as np
import os
import datetime
import json
import argparse

import anndata
import scanpy as sc

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/10.co_existence/marker_immune/out/immune.diff.gene")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/03.cluster_primary/pancancer.ref.0723.raw.h5ad")
cellist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/10.co_existence/marker_immune/filter.input")
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/10.co_existence/marker_immune/filter.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
sc.pp.normalize_total(adata, target_sum=1e4)
#sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
mf="cluster.marker.lsit.forplot.input"
marker_genes_df = pd.read_csv(mf, delimiter=',', header=None)
gene_ids = []
with open('cluster.marker.lsit.forplot.input', 'r') as file:
    for line in file:
        gene_id = line.strip().split(',')[1]
        gene_ids.append(gene_id)
adata_subset = adata[:, gene_ids]
sc.pp.scale(adata_subset, zero_center=True, max_value=3)

marker_genes_dict = {}
# 遍历数据框的行
for row in marker_genes_df.iterrows():
    key = row[1][0]  # 第一列的值作为键
    value = row[1][1]  # 第二列的值作为值
    # 跳过NaN值
    if pd.isna(value):
        continue
    # 初始化值为一个空列表，如果键不存在于字典中
    if key not in marker_genes_dict:
        marker_genes_dict[key] = []
    # 将值添加到键对应的列表中
    marker_genes_dict[key].append(value)

sc.pl.matrixplot(adata, marker_genes_dict, use_raw=False, groupby="icm",vmin=-1, vmax=1, cmap='RdBu_r',save=".icm.diff.marker.pdf",figsize=(8,2),categories_order=["ICM1", "ICM2", "ICM3"])

