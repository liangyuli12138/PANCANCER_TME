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
sc.pp.log1p(adata)

# 读取csv文件
df = pd.read_csv('cluster.marker.lsit.forplot.input', delimiter=',', header=None)

# 提取基因列表
gene_list = df[1].tolist()

# 根据基因id提取adata中的基因
genes = adata.var_names.intersection(gene_list)

# 计算zscore
adata_sub = adata[:, genes]
sc.pp.scale(adata_sub, zero_center=True, max_value=3)

# 计算zscore的平均值
mean_zscore = adata_sub.to_df().groupby(adata_sub.obs['icm']).mean()

# 输出结果为csv文件
mean_zscore.to_csv('output.csv', index=True)

