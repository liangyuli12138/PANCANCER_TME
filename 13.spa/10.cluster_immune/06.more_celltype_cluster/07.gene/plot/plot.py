import pandas as pd
import numpy as np
import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/07.gene/plot")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/07.gene/merge_25chip_immune_area.nor.h5ad")

sc.pp.scale(adata, max_value=10)


# 从文件中读取数据
mf="/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/07.gene/plot/merge.csv.filter.filter"
marker_genes_df = pd.read_csv(mf, delimiter=',', header=None)


# 创建一个空字典
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

sc.pl.matrixplot(adata, marker_genes_dict, groupby="LM",use_raw=False, vmin=-1, vmax=1, cmap='RdBu_r',save=".LM.test.pdf")
sc.pl.matrixplot(adata, marker_genes_dict, groupby="LM",use_raw=False, vmin=-0.5, vmax=0.5, cmap='RdBu_r',save=".LM.test.png")
sc.pl.matrixplot(adata, marker_genes_dict, groupby="LM",use_raw=False, vmin=-1, vmax=1, cmap='RdBu_r',save=".L.gene.pdf")
sc.pl.matrixplot(adata, marker_genes_dict, groupby="LM",use_raw=False, vmin=-1, vmax=1, cmap='RdBu_r',save=".M.gene.pdf")
sc.pl.heatmap(adata, marker_genes_dict, groupby="LM",use_raw=False, vmin=-1, vmax=1, cmap='RdBu_r',figsize=(11,2),save=".LM.test.hm.pdf")

# 调整绘图尺寸和布局
fig, ax = plt.subplots(figsize=(10, 3))
plt.subplots_adjust(left=0.2, bottom=0.3)

# 使用 sc.pl.matrixplot 函数绘制矩阵图
sc.pl.matrixplot(adata, marker_genes_dict, groupby="LM", use_raw=False, vmin=-0.8, vmax=0.8, cmap='RdBu_r', ax=ax)

# 调整每个基因行的宽度和旋转坐标轴标签
ax.set_xticks(range(len(adata.var_names)))
ax.set_xticklabels(adata.var_names, fontsize=6, rotation='vertical')

plt.savefig("gene.LM.diff.pdf")
