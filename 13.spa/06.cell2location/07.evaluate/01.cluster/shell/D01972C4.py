import os
import datetime
import json
import argparse
import pandas as pd

import anndata
import scanpy as sc
import harmonypy as hm

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/07.evaluate/01.cluster/cluster/")

adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/D01972C4/D01972C4_cellbin.final.h5ad")
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972C4_cellbin.final.celltype.obs",index_col=0)
adata_concat.obs = atlist

adata_concat.var['mt'] = adata_concat.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_concat, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata_concat = adata_concat[adata_concat.obs.pct_counts_mt < 15, :]

sc.pp.normalize_total(adata_concat, target_sum=1e4)
sc.pp.log1p(adata_concat)

sc.pp.highly_variable_genes(adata_concat, flavor='cell_ranger', n_top_genes=2000)
sc.pl.highly_variable_genes(adata_concat, save=".D01972C4.highly_variable_genes.pdf")

adata_concat.raw = adata_concat
adata_concat = adata_concat[:, adata_concat.var.highly_variable]

sc.pp.scale(adata_concat, max_value=10)

sc.pp.pca(adata_concat, n_comps=50, svd_solver='arpack', use_highly_variable=True)
sc.pl.pca(adata_concat, save= '.D01972C4.adata.pca.png')
sc.pl.pca_variance_ratio(adata_concat, log=True, save= ".D01972C4.adata.pcavar_ratio.png")

sc.pp.neighbors(adata_concat, use_rep = "X_pca")
adata_concat.write_h5ad("D01972C4_cellbin.final.neighbors.h5ad")
adata_concat.obs.to_csv("D01972C4_cellbin.final.neighbors.obs.csv")

sc.tl.leiden(adata_concat, key_added='groups',resolution=1)

sc.tl.umap(adata_concat, min_dist=0.4)

adata_concat.write_h5ad("D01972C4_cellbin.final.leiden.h5ad")
adata_concat.obs.to_csv("D01972C4_cellbin.final.leiden.obs.csv")
adata_concat.var.to_csv("D01972C4_cellbin.final.leiden.var.csv")
sc.pl.umap(adata_concat, color='groups', save= ".D01972C4.cluster.png")
sc.pl.umap(adata_concat, color='groups',  legend_loc='on data', save= ".D01972C4.on.cluster.png")

