import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import harmonypy as hm

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/06.K_and_L/01.cluster/cluster/")

adata_concat = sc.read_h5ad("aaaa_cellbin.final.h5ad")

adata_concat.var['mt'] = adata_concat.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_concat, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata_concat = adata_concat[adata_concat.obs.pct_counts_mt < 15, :]

sc.pp.normalize_total(adata_concat, target_sum=1e4)
sc.pp.log1p(adata_concat)

sc.pp.highly_variable_genes(adata_concat, flavor='cell_ranger', n_top_genes=2000)
sc.pl.highly_variable_genes(adata_concat, save=".aaaa.highly_variable_genes.pdf")

adata_concat.raw = adata_concat
adata_concat = adata_concat[:, adata_concat.var.highly_variable]

sc.pp.scale(adata_concat, max_value=10)

sc.pp.pca(adata_concat, n_comps=50, svd_solver='arpack', use_highly_variable=True)
sc.pl.pca(adata_concat, save= '.aaaa.adata.pca.png')
sc.pl.pca_variance_ratio(adata_concat, log=True, save= ".aaaa.adata.pcavar_ratio.png")

sc.pp.neighbors(adata_concat, use_rep = "X_pca")

sc.tl.leiden(adata_concat, key_added='groups',resolution=0.5)

sc.tl.umap(adata_concat, min_dist=0.4)

adata_concat.write_h5ad("aaaa_cellbin.final.leiden.h5ad")
adata_concat.obs.to_csv("aaaa_cellbin.final.leiden.obs.csv")
adata_concat.var.to_csv("aaaa_cellbin.final.leiden.var.csv")
sc.pl.umap(adata_concat, color='groups', save= ".aaaa.cluster.png")
sc.pl.umap(adata_concat, color='groups',  legend_loc='on data', save= ".aaaa.on.cluster.png")

