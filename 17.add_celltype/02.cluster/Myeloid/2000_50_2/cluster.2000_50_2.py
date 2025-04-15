import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import harmonypy as hm

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/17.add_celltype/02.cluster/Myeloid/2000_50_2")

adata_concat = sc.read_h5ad("pancancer.ref.0723.final.Myeloid.h5ad")

sc.pp.filter_cells(adata_concat, min_genes=100)
sc.pp.filter_genes(adata_concat, min_cells=30)
adata_concat.var['mt'] = adata_concat.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_concat, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata_concat = adata_concat[adata_concat.obs.pct_counts_mt < 20, :]

sc.pl.violin(adata_concat, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],jitter=False, multi_panel=True,save= ".counts.violin.pdf")
sc.pl.scatter(adata_concat, x='total_counts', y='pct_counts_mt',save= ".counts_mt.scatter.pdf")
sc.pl.scatter(adata_concat, x='total_counts', y='n_genes_by_counts',save= ".counts_genes.scatter.pdf")

#sc.pp.normalize_per_cell(adata_concat, counts_per_cell_after=1e4)
sc.pp.normalize_total(adata_concat, target_sum=1e4)
sc.pp.log1p(adata_concat)

#sc.pp.highly_variable_genes(adata_concat, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=3000)
sc.pp.highly_variable_genes(adata_concat, n_top_genes=2000)
sc.pl.highly_variable_genes(adata_concat, save=".highly_variable_genes.pdf")

adata_concat.raw = adata_concat
adata_concat = adata_concat[:, adata_concat.var.highly_variable]

#adata_concat = adata_concat[:, adata_concat.var.highly_variable]
sc.pp.regress_out(adata_concat, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata_concat, max_value=10)

sc.pp.pca(adata_concat, n_comps=50, svd_solver='arpack', use_highly_variable=True)
sc.pl.pca(adata_concat, save= '.adata.pca.png')
sc.pl.pca_variance_ratio(adata_concat, log=True, save= ".a.adata.pcavar_ratio.png")
sc.pl.pca(adata_concat, save= '.adata.pca.pdf')
sc.pl.pca_variance_ratio(adata_concat, log=True, save= ".a.adata.pcavar_ratio.pdf")

meta_data = adata_concat.obs
data_mat = adata_concat.obsm["X_pca"]

ho = hm.run_harmony(data_mat, meta_data, "batch")
adata_concat.obsm["X_harmony"] = ho.Z_corr.T

sc.pp.neighbors(adata_concat, use_rep = "X_harmony")

sc.tl.leiden(adata_concat, key_added='groups_secori',resolution=2)

sc.tl.umap(adata_concat, min_dist=0.4)

adata_concat.write_h5ad("pancancer.ref.0723.final.Myeloid.umap.h5ad")
adata_concat.obs.to_csv("pancancer.ref.0723.final.Myeloid.umap.obs.txt")
adata_concat.var.to_csv("pancancer.ref.0723.final.Myeloid.umap.var.txt")

sc.pl.umap(adata_concat, color='groups_secori', save= ".cluster.png")
sc.pl.umap(adata_concat, color='groups_secori',  legend_loc='on data', save= ".on.cluster.png")
sc.pl.umap(adata_concat, color="Tissue", save=".Tissue.cluster.png")
sc.pl.umap(adata_concat, color="Phenotype", save=".Phenotype.cluster.png")
sc.pl.umap(adata_concat, color='groups_secori',  legend_loc='on data', save= ".on.cluster.pdf")
sc.pl.umap(adata_concat, color="Tissue", save=".Tissue.cluster.pdf")
sc.pl.umap(adata_concat, color="Phenotype", save=".Phenotype.cluster.pdf")

