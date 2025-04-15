import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import harmonypy as hm
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/05.cluster/CD8_pca/aaaa/bbbb_cccc_dddd")

adata_concat = sc.read_h5ad("pancancer.ref.0804.final.CD8.umap.h5ad")

sc.pp.filter_cells(adata_concat, min_genes=100)
sc.pp.filter_genes(adata_concat, min_cells=30)
adata_concat.var['mt'] = adata_concat.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_concat, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata_concat = adata_concat[adata_concat.obs.pct_counts_mt < 20, :]

#sc.pp.normalize_per_cell(adata_concat, counts_per_cell_after=1e4)
sc.pp.normalize_total(adata_concat, target_sum=1e4)
sc.pp.log1p(adata_concat)

#sc.pp.highly_variable_genes(adata_concat, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=3000)
sc.pp.highly_variable_genes(adata_concat, n_top_genes=bbbb)

adata_concat.raw = adata_concat
adata_concat = adata_concat[:, adata_concat.var.highly_variable]

#adata_concat = adata_concat[:, adata_concat.var.highly_variable]
sc.pp.regress_out(adata_concat, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata_concat, max_value=10)

sc.pp.pca(adata_concat, n_comps=cccc, svd_solver='arpack', use_highly_variable=True)

meta_data = adata_concat.obs
data_mat = adata_concat.obsm["X_pca"]

ho = hm.run_harmony(data_mat, meta_data, "batch")
adata_concat.obsm["X_harmony"] = ho.Z_corr.T

sc.pp.neighbors(adata_concat, use_rep = "X_harmony")

sc.tl.leiden(adata_concat, key_added='groups_secori',resolution=dddd)

sc.tl.umap(adata_concat, min_dist=0.3)

adata_concat.write_h5ad("pancancer.ref.0807.final.aaaa.umap.h5ad")
adata_concat.obs.to_csv("pancancer.ref.0807.final.aaaa.umap.obs.txt")
adata_concat.var.to_csv("pancancer.ref.0807.final.aaaa.umap.var.txt")

sc.pl.umap(adata_concat, color='groups_secori', save= ".cluster.png")
sc.pl.umap(adata_concat, color='groups_secori',  legend_loc='on data', save= ".on.cluster.png")
sc.pl.umap(adata_concat, color="Tissue", save=".Tissue.cluster.png")
sc.pl.umap(adata_concat, color="Phenotype", save=".Phenotype.cluster.png")
sc.pl.umap(adata_concat, color='groups_secori',  legend_loc='on data', save= ".on.cluster.pdf")
sc.pl.umap(adata_concat, color="Tissue", save=".Tissue.cluster.pdf")
sc.pl.umap(adata_concat, color="Phenotype", save=".Phenotype.cluster.pdf")

sc.pl.umap(adata_concat, color=[
"CD8A","TCF7","SELL","CCR7","PRDM1","KLRB1","ITGA1","SEMA4A","CD244","KLRC4","LGR6","CD27","LIMD2","CNN2","EOMES","CCR4","DKK3","MX1","IFIT1","HSPA1A",
"BAG3","NR4A1","GNLY","GZMH","FGFBP2","CTLA4","LAG3","PDCD1","FASLG","FAS","CD69","CD44","PRF1","GZMB","GZMK"
                               ],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".CD8.umap.markergene.png")
sc.pl.umap(adata_concat, color=[
"CD8A","TCF7","SELL","CCR7","PRDM1","KLRB1","ITGA1","SEMA4A","CD244","KLRC4","LGR6","CD27","LIMD2","CNN2","EOMES","CCR4","DKK3","MX1","IFIT1","HSPA1A","BAG3","NR4A1","GNLY","GZMH","FGFBP2","CTLA4","LAG3","PDCD1","FASLG","FAS","CD69","CD44","PRF1","GZMB","GZMK"
                               ],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".CD8.umap.markergene.pdf")

