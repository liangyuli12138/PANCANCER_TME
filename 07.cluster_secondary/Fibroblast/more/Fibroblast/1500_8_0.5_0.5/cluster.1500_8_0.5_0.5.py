import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import harmonypy as hm
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/07.cluster_secondary/Fibroblast/more/Fibroblast/1500_8_0.5_0.5")

adata_concat = sc.read_h5ad("pancancer.ref.0723.final.h5ad")
adata_concat = adata_concat.raw.to_adata()

genelist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/07.cluster_secondary/Fibroblast/more/pancancer.final.0314.var.csv.input.add")
adata_concat = adata_concat[:,genelist["geneID"]]

cellist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/07.cluster_secondary/at_input/Fibroblast.input")
adata_concat = adata_concat[cellist["cell"],:]
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/07.cluster_secondary/at_input/Fibroblast.at",index_col=0)
adata_concat.obs = adata_concat.obs.join(atlist)

sc.pp.filter_cells(adata_concat, min_genes=100)
sc.pp.filter_genes(adata_concat, min_cells=30)
adata_concat.var['mt'] = adata_concat.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata_concat, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata_concat = adata_concat[adata_concat.obs.pct_counts_mt < 20, :]


#sc.pp.normalize_per_cell(adata_concat, counts_per_cell_after=1e4)
sc.pp.normalize_total(adata_concat, target_sum=1e4)
sc.pp.log1p(adata_concat)

#sc.pp.highly_variable_genes(adata_concat, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=3000)
sc.pp.highly_variable_genes(adata_concat, n_top_genes=1500)

adata_concat.raw = adata_concat
adata_concat = adata_concat[:, adata_concat.var.highly_variable]

#adata_concat = adata_concat[:, adata_concat.var.highly_variable]
sc.pp.regress_out(adata_concat, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata_concat, max_value=10)

sc.pp.pca(adata_concat, n_comps=8, svd_solver='arpack', use_highly_variable=True)

meta_data = adata_concat.obs
data_mat = adata_concat.obsm["X_pca"]

ho = hm.run_harmony(data_mat, meta_data, "batch", max_iter_harmony=30)
adata_concat.obsm["X_harmony"] = ho.Z_corr.T

sc.pp.neighbors(adata_concat, use_rep = "X_harmony")

sc.tl.leiden(adata_concat, key_added='groups_secori',resolution=0.5)

sc.tl.umap(adata_concat, min_dist=0.5)

adata_concat.write_h5ad("pancancer.ref.0807.final.Fibroblast.umap.h5ad")
adata_concat.obs.to_csv("pancancer.ref.0807.final.Fibroblast.umap.obs.txt")

sc.pl.umap(adata_concat, color='groups_secori', save= ".Fibroblast.1500_8_0.5_0.5.cluster.png")
sc.pl.umap(adata_concat, color='groups_secori',  legend_loc='on data', save= ".Fibroblast.1500_8_0.5_0.5.on.cluster.png")
sc.pl.umap(adata_concat, color="Tissue", save=".Fibroblast.1500_8_0.5_0.5.Tissue.cluster.png")
sc.pl.umap(adata_concat, color="Phenotype", save=".Fibroblast.1500_8_0.5_0.5.Phenotype.cluster.png")

sc.pl.umap(adata_concat, color=[
'CCL8','GJA4','MYH11','MCAM','RGS5','IL6',
'COL5A1','COL5A2','COL6A3','DCN','FN1','LUM','POSTN','VCAN',
'CD74','HLA-DRA','HLA-DRB1','CCL21','CXCL12',
'KRT19','KRT8','SAA1','SLPI',
'APOA2','FABP1','FABP4','FRZB','GPX3',
'CXCL1','C3','C7','FBLN1','IGFBP6',
],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".Fibroblast.1500_8_0.5_0.5.Fibroblast.umap.markergene.png")

sc.pl.umap(adata_concat, color=[
"POSTN","KRT19","IL6","FBLN2","CD74","KCNN3","SCN7A","P2RY1",
],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".Fibroblast.1500_8_0.5_0.5.less.Fibroblast.umap.markergene.png")


