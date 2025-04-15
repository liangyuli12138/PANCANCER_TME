import os
import datetime
import json
import argparse
import pandas as pd
import anndata
import scanpy as sc

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/07.gene")
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
sc.pp.log1p(adata)
adata.raw = adata.copy()


adata.write_h5ad("merge_25chip_immune_area.nor.h5ad")
adata.obs.to_csv("merge_25chip_immune_area.nor.obs")
adata.var.to_csv("merge_25chip_immune_area.nor.var")

