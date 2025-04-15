import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import harmonypy as hm
import pandas as pd


os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/04.data.filter_round3/CD8")

adata_concat = sc.read_h5ad("pancancer.ref.0723.final.T_NK.h5ad")
#data_concat = adata_concat.raw.to_adata()

cellist = pd.read_csv("CD8.input")
adata_concat = adata_concat[cellist["cell"],:]

adata_concat.write_h5ad("pancancer.ref.0804.final.CD8.umap.h5ad")
adata_concat.obs.to_csv("pancancer.ref.0804.final.CD8.umap.obs.txt")
adata_concat.var.to_csv("pancancer.ref.0804.final.CD8.umap.var.txt")

