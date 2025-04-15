import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/01.data/anndata_0.8")
adata = sc.read_h5ad("pancancer.concatenate.oridata.h5ad")

genelist = pd.read_csv("pancancer.final.0314.var.csv.input.add")
adata = adata[:,genelist["geneID"]]

cellist = pd.read_csv("pancancer.final.0315.obs.csv.tlike.input")
adata = adata[cellist["cell"],:]

atlist = pd.read_csv("pancancer.final.0315.obs.csv.tlike.at",index_col=0)
adata.obs = adata.obs.join(atlist)

adata.write_h5ad("pancancer.ref.0721.raw.h5ad")
adata.obs.to_csv("pancancer.ref.0721.raw.obs.csv")
adata.var.to_csv("pancancer.ref.0721.raw.var.csv")

sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

adata.write_h5ad("pancancer.ref.0721.nor.h5ad")
adata.obs.to_csv("pancancer.ref.0721.nor.obs.csv")
adata.var.to_csv("pancancer.ref.0721.nor.var.csv")

