import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/02.DoubletFind/filter/anndata_0.7")
adata = sc.read_h5ad("pancancer.concatenate.oridata.h5ad")

cellist = pd.read_csv("merge.doublet.csv.input")
adata = adata[cellist["cell"],:]

atlist = pd.read_csv("merge.doublet.csv.at",index_col=0)
adata.obs = adata.obs.join(atlist)

adata.write_h5ad("pancancer.ref.0723.raw.h5ad")
adata.obs.to_csv("pancancer.ref.0723.raw.obs.csv")
adata.var.to_csv("pancancer.ref.0723.raw.var.csv")

sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)

adata.write_h5ad("pancancer.ref.0723.nor.h5ad")
adata.obs.to_csv("pancancer.ref.0723.nor.obs.csv")
adata.var.to_csv("pancancer.ref.0723.nor.var.csv")

