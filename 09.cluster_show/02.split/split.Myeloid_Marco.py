import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import harmonypy as hm
import pandas as pd


os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split")

adata = sc.read_h5ad("pancancer.ref.0723.raw.h5ad")


cellist = pd.read_csv("split.Myeloid_Marco.input")
adata = adata[cellist["cell"],:]

atlist = pd.read_csv("split.Myeloid_Marco.at",index_col=0)
adata.obs = adata.obs.join(atlist)

adata.raw = adata

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)

adata.write_h5ad("pancancer.split.Myeloid_Marco.0905.h5ad")
adata.obs.to_csv("pancancer.split.Myeloid_Marco.0905.obs.csv")
adata.var.to_csv("pancancer.split.Myeloid_Marco.0905.var.csv")

