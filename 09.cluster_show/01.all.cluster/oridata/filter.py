import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/01.all.cluster/oridata")
adata = sc.read_h5ad("pancancer.concatenate.oridata.h5ad")

genelist = pd.read_csv("pancancer.final.0314.var.csv.input.add")
adata = adata[:,genelist["geneID"]]

cellist = pd.read_csv("merge.all.all.level.sev.input")
adata = adata[cellist["cell"],:]

atlist = pd.read_csv("merge.all.all.level.sev.at",index_col=0)
adata.obs = adata.obs.join(atlist)

adata.write_h5ad("pancancer.ref.raw.0905.h5ad")
adata.obs.to_csv("pancancer.ref.raw.0905.obs.csv")
adata.var.to_csv("pancancer.ref.raw.0905.var.csv")
