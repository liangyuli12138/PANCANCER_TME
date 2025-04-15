import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/03.cluster_primary/anndata_0.8")
adata = sc.read_h5ad("pancancer.ref.umap.0723.h5ad")

cellist = pd.read_csv("pancancer.ref.umap.0723.obs.txt.stat.rela.input")
adata = adata[cellist["cell"],:]

atlist = pd.read_csv("pancancer.ref.umap.0723.obs.txt.stat.rela.at",index_col=0)
adata.obs = adata.obs.join(atlist)

adata.write_h5ad("pancancer.ref.0723.final.h5ad")
adata.obs.to_csv("pancancer.ref.0723.final.obs.csv")
adata.var.to_csv("pancancer.ref.0723.final.var.csv")

