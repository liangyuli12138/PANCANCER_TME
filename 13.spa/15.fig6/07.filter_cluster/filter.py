import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/07.filter_cluster")
adata = sc.read_h5ad("pancancer.icar.all.cell.filter.h5ad")

cellist = pd.read_csv("filter.input")
atlist = pd.read_csv("filter.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)

adata.write_h5ad("pancancer.icar.all.cell.filmerge.h5ad")
adata.obs.to_csv("pancancer.icar.all.cell.filmerge.obs")


