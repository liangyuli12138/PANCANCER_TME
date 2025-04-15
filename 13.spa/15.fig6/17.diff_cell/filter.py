import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/17.diff_cell")

adata = sc.read_h5ad("pancancer.icar.all.cell.filter.h5ad")

atlist = pd.read_csv("filter.all.cluster.list.at2",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.write_h5ad("pancancer.icar.all.cell.mergetype2.h5ad")
adata.obs.to_csv("pancancer.icar.all.cell.mergetype2.obs")

