import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/12.data_new")

adata = sc.read_h5ad("pancancer.icar.background.cell.h5ad")

cellist = pd.read_csv("pancancer.icar.all.cell.filter.obs.filter.input")
atlist = pd.read_csv("pancancer.icar.all.cell.filter.obs.filter.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = atlist

adata.write_h5ad("pancancer.icar.all.cell.filter.0527.h5ad")
adata.obs.to_csv("pancancer.icar.all.cell.filter.0527.obs")

