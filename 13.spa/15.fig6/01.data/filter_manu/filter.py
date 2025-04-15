import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/01.data/filter_manu")

adata = sc.read_h5ad("pancancer.icar.all.cell.h5ad")
cellist = pd.read_csv("pancancer.icar.all.cell.obs.input")
adata = adata[cellist["cell"],:]
adata.write_h5ad("pancancer.icar.all.cell.filter.h5ad")
adata.obs.to_csv("pancancer.icar.all.cell.filter.obs")

