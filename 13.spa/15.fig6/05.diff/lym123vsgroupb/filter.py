import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/05.diff/lym123vsgroupb")

adata = sc.read_h5ad("pancancer.icar.all.cell.h5ad")

atlist = pd.read_csv("pancancer.icar.all.cell.obs.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.write_h5ad("pancancer.icar.all.cell.mergetype.h5ad")
adata.obs.to_csv("pancancer.icar.all.cell.mergetype.obs")

