import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/05.diff/B")

adata = sc.read_h5ad("pancancer.icar.all.cell.h5ad")

cellist = pd.read_csv("filter.B.input")
badata =  adata[cellist["cell"],:]
badata.write_h5ad("pancancer.icar.B.cell.h5ad")
badata.obs.to_csv("pancancer.icar.B.cell.obs")

cellist = pd.read_csv("filter.T.input")
tadata =  adata[cellist["cell"],:]
tadata.write_h5ad("pancancer.icar.T.cell.h5ad")
tadata.obs.to_csv("pancancer.icar.T.cell.obs")

cellist = pd.read_csv("filter.cDC.input")
cadata =  adata[cellist["cell"],:]
cadata.write_h5ad("pancancer.icar.cDC.cell.h5ad")
cadata.obs.to_csv("pancancer.icar.cDC.cell.obs")


