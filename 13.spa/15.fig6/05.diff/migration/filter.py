import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/05.diff/migration")

adata = sc.read_h5ad("pancancer.icar.all.cell.filter.h5ad")

cellist = pd.read_csv("filter.mig.input")
badata =  adata[cellist["cell"],:]
badata.write_h5ad("pancancer.icar.mig.cell.h5ad")
badata.obs.to_csv("pancancer.icar.mig.cell.obs")
