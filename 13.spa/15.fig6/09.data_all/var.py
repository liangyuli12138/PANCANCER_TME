import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/09.data_all")

adata = sc.read_h5ad("pancancer.icar.background.cell.h5ad")

adata.var.to_csv("pancancer.icar.background.cell.var")
