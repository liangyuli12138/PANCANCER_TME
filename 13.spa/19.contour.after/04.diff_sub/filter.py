import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/19.contour.after/04.diff_sub")

adata = sc.read_h5ad("pancancer.icar.contour.cell.ori.h5ad")

adata.obs['subdiff'] = adata.obs['celltype'].str.cat(adata.obs['diff'], sep='_')

adata.write_h5ad("pancancer.icar.contour.cell.h5ad")
adata.obs.to_csv("pancancer.icar.contour.cell.obs")

