import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/02.DoubletFind/stata_check/ec_check_map")
adata = sc.read_h5ad("sub_cluster_EC.2000.30.h5ad")

cellist = pd.read_csv("doublet.input")
atlist = pd.read_csv("doublet.at",index_col=0)

adata = adata[cellist["cell"],:]

adata.obs = adata.obs.join(atlist)

sc.pl.umap(adata, color='Doublets', sort_order = True,
           legend_fontweight="light",
           legend_fontsize=8,size=1, legend_fontoutline=0, save= "ec.doublet..pancancer.groups.umap.png")

