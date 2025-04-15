import os
import datetime
import json
import argparse
import pandas as pd

import anndata
import scanpy as sc
import harmonypy as hm

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data")

adata_concat = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/D01872C4/D01872C4_cellbin.final.h5ad")
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872C4_cellbin.final.celltype.obs",index_col=0)
adata_concat.obs = atlist

adata_concat.write_h5ad("D01872C4_cellbin.final.celltype.h5ad")
adata_concat.obs.to_csv("D01872C4_cellbin.final.celltype.obs.csv")
adata_concat.var.to_csv("D01872C4_cellbin.final.celltype.var.csv")

