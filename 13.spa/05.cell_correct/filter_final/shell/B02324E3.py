import sys
import scanpy as sc
import pandas as pd
import numpy as np
import os


adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/B02324E3/B02324E3_cellbin.filter.h5ad")
cellist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/B02324E3/B02324E3_cellbin.filter.gene.list.cellbin")
adata = adata[cellist.cellbin.astype(str), :]
adata.obs.index = adata.obs['center_x'].astype(str) + '_' + adata.obs['center_y'].astype(str)
adata.write_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/B02324E3/B02324E3_cellbin.final.h5ad")
adata.obs.to_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/B02324E3/B02324E3_cellbin.final.obs")

