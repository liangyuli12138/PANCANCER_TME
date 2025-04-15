import sys
import scanpy as sc
import pandas as pd
import numpy as np
import os


adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/B02324E5/B02324E5_cellbin.filter.h5ad")
cellist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/B02324E5/B02324E5_cellbin.filter.gene.list.cellbin")
adata = adata[cellist.cellbin.astype(str), :]
adata.write_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/B02324E5/B02324E5_cellbin.final.h5ad")
adata.write_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/B02324E5/B02324E5_cellbin.final.obs")

