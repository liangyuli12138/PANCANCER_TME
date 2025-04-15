import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/19.contour.after/01.data")

adata = sc.read_h5ad("pancancer.icar.background.cell.h5ad")
cellist = pd.read_csv("contour.input")
atlist = pd.read_csv("contour.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata.obs['typediff'] = adata.obs['celltype_merge'].str.cat(adata.obs['diff'], sep='_')

adata.write_h5ad("pancancer.icar.contour.cell.h5ad")
adata.obs.to_csv("pancancer.icar.contour.cell.obs")

epithelium_data = adata[adata.obs['celltype_merge'] == 'Epithelium', :]
epithelium_data.write_h5ad("pancancer.icar.contour.epithelium.cell.h5ad")
epithelium_data.obs.to_csv("pancancer.icar.contour.epithelium.cell.obs")

epithelium_data = adata[adata.obs['celltype_merge'].isin(['Epithelium', 'Fibroblast','Marcophage','NK_cell','B_cell']), :]
epithelium_data.write_h5ad("pancancer.icar.contour.filter.cell.h5ad")
epithelium_data.obs.to_csv("pancancer.icar.contour.filter.cell.obs")

