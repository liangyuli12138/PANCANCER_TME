import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/02.zoom/marker")

adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972D1_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("D01972D1_region.csv.at",index_col=0)
adata.obs = adata.obs.join(atlist)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)

marker_genes_dict = {
'immune': ["CD52","CD79A","CCR7","CCL19","CXCL13","IGHG1","IGHG2","IGHG3","IGHG4","IGHGP","IGKC","IGLC2","IGLC3","JCHAIN","MS4A1","MZB1","PIM2","PTGDS","SH2D1A","TRBC2"]
}

sc.pl.dotplot(adata, marker_genes_dict, groupby="region", dot_max=0.9, vmax=3, save=".cluster.markergene.pdf")
sc.pl.dotplot(adata, marker_genes_dict, groupby="region", dot_max=0.9, use_raw=False, colorbar_title="mean z-score", vmin=-2, vmax=2, cmap="RdBu_r",                  save=".cluster.markergene.zscore.pdf")
sc.pl.matrixplot(adata, marker_genes_dict, groupby="region",dendrogram=False, colorbar_title='mean z-score', use_raw=False, vmin=-2, vmax=2, cmap='RdBu_r',           save=".cluster.markergene.zscore.matrix2.pdf")
sc.pl.matrixplot(adata, marker_genes_dict, groupby="region",dendrogram=False, colorbar_title='mean z-score', use_raw=False, vmin=-1, vmax=1, cmap='RdBu_r',           save=".cluster.markergene.zscore.matrix1.pdf")
sc.pl.matrixplot(adata, marker_genes_dict, groupby="region",dendrogram=False, colorbar_title='mean z-score', use_raw=False, vmin=-1, vmax=2, cmap='RdBu_r',           save=".cluster.markergene.zscore.matrix3.pdf")

