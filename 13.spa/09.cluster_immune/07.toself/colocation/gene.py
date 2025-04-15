import scanpy as sc
import pandas as pd

adata = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972D6_cellbin.final.celltype.h5ad"
adata = sc.read_h5ad(adata)

sc.pl.umap(adata_concat,color="log",frameon=False, na_color="grey",save=".log.cell.png")

