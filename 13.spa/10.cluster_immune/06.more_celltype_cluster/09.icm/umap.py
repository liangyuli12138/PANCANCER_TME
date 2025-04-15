import scanpy as sc
import pandas as pd
import os

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/09.icm")

adata = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/05.l2m/immune.cluster.11.r1.5.new.h5ad"
adata = sc.read_h5ad(adata)
atlist = pd.read_csv("l2m.stat.merge.csv.at.icm",index_col=0)
adata.obs = adata.obs.join(atlist)


sc.pl.umap(adata,color="ICM1_percent",frameon=False, na_color="grey",save=".ICM1_percentage.umap.pdf")
sc.pl.umap(adata,color="ICM2_percent",frameon=False, na_color="grey",save=".ICM2_percentage.umap.pdf")
sc.pl.umap(adata,color="ICM3_percent",frameon=False, na_color="grey",save=".ICM3_percentage.umap.pdf")
