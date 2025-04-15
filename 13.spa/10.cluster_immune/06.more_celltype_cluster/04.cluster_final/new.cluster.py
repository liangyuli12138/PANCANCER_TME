import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/04.cluster_final")
adata = sc.read_h5ad("immune.cluster.11.r1.5.h5ad")
atlist = pd.read_csv("new.group.list.at",index_col=0)
adata.obs = adata.obs.join(atlist)

sc.pl.umap(adata, color=['new_groups'], size=40, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10,save= "immune.cluster.11.r1.5.png")
sc.pl.umap(adata, color=['new_groups'], size=40, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10,save= "immune.cluster.11.r1.5.pdf")
sc.pl.umap(adata, color=['new_groups'], size=40, color_map = 'RdPu', ncols = 2, legend_fontsize=10,save= "immune.cluster.11.r1.5.loc.png")
sc.pl.umap(adata, color=['new_groups'], size=40, color_map = 'RdPu', ncols = 2, legend_fontsize=10,save= "immune.cluster.11.r1.5.loc.pdf")

adata.write_h5ad("immune.cluster.11.r1.5.new.h5ad")
adata.obs.to_csv("immune.cluster.11.r1.5.new.obs")

