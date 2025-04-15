import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import harmonypy as hm
import pandas as pd

from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import matplotlib as mpl

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/05.cluster/plot_marker")

adata_concat = sc.read_h5ad("pancancer.ref.0807.final.CD8.umap.h5ad")
sc.pl.umap(adata_concat, color=[
                                                   'CD4','CD8A','CD8B','SLC4A10','DPP4','TRDV2','TRGV9','TRGV10','FCGR3A','NCAM1'
                                                     ],
               s=0.6, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".CD8.umap.markergene.png")
sc.pl.umap(adata_concat, color=[
                                                   'CD4','CD8A','CD8B','SLC4A10','DPP4','TRDV2','TRGV9','TRGV10','FCGR3A','NCAM1'
                                                     ],
               s=0.6, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".CD8.umap.markergene.pdf")

adata_concat = sc.read_h5ad("pancancer.ref.0807.final.CD4.umap.h5ad")
sc.pl.umap(adata_concat, color=[
                                                   'CD4','CD8A','CD8B','SLC4A10','DPP4','TRGV9','TRGV10','FCGR3A','NCAM1'
                                                     ],
               s=0.6, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".CD4.umap.markergene.png")
sc.pl.umap(adata_concat, color=[
                                                   'CD4','CD8A','CD8B','SLC4A10','DPP4','TRGV9','TRGV10','FCGR3A','NCAM1'
                                                     ],
               s=0.6, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".CD4.umap.markergene.pdf")

adata_concat = sc.read_h5ad("pancancer.ref.0807.final.NK.umap.h5ad")
sc.pl.umap(adata_concat, color=[
                                                   'CD4','CD8A','CD8B','SLC4A10','DPP4','TRDV2','TRGV9','TRGV10','FCGR3A','NCAM1'
                                                     ],
               s=0.6, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".NK.umap.markergene.png")
sc.pl.umap(adata_concat, color=[
                                                   'CD4','CD8A','CD8B','SLC4A10','DPP4','TRDV2','TRGV9','TRGV10','FCGR3A','NCAM1'
                                                     ],
               s=0.6, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".NK.umap.markergene.pdf")


