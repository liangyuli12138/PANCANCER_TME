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


os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/07.cluster_new/NK")
adata_concat = sc.read_h5ad("pancancer.ref.0807.final.NK.umap.h5ad")

cellist = pd.read_csv("NK.input")
atlist = pd.read_csv("NK.at",index_col=0)

adata_concat = adata_concat[cellist["cell"],:]
adata_concat.obs = adata_concat.obs.join(atlist)

adata_concat.write_h5ad("pancancer.ref.final.final.NK.umap.h5ad")
adata_concat.obs.to_csv("pancancer.ref.final.final.NK.umap.obs.csv")
adata_concat.var.to_csv("pancancer.ref.final.final.NK.umap.var.csv")

sc.pl.umap(adata_concat, color='groups_ref', size=1.5, save= ".cluster.png")
sc.pl.umap(adata_concat, color='groups_ref', legend_fontweight="light",legend_fontsize=10,size=1.5, legend_fontoutline=0, 
           legend_loc='on data', save= ".on.cluster.png")
sc.pl.umap(adata_concat, color="Tissue", size=1.5, save=".Tissue.cluster.png")
sc.pl.umap(adata_concat, color="Phenotype", size=1.5, save=".Phenotype.cluster.png")
sc.pl.umap(adata_concat, color='groups_ref', legend_fontweight="light",legend_fontsize=10,size=1.5, legend_fontoutline=0,            
           legend_loc='on data', save= ".on.cluster.pdf")
sc.pl.umap(adata_concat, color="Tissue", size=1.5, save=".Tissue.cluster.pdf")
sc.pl.umap(adata_concat, color="Phenotype", size=1.5, save=".Phenotype.cluster.pdf")

fig, ([axLymphoid_ILC,axLymphoid_MAIT,axLymphoid_NK_CD16__],[axLymphoid_NK_CD56__,axLymphoid_NKT,xxx]      
     ) = plt.subplots(2, 3, figsize=(10,6),
     )

sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_ILC'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_ILC)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_MAIT'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_MAIT)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_NK_CD16+'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_NK_CD16__)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_NK_CD56+'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_NK_CD56__)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_NKT'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_NKT)

plt.savefig("figures/sub.celltype.0809.pdf")
plt.savefig("figures/sub.celltype.0809.png")

sc.pl.umap(adata_concat, color=[
                                                   'FCGR3A','CD8A','NCAM1','SLC4A10','DPP4'
                                                     ],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".NK.umap.markergene.png")
sc.pl.umap(adata_concat, color=[
                                                   'FCGR3A','CD8A','NCAM1','SLC4A10','DPP4'
                                                     ],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".NK.umap.markergene.pdf")

