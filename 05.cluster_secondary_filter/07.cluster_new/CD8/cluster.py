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


os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/07.cluster_new/CD8")
adata_concat = sc.read_h5ad("pancancer.ref.0807.final.CD8_pca.umap.h5ad")

cellist = pd.read_csv("CD8.input")
atlist = pd.read_csv("CD8.at",index_col=0)

adata_concat = adata_concat[cellist["cell"],:]
adata_concat.obs = adata_concat.obs.join(atlist)

adata_concat.write_h5ad("pancancer.ref.final.final.CD8.umap.h5ad")
adata_concat.obs.to_csv("pancancer.ref.final.final.CD8.umap.obs.csv")
adata_concat.var.to_csv("pancancer.ref.final.final.CD8.umap.var.csv")

sc.pl.umap(adata_concat, color='groups_ref', size=1.5, save= ".cluster.png")
sc.pl.umap(adata_concat, color='groups_ref', legend_fontweight="light",legend_fontsize=10,size=1.5, legend_fontoutline=0, 
           legend_loc='on data', save= ".on.cluster.png")
sc.pl.umap(adata_concat, color="Tissue", size=1.5, save=".Tissue.cluster.png")
sc.pl.umap(adata_concat, color="Phenotype", size=1.5, save=".Phenotype.cluster.png")
sc.pl.umap(adata_concat, color='groups_ref', legend_fontweight="light",legend_fontsize=10,size=1.5, legend_fontoutline=0,            
           legend_loc='on data', save= ".on.cluster.pdf")
sc.pl.umap(adata_concat, color="Tissue", size=1.5, save=".Tissue.cluster.pdf")
sc.pl.umap(adata_concat, color="Phenotype", size=1.5, save=".Phenotype.cluster.pdf")

fig, ([axLymphoid_CD8_Tn,axLymphoid_CD8_Tm,axLymphoid_CD8_Teff],[axLymphoid_CD8_Tex,axLymphoid_CD8_Tisg,axLymphoid_CD8_Tstr] 
     ) = plt.subplots(2, 3, figsize=(10,6),
     )

sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_CD8_Tn'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_CD8_Tn)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_CD8_Tm'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_CD8_Tm)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_CD8_Teff'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_CD8_Teff)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_CD8_Tex'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_CD8_Tex)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_CD8_Tisg'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_CD8_Tisg)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_CD8_Tstr'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_CD8_Tstr)



plt.savefig("figures/sub.celltype.0809.pdf")
plt.savefig("figures/sub.celltype.0809.png")

sc.pl.umap(adata_concat, color=[
                                                   'CD8A','CD8B','TCF7','CCR7','SELL','ITGA1','GZMK','CD69','PDCD1','CTLA4','LAG3','FOSB','JUN','IFIT1','MX1'
                                                     ],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".CD8.umap.markergene.png")
sc.pl.umap(adata_concat, color=[
                                                   'CD8A','CD8B','TCF7','CCR7','SELL','ITGA1','GZMK','CD69','PDCD1','CTLA4','LAG3','FOSB','JUN','IFIT1','MX1'
                                                     ],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".CD8.umap.markergene.pdf")
