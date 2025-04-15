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


os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/07.cluster_new/CD4")
adata_concat = sc.read_h5ad("pancancer.ref.0904.final.CD4_pca.umap.h5ad")

cellist = pd.read_csv("CD4.input")
atlist = pd.read_csv("CD4.at",index_col=0)

adata_concat = adata_concat[cellist["cell"],:]
adata_concat.obs = adata_concat.obs.join(atlist)

adata_concat.write_h5ad("pancancer.ref.final.final.CD4.umap.h5ad")
adata_concat.obs.to_csv("pancancer.ref.final.final.CD4.umap.obs.csv")
adata_concat.var.to_csv("pancancer.ref.final.final.CD4.umap.var.csv")

sc.pl.umap(adata_concat, color='groups_ref', size=1.5, save= ".cluster.png")
sc.pl.umap(adata_concat, color='groups_ref', legend_fontweight="light",legend_fontsize=10,size=1.5, legend_fontoutline=0, 
           legend_loc='on data', save= ".on.cluster.png")
sc.pl.umap(adata_concat, color="Tissue", size=1.5, save=".Tissue.cluster.png")
sc.pl.umap(adata_concat, color="Phenotype", size=1.5, save=".Phenotype.cluster.png")
sc.pl.umap(adata_concat, color='groups_ref', legend_fontweight="light",legend_fontsize=10,size=1.5, legend_fontoutline=0,            
           legend_loc='on data', save= ".on.cluster.pdf")
sc.pl.umap(adata_concat, color="Tissue", size=1.5, save=".Tissue.cluster.pdf")
sc.pl.umap(adata_concat, color="Phenotype", size=1.5, save=".Phenotype.cluster.pdf")

fig, ([axLymphoid_CD4_Tn,axLymphoid_CD4_Tcm,axLymphoid_CD4_Tfh],[axLymphoid_CD4_Th17,axLymphoid_CD4_CTL,axLymphoid_CD4_Treg],[axLymphoid_CD4_Tstr,xxx,yyy] 
     ) = plt.subplots(3, 3, figsize=(10,9),
     )

sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_CD4_Tn'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_CD4_Tn)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_CD4_Tcm'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_CD4_Tcm)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_CD4_Tfh'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_CD4_Tfh)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_CD4_Th17'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_CD4_Th17)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_CD4_CTL'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_CD4_CTL)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_CD4_Treg'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_CD4_Treg)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Lymphoid_CD4_Tstr'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_CD4_Tstr)


plt.savefig("figures/sub.celltype.0809.pdf")
plt.savefig("figures/sub.celltype.0809.png")

sc.pl.umap(adata_concat, color=[
                                                   'CD4','GPR183','KLRB1','CXCL13','TOX','IFNG','GZMK','FOXP3','IL2RA','CTLA4','FOS','JUN','GZMA','IL17A','RORA',
                                                     ],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".CD4.umap.markergene.png")
sc.pl.umap(adata_concat, color=[
                                                   'CD4','GPR183','KLRB1','CXCL13','TOX','IFNG','GZMK','FOXP3','IL2RA','CTLA4','FOS','JUN','GZMA','IL17A','RORA',
                                                     ],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".CD4.umap.markergene.pdf")

