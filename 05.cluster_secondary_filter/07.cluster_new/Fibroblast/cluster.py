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


os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/07.cluster_new/Fibroblast")
adata_concat = sc.read_h5ad("pancancer.ref.0807.final.Fibroblast.umap.h5ad")

cellist = pd.read_csv("Fibroblast.input")
atlist = pd.read_csv("Fibroblast.at",index_col=0)

adata_concat = adata_concat[cellist["cell"],:]
adata_concat.obs = adata_concat.obs.join(atlist)

adata_concat.write_h5ad("pancancer.ref.final.final.Fibroblast.umap.h5ad")
adata_concat.obs.to_csv("pancancer.ref.final.final.Fibroblast.umap.obs.csv")
adata_concat.var.to_csv("pancancer.ref.final.final.Fibroblast.umap.var.csv")

sc.pl.umap(adata_concat, color='groups_refs', palette = 'Paired_r', size=1.5, save= ".cluster.png")
sc.pl.umap(adata_concat, color='groups_refs', palette = 'Paired_r', legend_fontweight="light",legend_fontsize=10,size=1.5, legend_fontoutline=0, 
           legend_loc='on data', save= ".on.cluster.png")
sc.pl.umap(adata_concat, color="Tissue", size=1.5, save=".Tissue.cluster.png")
sc.pl.umap(adata_concat, color="Phenotype", size=1.5, save=".Phenotype.cluster.png")
sc.pl.umap(adata_concat, color='groups_refs', palette = 'Paired_r', legend_fontweight="light",legend_fontsize=10,size=1.5, legend_fontoutline=0,            
           legend_loc='on data', save= ".on.cluster.pdf")
sc.pl.umap(adata_concat, color="Tissue", size=1.5, save=".Tissue.cluster.pdf")
sc.pl.umap(adata_concat, color="Phenotype", size=1.5, save=".Phenotype.cluster.pdf")

fig, (
[axFibroblast_iCAF_IGFBP6,axFibroblast_iCAF_IL6,axFibroblast_iCAF_KCNN3],[axFibroblast_mCAF_KRT19,axFibroblast_mCAF_MMP11,axFibroblast_mCAF_WNT5A],[axFibroblast_apCAF_CD74,xxx,yyy]
) = plt.subplots(3, 3, figsize=(12,9),

     )

sc.pl.umap(adata_concat, color='groups_refs', groups=['Fibroblast_iCAF_IGFBP6'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axFibroblast_iCAF_IGFBP6)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Fibroblast_iCAF_IL6'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axFibroblast_iCAF_IL6)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Fibroblast_iCAF_KCNN3'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axFibroblast_iCAF_KCNN3)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Fibroblast_mCAF_KRT19'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axFibroblast_mCAF_KRT19)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Fibroblast_mCAF_POSTN'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axFibroblast_mCAF_MMP11)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Fibroblast_mCAF_WNT5A'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axFibroblast_mCAF_WNT5A)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Fibroblast_apCAF_CD74'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axFibroblast_apCAF_CD74)

plt.savefig("figures/sub.celltype.0809.pdf")
plt.savefig("figures/sub.celltype.0809.png")

sc.pl.umap(adata_concat, color=[
"CFD","IGFBP6","MFAP5","FBLN2","CD34","TGFBI","FHL2","WNT5A","ITGA8","POSTN","BGN","COL3A1","INHBA","COL10A1","MMP11",
"C7","KCNN3","DPT","SCN7A","P2RY1","IL6","CXCL2","APOC1",
"APOE","CFD","CD74","HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","MYLK","KRT19","MYH11","NCAM1","LEFTY2"
],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".Fibroblast.umap.markergene.png")
