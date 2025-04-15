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


os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/07.cluster_new/EC")
adata_concat = sc.read_h5ad("pancancer.ref.0807.final.EC_pca.umap.h5ad")

cellist = pd.read_csv("EC.input")
atlist = pd.read_csv("EC.at",index_col=0)

adata_concat = adata_concat[cellist["cell"],:]
adata_concat.obs = adata_concat.obs.join(atlist)

adata_concat.write_h5ad("pancancer.ref.final.final.EC.umap.h5ad")
adata_concat.obs.to_csv("pancancer.ref.final.final.EC.umap.obs.csv")
adata_concat.var.to_csv("pancancer.ref.final.final.EC.umap.var.csv")

sc.pl.umap(adata_concat, color='groups_ref', palette = 'Paired_r', size=1.5, save= ".cluster.png")
sc.pl.umap(adata_concat, color='groups_ref', palette = 'Paired_r', legend_fontweight="light",legend_fontsize=10,size=1.5, legend_fontoutline=0, 
           legend_loc='on data', save= ".on.cluster.png")
sc.pl.umap(adata_concat, color="Tissue", size=1.5, save=".Tissue.cluster.png")
sc.pl.umap(adata_concat, color="Phenotype", size=1.5, save=".Phenotype.cluster.png")
sc.pl.umap(adata_concat, color='groups_ref', palette = 'Paired_r', legend_fontweight="light",legend_fontsize=10,size=1.5, legend_fontoutline=0,            
           legend_loc='on data', save= ".on.cluster.pdf")
sc.pl.umap(adata_concat, color="Tissue", size=1.5, save=".Tissue.cluster.pdf")
sc.pl.umap(adata_concat, color="Phenotype", size=1.5, save=".Phenotype.cluster.pdf")

fig, (
     [axEC_Vein,axEC_Capillary,axEC_Angiogenic],[axEC_Artery,axEC_Alveolar,axEC_Glomerular],[axEC_Sinusoidal,axEC_Lymph,xxx] 
     ) = plt.subplots(3, 3, figsize=(10,9),
     )
sc.pl.umap(adata_concat, color='groups_ref', groups=['EC_Vein'],  legend_loc='on data', palette = 'Paired_r', 
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axEC_Vein)


sc.pl.umap(adata_concat, color='groups_ref', groups=['EC_Capillary'],  legend_loc='on data', palette = 'Paired_r', 
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axEC_Capillary)


sc.pl.umap(adata_concat, color='groups_ref', groups=['EC_Angiogenic'],  legend_loc='on data', palette = 'Paired_r', 
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axEC_Angiogenic)


sc.pl.umap(adata_concat, color='groups_ref', groups=['EC_Artery'],  legend_loc='on data', palette = 'Paired_r', 
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axEC_Artery)


sc.pl.umap(adata_concat, color='groups_ref', groups=['EC_Alveolar'],  legend_loc='on data', palette = 'Paired_r', 
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axEC_Alveolar)


sc.pl.umap(adata_concat, color='groups_ref', groups=['EC_Glomerular'],  legend_loc='on data', palette = 'Paired_r', 
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axEC_Glomerular)


sc.pl.umap(adata_concat, color='groups_ref', groups=['EC_Sinusoidal'],  legend_loc='on data', palette = 'Paired_r', 
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axEC_Sinusoidal)


sc.pl.umap(adata_concat, color='groups_ref', groups=['EC_Lymph'],  legend_loc='on data', palette = 'Paired_r', 
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axEC_Lymph)


plt.savefig("figures/sub.celltype.0809.pdf")
plt.savefig("figures/sub.celltype.0809.png")

sc.pl.umap(adata_concat, color=[
                               'ACKR1','VCAM1','CA4','CD36','INSR','ESM1','CXCL12','IGFBP3'
                                                     ],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".EC.umap.markergene.png")
sc.pl.umap(adata_concat, color=[
                                'ACKR1','VCAM1','CA4','CD36','INSR','ESM1','CXCL12','IGFBP3'                     
                                ],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".EC.umap.markergene.pdf")
