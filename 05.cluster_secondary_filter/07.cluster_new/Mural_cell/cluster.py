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


os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/07.cluster_new/Mural_cell")
adata_concat = sc.read_h5ad("pancancer.ref.0723.final.Mural_cell.umap.h5ad")

cellist = pd.read_csv("Mural_cell.input")
atlist = pd.read_csv("Mural_cell.at",index_col=0)

adata_concat = adata_concat[cellist["cell"],:]
adata_concat.obs = adata_concat.obs.join(atlist)

adata_concat.write_h5ad("pancancer.ref.final.final.Mural_cell.umap.h5ad")
adata_concat.obs.to_csv("pancancer.ref.final.final.Mural_cell.umap.obs.csv")
adata_concat.var.to_csv("pancancer.ref.final.final.Mural_cell.umap.var.csv")

sc.pl.umap(adata_concat, color='groups_ref',  size=1.5, save= ".cluster.png")
sc.pl.umap(adata_concat, color='groups_ref',  legend_fontweight="light",legend_fontsize=10,size=1.5, legend_fontoutline=0, 
           legend_loc='on data', save= ".on.cluster.png")
sc.pl.umap(adata_concat, color="Tissue", size=1.5, save=".Tissue.cluster.png")
sc.pl.umap(adata_concat, color="Phenotype", size=1.5, save=".Phenotype.cluster.png")
sc.pl.umap(adata_concat, color='groups_ref',  legend_fontweight="light",legend_fontsize=10,size=1.5, legend_fontoutline=0,            
           legend_loc='on data', save= ".on.cluster.pdf")
sc.pl.umap(adata_concat, color="Tissue", size=1.5, save=".Tissue.cluster.pdf")
sc.pl.umap(adata_concat, color="Phenotype", size=1.5, save=".Phenotype.cluster.pdf")

fig, (
     [axMural_cell_SMC1,axMural_cell_SMC2],[axMural_cell_Pericyte1,axMural_cell_Pericyte2]
     ) = plt.subplots(2, 2, figsize=(10,8),
     )

sc.pl.umap(adata_concat, color='groups_ref', groups=['Mural_cell_SMC1'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axMural_cell_SMC1)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Mural_cell_SMC2'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axMural_cell_SMC2)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Mural_cell_Pericyte1'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axMural_cell_Pericyte1)


sc.pl.umap(adata_concat, color='groups_ref', groups=['Mural_cell_Pericyte2'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axMural_cell_Pericyte2)

plt.savefig("figures/sub.celltype.0809.pdf")
plt.savefig("figures/sub.celltype.0809.png")

sc.pl.umap(adata_concat, color=[
      "MYH11","NDUFA4","MUSTN1","PLN","RERGL","CASQ2","ACTG2","DES","CLU","RAMP1","NCAM1","FBLN1","NOTCH3","PDGFRB","HIGD1B","MYO1B","COL4A1","THY1","CCL2",
      "CD44","CXCL12","C1R","C1S","CFD"
                                                     ],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".Mural_cell.umap.markergene.png")
sc.pl.umap(adata_concat, color=[
      "MYH11","NDUFA4","MUSTN1","PLN","RERGL","CASQ2","ACTG2","DES","CLU","RAMP1","NCAM1","FBLN1","NOTCH3","PDGFRB","HIGD1B","MYO1B","COL4A1","THY1","CCL2",      "CD44","CXCL12","C1R","C1S","CFD"
                                ],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".Mural_cell.umap.markergene.pdf")
