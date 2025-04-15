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


os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/07.cluster_new/B")
adata_concat = sc.read_h5ad("pancancer.ref.0807.final.B.umap.h5ad")

cellist = pd.read_csv("B.input")
atlist = pd.read_csv("B.at",index_col=0)

adata_concat = adata_concat[cellist["cell"],:]
adata_concat.obs = adata_concat.obs.join(atlist)

adata_concat.write_h5ad("pancancer.ref.final.final.B.umap.h5ad")
adata_concat.obs.to_csv("pancancer.ref.final.final.B.umap.obs.csv")
adata_concat.var.to_csv("pancancer.ref.final.final.B.umap.var.csv")

sc.pl.umap(adata_concat, color='groups_refs', size=1.5, save= ".cluster.png")
sc.pl.umap(adata_concat, color='groups_refs', legend_fontweight="light",legend_fontsize=10,size=1.5, legend_fontoutline=0, 
           legend_loc='on data', save= ".on.cluster.png")
sc.pl.umap(adata_concat, color="Tissue", size=1.5, save=".Tissue.cluster.png")
sc.pl.umap(adata_concat, color="Phenotype", size=1.5, save=".Phenotype.cluster.png")
sc.pl.umap(adata_concat, color='groups_refs', legend_fontweight="light",legend_fontsize=10,size=1.5, legend_fontoutline=0,            
           legend_loc='on data', save= ".on.cluster.pdf")
sc.pl.umap(adata_concat, color="Tissue", size=1.5, save=".Tissue.cluster.pdf")
sc.pl.umap(adata_concat, color="Phenotype", size=1.5, save=".Phenotype.cluster.pdf")

fig, (
[axLymphoid_B_memory,axLymphoid_B_naive],[axLymphoid_Plamsa_IGLC,axLymphoid_Plamsa_IGKC]     
) = plt.subplots(2, 2, figsize=(8,6),
     )

sc.pl.umap(adata_concat, color='groups_refs', groups=['Lymphoid_B_memory'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_B_memory)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Lymphoid_B_naive'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_B_naive)

sc.pl.umap(adata_concat, color='groups_refs', groups=['Lymphoid_Plamsa_IGLC'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_Plamsa_IGLC)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Lymphoid_Plamsa_IGKC'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axLymphoid_Plamsa_IGKC)


plt.savefig("figures/sub.celltype.0809.pdf")
plt.savefig("figures/sub.celltype.0809.png")

sc.pl.umap(adata_concat, color=[
"YBX3","CD37","TCL1A","KLF2","PLCG2","CD79B","CD74","IGHD","SELL","YBX3","IL4R","FCER2","TCL1A","IL4R","FCER2","LTB","GPR183","HLA-DRA","HLA-DPA1","HLA-DPB1","HLA-DRB1","BCL2A1","DUSP2","MYC","CD83","EGR3","REL","MIR155HG","CRIP1","ANXA2","S100A10","VIM","TNFRSF13B","HSPA1A","HSPA1B","DNAJA1","DNAJB1","NR4A2","ZBTB32","MYC","IL2RA","IRF4","TNFRSF1B","RGS13","NEIL1","RFTN1","CDCA7","AICDA","IGHA1","IGHA2","IGHM","IGLC2","IGLC3","IGLL5","JCHAIN","IGHG1","IGHG2","IGHG3","IGHG4","IGHGP","IGLC7","IGKC","IGHA1","IGHA2","PRDM1","XBP1","MZB1","SDC1","IGHG1","IGHG2","PRDM1","XBP1","MZB1","SDC1","IGHG1","IGHG4","MZB1","XBP1","JCHAIN"                                                     
],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".B.umap.markergene.png")
sc.pl.umap(adata_concat, color=[
"YBX3","CD37","TCL1A","KLF2","PLCG2","CD79B","CD74","IGHD","SELL","YBX3","IL4R","FCER2","TCL1A","IL4R","FCER2","LTB","GPR183","HLA-DRA","HLA-DPA1","HLA-DPB1","HLA-DRB1","BCL2A1","DUSP2","MYC","CD83","EGR3","REL","MIR155HG","CRIP1","ANXA2","S100A10","VIM","TNFRSF13B","HSPA1A","HSPA1B","DNAJA1","DNAJB1","NR4A2","ZBTB32","MYC","IL2RA","IRF4","TNFRSF1B","RGS13","NEIL1","RFTN1","CDCA7","AICDA","IGHA1","IGHA2","IGHM","IGLC2","IGLC3","IGLL5","JCHAIN","IGHG1","IGHG2","IGHG3","IGHG4","IGHGP","IGLC7","IGKC","IGHA1","IGHA2","PRDM1","XBP1","MZB1","SDC1","IGHG1","IGHG2","PRDM1","XBP1","MZB1","SDC1","IGHG1","IGHG4","MZB1","XBP1","JCHAIN"                     
],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".B.umap.markergene.pdf")
