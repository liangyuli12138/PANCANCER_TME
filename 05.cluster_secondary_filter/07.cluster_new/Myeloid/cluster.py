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


os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/05.cluster_secondary_filter/07.cluster_new/Myeloid")
adata_concat = sc.read_h5ad("pancancer.ref.final.final.Myeloid.umap.old.h5ad")
adata_concat = adata_concat.raw.to_adata()

cellist = pd.read_csv("Myeloid.input")
atlist = pd.read_csv("Myeloid.at",index_col=0)
adata_concat = adata_concat[cellist["cell"],:]
adata_concat.obs = adata_concat.obs.join(atlist)

adata_concat.write_h5ad("pancancer.ref.final.final.Myeloid.umap.h5ad")
adata_concat.obs.to_csv("pancancer.ref.final.final.Myeloid.umap.obs.csv")
adata_concat.var.to_csv("pancancer.ref.final.final.Myeloid.umap.var.csv")


sc.pl.umap(adata_concat, color='groups_refs', size=1.5, save= ".cluster.png")
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
[axMyeloid_Mono,axMyeloid_cDC1,axMyeloid_cDC2],[axMyeloid_cDC3,axMyeloid_pDC,axMyeloid_Marco_LYVE1],[axMyeloid_Marco_C1QC,axMyeloid_Marco_SPP1,axMyeloid_Mast]     
) = plt.subplots(3, 3, figsize=(12,9),
     )

sc.pl.umap(adata_concat, color='groups_refs', groups=['Myeloid_Mono'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axMyeloid_Mono)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Myeloid_cDC1'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axMyeloid_cDC1)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Myeloid_cDC2'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axMyeloid_cDC2)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Myeloid_cDC3'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axMyeloid_cDC3)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Myeloid_pDC'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axMyeloid_pDC)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Myeloid_Marco_LYVE1'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axMyeloid_Marco_LYVE1)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Myeloid_Marco_C1QC'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axMyeloid_Marco_C1QC)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Myeloid_Marco_SPP1'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axMyeloid_Marco_SPP1)


sc.pl.umap(adata_concat, color='groups_refs', groups=['Myeloid_Mast'],  legend_loc='on data',
           legend_fontweight="light",legend_fontsize=8,size=1.5, legend_fontoutline=0, ax=axMyeloid_Mast)



plt.savefig("figures/sub.celltype.0809.pdf")
plt.savefig("figures/sub.celltype.0809.png")

sc.pl.umap(adata_concat, color=[
"LYVE1","FOLR2","NLRP3","CD163","MRC1","TIMD4","APOE","APOC1","ACP5","FABP4","FABP5","LIPA","GPNMB","ISG15","CXCL10","IL1B","CCL3L1","CCL3L1","CXCL1","CXCL3","CXCL8",
"C1QA","C1QB","C1QC","MARCO","TREM2","VEGFA","SPP1","VCAN","FCN1","FN1","INHBA","THBS1","CLEC9A","CADM1","CLNK","CD1C","CLEC10A","FCER1A","LAMP3","EBI3","CCL19",
"GZMB","LILRA4","CLEC4C","FCGR3A","EREG","S100A8","S100A9"
],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".Myeloid.umap.markergene.png")
sc.pl.umap(adata_concat, color=[
"LYVE1","FOLR2","NLRP3","CD163","MRC1","TIMD4","APOE","APOC1","ACP5","FABP4","FABP5","LIPA","GPNMB","ISG15","CXCL10","IL1B","CCL3L1","CCL3L1","CXCL1","CXCL3","CXCL8",
"C1QA","C1QB","C1QC","MARCO","TREM2","VEGFA","SPP1","VCAN","FCN1","FN1","INHBA","THBS1","CLEC9A","CADM1","CLNK","CD1C","CLEC10A","FCER1A","LAMP3","EBI3","CCL19",
"GZMB","LILRA4","CLEC4C","FCGR3A","EREG","S100A8","S100A9"                     
],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".Myeloid.umap.markergene.pdf")
