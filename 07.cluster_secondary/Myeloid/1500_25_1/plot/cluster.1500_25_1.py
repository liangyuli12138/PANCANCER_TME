import os
import datetime
import json
import argparse

import anndata
import scanpy as sc
import harmonypy as hm
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/07.cluster_secondary/Myeloid/1500_25_1/plot")

adata_concat = sc.read_h5ad("pancancer.ref.0807.final.Myeloid.umap.h5ad")


sc.pl.umap(adata_concat, color='groups_secori',  legend_loc='on data', save= ".Myeloid.1500_25_1.on.cluster.png")
sc.pl.umap(adata_concat, color="Tissue", save=".Myeloid.1500_25_1.Tissue.cluster.png")
sc.pl.umap(adata_concat, color="Phenotype", save=".Myeloid.1500_25_1.Phenotype.cluster.png")

sc.pl.umap(adata_concat, color=[
"LYVE1","FOLR2","NLRP3","CD163","MRC1","TIMD4","APOE","APOC1","ACP5","FABP4","FABP5","LIPA","GPNMB","ISG15","CXCL10","IL1B","CCL3L1","CCL3L1","CXCL1","CXCL3","CXCL8","C1QA","C1QB","C1QC","MARCO","TREM2","VEGFA","SPP1","VCAN","FCN1","FN1","INHBA","THBS1"
],
               s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save=".Myeloid.1500_25_1.Myeloid.umap.markergene.png")
