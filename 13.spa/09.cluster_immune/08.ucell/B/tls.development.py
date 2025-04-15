import scanpy as sc
import pandas as pd
import os

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/08.ucell/B")
adata = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/05.stat/leiden/immune.cluster.r0.5.h5ad"
adata = sc.read_h5ad(adata)
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/08.ucell/at/tls.development.gene.at.B",index_col=0)
adata.obs = adata.obs.join(atlist)


sc.pl.umap(adata,color="early",frameon=False, na_color="grey",save=".B.early.gene.png")
sc.pl.umap(adata,color="intermediate",frameon=False, na_color="grey",save=".B.intermediate.gene.png")
sc.pl.umap(adata,color="Late",frameon=False, na_color="grey",save=".B.Late.gene.png")
sc.pl.umap(adata,color="GO0002317",frameon=False, na_color="grey",save=".B.GO0002317.gene.png")

import matplotlib.pyplot as plt
sc.pl.violin(adata, ["early","intermediate","Late","GO0002317"], groupby='leiden')
plt.savefig('filter/B.violin_plot.png')

