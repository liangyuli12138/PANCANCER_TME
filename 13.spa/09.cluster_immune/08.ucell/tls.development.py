import scanpy as sc
import pandas as pd

adata = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/05.stat/leiden/immune.cluster.r0.5.h5ad"
adata = sc.read_h5ad(adata)
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/08.ucell/at/tls.development.gene.at",index_col=0)
adata.obs = adata.obs.join(atlist)


sc.pl.umap(adata,color="early",frameon=False, na_color="grey",save=".early.gene.png")
sc.pl.umap(adata,color="intermediate",frameon=False, na_color="grey",save=".intermediate.gene.png")
sc.pl.umap(adata,color="Late",frameon=False, na_color="grey",save=".Late.gene.png")
sc.pl.umap(adata,color="GO0002317",frameon=False, na_color="grey",save=".GO0002317.gene.png")

import matplotlib.pyplot as plt
sc.pl.violin(adata, ["early","intermediate","Late","GO0002317"], groupby='leiden')
plt.savefig('violin_plot.png')

