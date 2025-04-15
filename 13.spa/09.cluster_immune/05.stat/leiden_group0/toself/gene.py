import scanpy as sc
import pandas as pd

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/05.stat/leiden_group0/toself")
adata = "immune.group0.r0.5.h5ad"
adata = sc.read_h5ad(adata)
atlist = pd.read_csv("all.toself.stat.csv.log.group0.at",index_col=0)
adata.obs = adata.obs.join(atlist)

sc.pl.umap(adata,color="-log10",frameon=False, na_color="grey",save=".group0.completeness.cell.png")

