import scanpy as sc
import pandas as pd

adata = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/05.stat/leiden/immune.cluster.r0.5.h5ad"
adata = sc.read_h5ad(adata)
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/05.stat/l2m/l2m.stat.merge.csv.at",index_col=0)
adata.obs = adata.obs.join(atlist)


sc.pl.umap(adata,color="Lymphoid_percentage",frameon=False, na_color="grey",save=".Lymphoid_percentage.umap.png")
sc.pl.umap(adata,color="Myeloid_percentage",frameon=False, na_color="grey",save=".Myeloid_percentage.umap.png")
sc.pl.umap(adata,color="log10_Lym_Mye_ratio",frameon=False, na_color="grey",vmin=-2, vmax=2,save=".log10_Lym_Mye_ratio.umap.png")
sc.pl.umap(adata,color="density",frameon=False, na_color="grey",save=".density.umap.png")
sc.pl.umap(adata,color="area",frameon=False, na_color="grey",vmax=0.5,save=".area.umap.png")
sc.pl.umap(adata,color="elongation",frameon=False, na_color="grey",save=".elongation.umap.png")
