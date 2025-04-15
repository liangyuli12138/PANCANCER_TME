import scanpy as sc
import pandas as pd
import os

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/05.l2m")

adata = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/05.l2m/immune.cluster.11.r1.5.new.h5ad"
adata = sc.read_h5ad(adata)
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/05.l2m/l2m.stat.merge.csv.at",index_col=0)
adata.obs = adata.obs.join(atlist)


sc.pl.umap(adata,color="Lymphoid_percentage",frameon=False, na_color="grey",save=".Lymphoid_percentage.umap.pdf")
sc.pl.umap(adata,color="Myeloid_percentage",frameon=False, na_color="grey",save=".Myeloid_percentage.umap.pdf")
sc.pl.umap(adata,color="log10_Lym2Mye_ratio",frameon=False, na_color="grey",vmin=-2, vmax=2,save=".log10_Lym2Mye_ratio.umap.pdf")
sc.pl.umap(adata,color="density",frameon=False, na_color="grey",save=".density.umap.pdf")
sc.pl.umap(adata,color="area",frameon=False, na_color="grey",vmax=0.5,save=".area.umap.pdf")
sc.pl.umap(adata,color="elongation",frameon=False, na_color="grey",save=".elongation.umap.pdf")
sc.pl.umap(adata,color="distance",frameon=False, na_color="grey",vmin=-4000, vmax=4000,save=".distance.umap.pdf")
