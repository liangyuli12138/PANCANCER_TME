import scanpy as sc
import pandas as pd

adata = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/05.stat/leiden/immune.cluster.r0.5.h5ad"
adata = sc.read_h5ad(adata)
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/gene/merge.gene.csv.at",index_col=0)
adata.obs = adata.obs.join(atlist)

#sc.pl.umap(adata,color="log",frameon=False, na_color="grey",save=".log.cell.png")

sc.pl.umap(adata,color="CD19",frameon=False, na_color="grey",save=".CD19.gene.png")
sc.pl.umap(adata,color="MS4A1",frameon=False, na_color="grey",save=".MS4A1.gene.png")
sc.pl.umap(adata,color="CD3D",frameon=False, na_color="grey",save=".CD3D.gene.png")
sc.pl.umap(adata,color="CXCL9",frameon=False, na_color="grey",save=".CXCL9.gene.png")
sc.pl.umap(adata,color="CXCL10",frameon=False, na_color="grey",save=".CXCL10.gene.png")
sc.pl.umap(adata,color="CXCL11",frameon=False, na_color="grey",save=".CXCL11.gene.png")
sc.pl.umap(adata,color="CXCL12",frameon=False, na_color="grey",save=".CXCL12.gene.png")
sc.pl.umap(adata,color="CXCL13",frameon=False, na_color="grey",save=".CXCL13.gene.png")
sc.pl.umap(adata,color="CXCR5",frameon=False, na_color="grey",save=".CXCR5.gene.png")
sc.pl.umap(adata,color="ICOS",frameon=False, na_color="grey",save=".ICOS.gene.png")
sc.pl.umap(adata,color="CTLA4",frameon=False, na_color="grey",save=".CTLA4.gene.png")
sc.pl.umap(adata,color="PDCD1",frameon=False, na_color="grey",save=".PDCD1.gene.png")
sc.pl.umap(adata,color="LTA",frameon=False, na_color="grey",save=".LTA.gene.png")
sc.pl.umap(adata,color="IL21",frameon=False, na_color="grey",save=".IL21.gene.png")
sc.pl.umap(adata,color="IL6",frameon=False, na_color="grey",save=".IL6.gene.png")
sc.pl.umap(adata,color="IL17A",frameon=False, na_color="grey",save=".IL17A.gene.png")
sc.pl.umap(adata,color="FCER2",frameon=False, na_color="grey",save=".FCER2.gene.png")

