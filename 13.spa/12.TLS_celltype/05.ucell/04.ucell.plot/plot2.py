import scanpy as sc
import os
import anndata as ad
import pandas as pd
import numpy as np
adata = sc.read("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/06.stat/immune.cluster.11.r1.5.new.h5ad")
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/05.ucell/04.ucell.plot/merge.all.ucell.list.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata_new = sc.AnnData(obs=adata.obs)
adata_new.obsm['X_umap'] = adata.obsm['X_umap']
adata_concat=adata_new
sc.pl.umap(adata_concat,color=["FDC_ucell"],frameon=False, na_color="grey",save=".Lymphoid_FDC_score.cell.pdf")

