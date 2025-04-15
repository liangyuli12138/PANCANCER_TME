import scanpy as sc
import os
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.sparse import csr_matrix
import glob
from anndata import AnnData
from adjustText import adjust_text
from adjustText import adjust_text
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch


adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/05.cluster_by_gene/cluster/immune.cluster.r0.5.h5ad")

plt.figure(figsize=(20, 20))
sc.pl.umap(adata, color=['region_cluster'], size=50, color_map = 'RdPu', ncols = 2, legend_loc='right margin',legend_fontsize=20)

texts = []
for i, txt in enumerate(adata.obs_names):
     texts.append(plt.text(adata.obsm['X_umap'][i, 0], adata.obsm['X_umap'][i, 1], txt, fontsize=1))


#adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey'))
#adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey', lw=0.5))
#adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey'))
adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey', alpha=0))

plt.savefig("immune.cluster.id.r0.5.png",dpi=2000, bbox_inches='tight')
plt.savefig("immune.cluster.id.r0.5.pdf",dpi=2000, bbox_inches='tight')

