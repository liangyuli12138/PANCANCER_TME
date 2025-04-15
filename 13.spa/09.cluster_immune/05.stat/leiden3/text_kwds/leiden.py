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

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/05.stat/leiden3/text_kwds")
adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/05.stat/leiden3/immune.cluster.r2.h5ad")
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/05.stat/leiden3/text_kwds/list.at",index_col=0)
adata.obs = adata.obs.join(atlist)

plt.figure(figsize=(20, 20))
sc.pl.umap(adata, color=['region_cluster'], size=50, color_map = 'RdPu', ncols = 2, legend_loc='right margin',legend_fontsize=20)

texts = []
for i, (txt, color) in enumerate(zip(adata.obs['TLS'], adata.obs['colour'])):
     texts.append(plt.text(adata.obsm['X_umap'][i, 0], adata.obsm['X_umap'][i, 1], txt, fontsize=1, color=color))

adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey', alpha=0))

plt.savefig("immune.cluster2.TLS.r2.png",dpi=2000, bbox_inches='tight')
plt.savefig("immune.cluster2.TLS.r2.pdf",dpi=2000, bbox_inches='tight')

