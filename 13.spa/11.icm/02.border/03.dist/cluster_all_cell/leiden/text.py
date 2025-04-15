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

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/11.icm/02.border/03.dist/cluster/leiden")
adata = sc.read_h5ad("immune.cluster.r1.h5ad")

plt.figure(figsize=(50, 50))
sc.pl.umap(adata, color=['region_cluster'], size=300, color_map = 'RdPu', ncols = 2, legend_loc='on data',legend_fontsize=10)

texts = []
for i, (txt) in enumerate(zip(adata.obs_names)):
     texts.append(plt.text(adata.obsm['X_umap'][i, 0], adata.obsm['X_umap'][i, 1], txt, fontsize=5))

adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey', alpha=0))

plt.savefig("leiden/immune.cluster.png",dpi=300, bbox_inches='tight')
plt.savefig("immune.cluster.pdf",dpi=300, bbox_inches='tight')

