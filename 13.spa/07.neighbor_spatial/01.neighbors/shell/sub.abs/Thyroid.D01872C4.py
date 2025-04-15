import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import os
import anndata as ad
from tqdm import tqdm
import scipy.stats
from statsmodels.stats.multitest import multipletests
from sklearn.neighbors import KDTree
import multiprocessing
from joblib import Parallel, delayed
from plotnine import *
import pandas as pd

def count_nearest_neighbors(X,Y,dist_thresh):
    if X.shape[0] > 0 and Y.shape[0] > 0:
        kdtree = KDTree(Y)
        idx, dists = kdtree.query_radius(X, r=dist_thresh, count_only=False, return_distance=True)
        dists = np.hstack(dists)
        return len(dists[dists>0])
    else:
        return 0

def _compute_neighborhood(pos, labels, celltypes, radius):
    neighbors = np.zeros((len(celltypes), len(celltypes)))

    for i, c1 in enumerate(celltypes):
        curr_X = pos[labels==c1]
        #print(c1, curr_X.shape[0])
        for j, c2 in enumerate(celltypes):
            curr_Y = pos[labels==c2]
            if i <= j:
                neighbors[i,j] = np.sum(count_nearest_neighbors(curr_X, curr_Y, dist_thresh=radius))#/curr_X.shape[0]
                neighbors[j,i] = neighbors[i,j]
    return neighbors

adata = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872C4_cellbin.final.celltype.h5ad"
groupby = "celltype"
A = sc.read_h5ad(adata)

celltype_key = groupby
celltypes = list(sorted(A.obs[celltype_key].unique()))
pos = A.obsm['spatial']
labels = A.obs[celltype_key]
radius = 40
neighbors = _compute_neighborhood(pos, labels, celltypes, radius)
neighbors = [[int(element) for element in row] for row in neighbors]
c_list = list(sorted(A.obs[groupby].unique()))
neighbors = pd.DataFrame(neighbors, columns=c_list, index=c_list)
neighbors.to_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/07.neighbor_spatial/01.neighbors/zscore.out.sub.abs/Thyroid.D01872C4.out.sub.abs.csv')
