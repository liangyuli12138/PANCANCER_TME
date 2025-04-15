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

def _compute_global_shuffled(pos, labels, celltypes, radius):
    # this is global permutation
    labels = labels[np.random.permutation(len(labels))]#[labels[i] for i in np.random.choice(len(labels),len(labels))]
    return _compute_neighborhood(pos, labels, celltypes, radius)


adata = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872C4_cellbin.final.celltype.h5ad"
groupby = "tls_type"
A = sc.read_h5ad(adata)

cellist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/at/D01872C4_Cluster33.input")
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/at/D01872C4_Cluster33.at",index_col=0)
A = A[cellist["cell"],:]
A.obs = A.obs.join(atlist)

celltype_key = groupby
celltypes = list(sorted(A.obs[celltype_key].unique()))
pos = A.obsm['spatial']
labels = A.obs[celltype_key]
radius = 30
neighbors = _compute_neighborhood(pos, labels, celltypes, radius)
neighbors = [[int(element) for element in row] for row in neighbors]

c_list = list(sorted(A.obs[groupby].unique()))
neighbors = pd.DataFrame(neighbors, columns=c_list, index=c_list)
outfile = '/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/out/D01872C4_Cluster33.obs.csv'
neighbors.to_csv(outfile)

n_jobs=1
niter=100
iterations = tqdm(range(niter))
random_freq = Parallel(n_jobs=n_jobs)(delayed(_compute_global_shuffled)(pos, labels, celltypes, radius) for i in iterations)
shuffled_mean = np.dstack(random_freq).mean(2)
shuffled_mean = [[int(element) for element in row] for row in shuffled_mean]
shuffled_mean = pd.DataFrame(shuffled_mean, columns=c_list, index=c_list)
outfile = '/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/out/D01872C4_Cluster33.shuffled.csv'
shuffled_mean.to_csv(outfile)

