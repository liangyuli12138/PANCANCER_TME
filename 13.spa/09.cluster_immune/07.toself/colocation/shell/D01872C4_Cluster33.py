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

def calc_pval_onesided(obs, rand):
    z = (obs-np.mean(rand))/np.std(rand)
    if z > 0:
        return scipy.stats.norm.sf(abs(z))
    else:
        return 1
    
    
def calc_pval(obs, rand, empirical=False):
    if empirical:
        return np.sum(obs <= np.array(rand))/len(rand)
    else:
        z = (obs - np.mean(rand))/np.std(rand)
        return scipy.stats.norm.sf(abs(z))*2


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
   
def _comput_local_shuffled(pos, labels, celltypes, radius, local_permute_radius=300):
    # this is local permutation
    kdtree = KDTree(pos)
    neighbors = kdtree.query_radius(pos, r=local_permute_radius, count_only=False, return_distance=False)
    
    pre_order = np.arange(len(labels))
    remaining = set(np.arange(len(labels)))
    
    # extract neighbors
    for i in range(len(neighbors)):
        if i in remaining:
            j = np.random.choice(neighbors[i])
            pre_order[i] = j
            pre_order[j] = i
            
            remaining.discard(i)
            remaining.discard(j)
            
    labels = labels[pre_order]
    return _compute_neighborhood(pos, labels, celltypes, radius)

def compute_celltype_neighborhood(A, celltype_key, celltypes=None, radius=40, 
                                  niter=1000, permutation_method='global', oneside=False,
                                  local_permute_radius=300, n_jobs=-1):
    from functools import partial
    assert permutation_method in ['global', 'local']
    
    if celltypes is None:
        celltypes = list(sorted(A.obs[celltype_key].unique()))
    pos = A.obsm['spatial']
    labels = A.obs[celltype_key]
    neighbors = _compute_neighborhood(pos, labels, celltypes, radius)
    iterations = tqdm(range(niter))
    # for each iteration, shuffle celltype labels
    # num_cores = multiprocessing.cpu_count()
    compute_shuff = (_compute_global_shuffled 
                     if permutation_method=='global' 
                     else partial(_comput_local_shuffled, 
                                  local_permute_radius=local_permute_radius)
                     )
    random_freq = Parallel(n_jobs=n_jobs)(delayed(compute_shuff)(pos, labels, celltypes, radius) for i in iterations)    
    #print(len(random_freq))
    # z score
    zs = np.zeros_like(neighbors)
    pval = np.zeros_like(neighbors)

    shuffled_mean = np.dstack(random_freq).mean(2)
    shuffled_std = np.std(np.dstack(random_freq),2)
    for i in range(neighbors.shape[0]):
        for j in range(neighbors.shape[1]):
            zs[i,j] = (neighbors[i,j] - shuffled_mean[i,j])/shuffled_std[i,j]
            if not oneside:
                pval[i,j] = calc_pval(neighbors[i,j],  np.dstack(random_freq)[i,j,:])#np.sum(neighbors[i,j] <= np.dstack(random_freq)[i,j,:])/niter#np.sum(neighbors[i,j] <= np.dstack(random_freq)[i,j,:])/niter #calc_pval(neighbors[i,j],  np.dstack(random_freq)[i,j,:])#np.sum(neighbors[i,j] <= np.dstack(random_freq)[i,j,:])/niter
            else:
                pval[i,j] = calc_pval_onesided(neighbors[i,j], np.dstack(random_freq)[i,j,:])
                
    for i in range(pval.shape[0]):
        pval[i,:] = multipletests(pval[i,:], method='fdr_bh')[1]
        
    pval = pval.reshape(neighbors.shape)
    return neighbors, shuffled_mean, zs, pval

def _set_align(rotation):
    if rotation in [0, 90, -90]:
        return {'ha': 'center', 'va': 'top'}

    else:
        return {'ha': 'right', 'va': 'top'}

def dot_plot(data, x, y, size, color, size_title, colorbar_title,
              cmap='Purples', scale_range=(0, 6), fontsize=11,
              legend_position='right', rotation=0, grid_line=False, 
              linetype='dashed', grid_color='0.8', ticks=False,
              fontstyle='italic', figsize=(8, 8), **kwargs):
    """dot plot"""

    dot_theme = (theme_classic() +
        theme(
            legend_position=legend_position, 
            text=element_text(color='black'),
            legend_title=element_text(size=fontsize), 
            legend_text=element_text(size=fontsize),
            axis_title_x=element_blank(),
            axis_text_x=element_text(rotation=rotation, 
                                     size=fontsize, 
                                     fontstyle=fontstyle,
                                     **_set_align(rotation)
                                     ),
            axis_text_y=element_text(size=fontsize),
            axis_title_y=element_blank(),
            axis_line_x=element_line(color="black", size=1.2),
            panel_grid_major=(element_line(color=grid_color, 
                                           linetype=linetype)
                              if grid_line else 
                              element_blank()),
            panel_ontop=False,
            figure_size=figsize, 
            )
        )
    
    p = (
        ggplot() + 
        dot_theme +
        geom_point(data, 
                    aes(x=x, y=y, size=size, color=color),
                    ) + 
        xlab('') + labs(color=colorbar_title, 
                        size=size_title) + 
        scale_colour_cmap(
            cmap_name=cmap, 
            guide=guide_colorbar(
                ticks=False, order=1, raster=True,
                barwidth=11, barheight=8.5),
            **kwargs) + 
        scale_size(range=scale_range)
    )
    return p

adata = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01872C4_cellbin.final.celltype.h5ad"
adata = sc.read_h5ad(adata)
cellist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/at/D01872C4_Cluster33.input")
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/at/D01872C4_Cluster33.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)

groupby = "tls_type"

neighbors, shuffled_mean, zscore, nbor_pvals = compute_celltype_neighborhood(adata, groupby, niter=500, radius=30, permutation_method='global', oneside=False, n_jobs=2)

logfc = np.log2(neighbors / shuffled_mean)
# nbor_pvals = nbor_pvals[1].reshape(neighbors.shape)
c_list = list(sorted(adata.obs[groupby].unique()))
df = pd.DataFrame(nbor_pvals, columns=c_list, index=c_list)
df_to_merge = pd.DataFrame(zscore, columns=c_list, index=c_list)

df.to_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/out/D01872C4_Cluster33.out.global.nbor_pvals.csv')
df_to_merge.to_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/out/D01872C4_Cluster33.out.global.zscore.csv')

