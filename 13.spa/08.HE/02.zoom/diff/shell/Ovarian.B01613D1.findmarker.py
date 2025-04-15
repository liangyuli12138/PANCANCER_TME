import os
import datetime
import pandas as pd
import scanpy as sc
from matplotlib.pyplot import rc_context
import sys

def _log(info):
    current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(' '.join([current_time, info]))


if __name__ == '__main__':
    h5ad_f  = sys.argv[1]
    cluster = sys.argv[2]
    outf    = sys.argv[3]

    #wkd = '/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/liaoqijun/project/pancancer/sc/s2.clustering/combined/pancancer' 
    #wkd = '/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/liaoqijun/project/pancancer/sc/test1/test3' 
    _log('Loading data...')
    #h5_file = os.path.join(wkd,  'J-1-C_dubs.h5ad')
    adata = sc.read_h5ad(h5ad_f)
    #adata.uns['log1p']["base"] = None
    cellist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/02.zoom/at/B01613D1.input")
    atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/02.zoom/at/B01613D1.at",index_col=0)
    adata = adata[cellist["cell"],:]
    adata.obs = adata.obs.join(atlist)
    
    '''
    sc.tl.rank_genes_groups(adata, 'groups', method='wilcoxon')
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    marker=pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
        for group in groups for key in ['names', 'pvals']}) 
    markers_file = os.path.join(wkd,  'apan_all_concat_harmony_leiden_markers.csv')
    marker.to_csv(markers_file)
    '''
    _log('Find_all_markers...')
    sc.tl.rank_genes_groups(adata, pts=True,  groupby='region', groups = [cluster], reference='rest', method='t-test')
    ret = sc.get.rank_genes_groups_df(adata, group=None)   #.to_csv('marker_list.csv')
    ret.to_csv(outf)

