# %%time
from shapely import MultiPolygon
import geopandas as gpd
import pandas as pd
import scanpy as sc
import numpy as np
import os

#chip_ids = ['B01317E6', 'B01613A5', 'B01613B1', 'B01613B2', 'B01613C6',
#       'B01613D1', 'B01615B1', 'B01615B2', 'B01615B3', 'B01615B5',
#       'B02324A1', 'B02324A5', 'B02324B5', 'B02324E3', 'B02324E4',
#       'B02324E5', 'D01872C4', 'D01872C5', 'D01872C6', 'D01872D1',
#       'D01872D2', 'D01872D3', 'D01872D4', 'D01972B6', 'D01972C1',
#       'D01972C4', 'D01972D1', 'D01972D2', 'D01972D6', 'SS200000929BL']

meta = pd.read_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/18.contour/Lymphoid.groups.csv', sep='\t', names=['type', 'group']).set_index('group')

meta['sample'] = meta.index.str.replace('Cluster\d+', '')
chip_ids = meta['sample'].unique()
#bdata = sc.read(f'/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/{chip}_cellbin.final.celltype.h5ad')
all_distances = pd.DataFrame()

for chip in chip_ids:
    filepath = f'/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/03.net/out/{chip}.original_polygons.shp'
    if not os.path.exists(filepath):
        continue
    adata = sc.read(f'/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/{chip}_cellbin.final.celltype.h5ad')
    gdf1 = gpd.read_file(f'/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/03.net/out/{chip}.original_polygons.shp')
    gdf1['group'] = chip + 'Cluster' + (gdf1['FID'] + 1).astype(str)
    inter = np.intersect1d(gdf1.group, meta.index)
    gdf1 = gdf1.loc[gdf1.group.isin(inter)]
    #gdf1['type'] = meta.loc[gdf1.group, 'type'].values
    gdf1['type'] = gdf1['group']
    geo = MultiPolygon(gdf1.geometry.values)
    geos = gdf1.groupby('type').apply(lambda x: MultiPolygon(list(x['geometry']))).to_dict()
    coords = pd.DataFrame(data=adata.obsm['spatial'], index=adata.obs.index, columns=['x', 'y'])
    gdf2 = gpd.points_from_xy(coords['x'], coords['y'])
    # mask = gdf.within(geo)
    # dist = gdf2.distance(geo)
    res = pd.DataFrame(index=coords.index)
    for k, v in geos.items():
        res[k] = gdf2.distance(v)
    res['type'] = res.idxmin(axis=1)
    res['dist'] = res[geos.keys()].min(axis=1)
    conditions = [
        (res['dist'] == 0),
        (res['dist'] <= 100),
        (res['dist'] <= 300),
        (res['dist'] <= 500)
    ]
    choices = ['0', '0-100', '100-300', '300-500']
    res['class'] = np.select(conditions, choices, '>500')
    res.index = res['type'] + "_" + res.index
    res.to_csv(f'./out_file/{chip}_distances.csv', index=True, header=True)
    all_distances = pd.concat([all_distances, res[['type', 'dist', 'class']]])
all_distances.to_csv('all.merge.distances.csv', index=True, header=True)
