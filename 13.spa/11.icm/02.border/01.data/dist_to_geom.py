from typing import Dict, List, Optional, Tuple, Union, Sequence, Callable, Literal
import numpy as np
import pandas as pd
from pandarallel import pandarallel
from shapely import Polygon, Point, MultiPolygon, LineString
from anndata import AnnData
import scanpy as sc


def set_value(value, default):
    """set value when value is not None else set default"""
    return default if value is None else value

def _parse_coord(
    coords: Union[np.ndarray, pd.Series, pd.DataFrame]
    , is_contour=False):
    """ Parse the coordinates. """
    if isinstance(coords, (pd.Series, pd.DataFrame)):
        coords = coords.values
    else:
        if not hasattr(coords, 'shape'):
            raise ValueError(
                'Coordinates must be a numpy array '
                'or pandas Series/DataFrame')
        
    if is_contour:
        from geopandas import GeoSeries
        
        try:
            coords = GeoSeries(coords.flatten()).values
        except TypeError:
            coords = GeoSeries.from_wkt(coords.flatten()).values
        
        colnames = 'contour'

    else:
        if coords.shape[1] == 1:
            raise ValueError(
                'Coordinates must be 2D or 3D')
            
        colnames = ['x', 'y', 'z'][: coords.shape[1]]
        
    return coords, colnames

def dist_to_geom(
    adata, geom, key_added=None, inplace=True, 
    progress_bar=False, n_workers=32):
    
    points = get_cell_embeddings(adata, reduction='spatial')
    pandarallel.initialize(
        progress_bar=progress_bar, nb_workers=n_workers)   
    
    dist = points.parallel_apply(
        lambda x: geom.distance(Point(x['x'], x['y'])), axis=1)
    key_added = set_value(key_added, 'dist_to_geom')
    
    if inplace:
        adata.obs[key_added] = dist
    else:
        return dist
    
def get_cell_embeddings(
    adata: AnnData, 
    reduction: Optional[str] = 'pca', #TODO set default to None
    contour_col: str = 'contour',  # as_geo
    n_components: Optional[Union[int, Sequence]] = None,
    as_frame: bool = True, 
    ):  
    
    if reduction not in adata.obsm_keys():
        if reduction == 'contour':  
            import geopandas as gpd
            data, _ = _parse_coord(adata.obs[contour_col], True)
            data = gpd.GeoDataFrame(
                {'contour': data}, index=adata.obs.index,
                geometry='contour')
            n_components = None
        else:
            try:
                data = adata.obsm[f'X_{reduction}']
                
            except KeyError:
                raise ValueError(
                    f'{reduction} or X_{reduction} is not in adata.obsm_keys()')
    else:   
        data= adata.obsm[reduction]
    
    data = pd.DataFrame(data, index=adata.obs.index)
    
    N_COMP = data.shape[1]
    if n_components is not None:
        if isinstance(n_components, int):
            data = data.iloc[:, :min(n_components, N_COMP)]
        elif isinstance(n_components, Sequence):
            if N_COMP < len(n_components):
                n_components = slice(None, N_COMP)
            data = data.iloc[:, list(n_components)]
        else:
            raise ValueError(
                'n_component must be int or Sequence[int]')
    else:
        data = data.iloc[:, :min(2, N_COMP)]
    
    if data.shape[1] == 2:
        data.columns = ['x', 'y']   
    elif data.shape[1] == 3:
        data.columns = ['x', 'y', 'z']
    else:
        pass
    
    return data if as_frame else data.values 

if __name__ == '__main__':
    import argparse
    import geopandas as gpd
    from os.path import join
    parser = argparse.ArgumentParser()
    
    reduction = 'spatial'
    
    input_path = parser.add_argument(
        'input_path', type=str, help='input_path')
    chip_id = parser.add_argument(
        'chip_id', type=str, help='chip_id')
    
    adata_path = join(input_path, f'{chip_id}_cellbin.final.h5ad')
    gdf_path = join(input_path, f'{chip_id}_shape.shp')
    
    adata = sc.read(adata_path)
    gdf = gpd.read_file(gdf_path)
    gdf_immue = gdf[gdf['type'].str.contains('malignant')]
    
    def _dist_to_border(df):
        point = Point(df['center_x'], df['center_y'])
        pos_neg = -1 if gdf_immue.contains(point) else 1
        return gdf_immue.explode().exterior.distance(point).min() * pos_neg
    
    dist_to_border = adata.obs.apply(_dist_to_border, axis=1)
    dist_to_border.to_csv(join(input_path, f'{chip_id}_dist_to_border.csv'))