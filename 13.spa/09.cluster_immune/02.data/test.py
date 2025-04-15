import scanpy as sc
import numpy as np
import networkx as nx
from scipy.spatial import KDTree
import pandas as pd
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon, Point
adata = sc.read_h5ad('celltype.h5ad')
celltype_column = 'region' 
x_column = 'x' 
y_column = 'y' 

celltypes = adata.obs[celltype_column].values
x_coords = adata.obs[x_column].values
y_coords = adata.obs[y_column].values

points = np.column_stack((x_coords, y_coords))
kdtree = KDTree(points)

G = nx.Graph()

for i, (celltype, x, y) in enumerate(zip(celltypes, x_coords, y_coords)):
    if celltype == 'Immune':
        indices = kdtree.query_ball_point([x, y], r=50)
        indices = [idx for idx in indices if idx != i]
        for idx in indices:
            neighbor_celltype = celltypes[idx]
            if neighbor_celltype == celltype:
                G.add_edge(i, idx)
clusters = list(nx.connected_components(G))

cluster_sizes = [len(cluster) for cluster in clusters]
filtered_clusters = [cluster for cluster, size in zip(clusters, cluster_sizes) if size > 100]
clusters = filtered_clusters
cluster_sizes = [len(cluster) for cluster in clusters]

expanded_polygons = []
original_polygons = []
for cluster in clusters:
    cluster_coords = points[list(cluster)]
    hull = ConvexHull(cluster_coords)
    polygon = Polygon(cluster_coords[hull.vertices])
    expanded_polygon = polygon.buffer(200)
    expanded_polygons.append(expanded_polygon)
    original_polygons.append(polygon)

output_file = 'output.txt'

with open(output_file, 'w') as f:
    for i, (cluster, size, original_polygon, expanded_polygon) in enumerate(zip(clusters, cluster_sizes, original_polygons, expanded_polygons)):
        cell_ids = adata.obs.index[list(cluster)]  # 替换为obs的index列或第一列id的名称
        expanded_indices = [idx for idx in list(cluster) if expanded_polygon.contains(Point(points[idx]))]
        expanded_cell_ids = adata.obs.index[expanded_indices]
        original_indices = [idx for idx in list(cluster) if original_polygon.contains(Point(points[idx]))]
        original_cell_ids = adata.obs.index[original_indices]
        f.write(f"Cluster {i+1} (Size: {size}): {', '.join(cell_ids)}\n")
        f.write(f"Cluster {i+1} Original Boundary: {original_polygon}\n")
        f.write(f"Cluster {i+1} Expanded Boundary: {expanded_polygon}\n")
        f.write(f"Cluster_ex {i+1} (Size: {len(expanded_indices)}): {', '.join(expanded_cell_ids)}\n")
        f.write(f"Cluster_ori {i+1} (Size: {len(original_indices)}): {', '.join(original_cell_ids)}\n")
