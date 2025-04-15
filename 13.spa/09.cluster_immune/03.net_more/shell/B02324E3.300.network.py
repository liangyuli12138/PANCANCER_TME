import scanpy as sc
import numpy as np
import networkx as nx
from scipy.spatial import KDTree
import pandas as pd
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon, Point
import geopandas as gpd
from shapely import affinity
from shapely.geometry import asPolygon
from shapely.geometry import shape
import shapely.wkt as wkt
from shapely.ops import unary_union, transform

# 读取h5ad文件
adata = sc.read_h5ad('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B02324E3_cellbin.final.celltype.h5ad')
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/02.data/at/B02324E3.at",index_col=0)
adata.obs = adata.obs.join(atlist)


# 获取细胞类型、x和y的列名称
celltype_column = 'region'  # 替换为细胞类型列的名称
x_column = 'x'  # 替换为x坐标列的名称
y_column = 'y'  # 替换为y坐标列的名称

# 获取细胞类型、x和y的数据
celltypes = adata.obs[celltype_column].values
x_coords = adata.obs[x_column].values
y_coords = adata.obs[y_column].values

# 构造KDTree来加速邻近点的查找
points = np.column_stack((x_coords, y_coords))
kdtree = KDTree(points)

# 创建一个空的无向图来存储cluster
G = nx.Graph()

# 遍历所有的点
for i, (celltype, x, y) in enumerate(zip(celltypes, x_coords, y_coords)):
    if celltype == 'Immune':
        # 在给定半径内查找相邻点
        indices = kdtree.query_ball_point([x, y], r=50)

        # 过滤掉自身点
        indices = [idx for idx in indices if idx != i]

        # 遍历相邻点，添加到cluster
        for idx in indices:
            neighbor_celltype = celltypes[idx]

            # 如果相邻点具有相同的细胞类型，连接它们
            if neighbor_celltype == celltype:
                G.add_edge(i, idx)

# 获取所有的连通分量（clusters）
clusters = list(nx.connected_components(G))

# 统计每个cluster内的cell数目
cluster_sizes = [len(cluster) for cluster in clusters]

filtered_clusters = [cluster for cluster, size in zip(clusters, cluster_sizes) if size > 100]

clusters = filtered_clusters

cluster_sizes = [len(cluster) for cluster in clusters]

# 获取每个cluster的多边形外边界并向外扩张300
expanded_polygons = []
original_polygons = []
for cluster in clusters:
    # 获取cluster的x和y坐标
    cluster_coords = points[list(cluster)]

    # 计算cluster的凸包
    hull = ConvexHull(cluster_coords)

    # 构造cluster的多边形
    polygon = Polygon(cluster_coords[hull.vertices])

    # 向外扩张100
    # 向外扩张100
    def round_coordinates(x, y, z=None):
        return round(x), round(y)
    # 将原多边形的坐标进行四舍五入
    original_polygon = transform(round_coordinates, unary_union(polygon.buffer(10)))
    # 将原多边形保存为shp文件
    original_polygons.append(original_polygon)
    
    expanded_polygon = transform(round_coordinates, unary_union(polygon.buffer(300)))
    # 将扩展多边形保存为shp文件
    expanded_polygons.append(expanded_polygon)


# 创建包含扩展多边形的GeoDataFrame
expanded_df = gpd.GeoDataFrame(geometry=expanded_polygons)

# 将扩展多边形保存为shapefile文件
expanded_df.to_file('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/03.net_more/out300/B02324E3.expanded_polygons.shp')

# 创建包含原始多边形的GeoDataFrame
original_df = gpd.GeoDataFrame(geometry=original_polygons)

# 将原始多边形保存为shapefile文件
original_df.to_file('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/03.net_more/out300/B02324E3.original_polygons.shp')


# 将结果写入文件
output_file = '/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/03.net_more/out300/B02324E3.output.txt'  # 替换为输出文件的路径和名称

with open(output_file, 'w') as f:
    for i, (cluster, size, original_polygon, expanded_polygon) in enumerate(zip(clusters, cluster_sizes, original_polygons, expanded_polygons)):
        # 获取cluster的obs index或第一列id
        cell_ids = adata.obs.index[list(cluster)]  # 替换为obs的index列或第一列id的名称

        # 获取在扩展范围内的细胞
        expanded_indices = [idx for idx in list(cluster) if expanded_polygon.contains(Point(points[idx]))]
        expanded_cell_ids = adata.obs.index[expanded_indices]

        # 获取原多边形内的全部细胞
        original_indices = [idx for idx in list(cluster) if original_polygon.contains(Point(points[idx]))]
        original_cell_ids = adata.obs.index[original_indices]

        f.write(f"Cluster {i+1} (Size: {size}): {', '.join(cell_ids)}\n")
        f.write(f"Cluster {i+1} Original Boundary: {original_polygon}\n")
        f.write(f"Cluster {i+1} Expanded Boundary: {expanded_polygon}\n")
#       f.write(f"Cluster_ex {i+1} (Size: {len(expanded_indices)}): {', '.join(expanded_cell_ids)}\n")
#       f.write(f"Cluster_ori {i+1} (Size: {len(original_indices)}): {', '.join(original_cell_ids)}\n")


output_file = '/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/03.net_more/out300/B02324E3.output.ex.txt'
with open(output_file, 'w') as f:
    for i, (cluster, size, original_polygon, expanded_polygon) in enumerate(zip(clusters, cluster_sizes, original_polygons, expanded_polygons)):
        # 获取cluster的obs index或第一列id
        cell_ids = []
        for j in range(len(adata)):
            cell = adata.obs.index[j]  # 替换为obs的index列或第一列id的名称
            cell_point = Point(points[j])
            if expanded_polygon.contains(cell_point):
                cell_ids.append(cell)

        # 获取在扩展范围内的细胞
        expanded_indices = [j for j in range(len(adata)) if expanded_polygon.contains(Point(points[j]))]
        expanded_cell_ids = [adata.obs.index[j] for j in expanded_indices]

        # 获取原多边形内的全部细胞
        original_indices = [j for j in range(len(adata)) if original_polygon.contains(Point(points[j]))]
        original_cell_ids = [adata.obs.index[j] for j in original_indices]

        f.write(f"Cluster_ori {i+1} (Size: {len(original_indices)}): {', '.join(original_cell_ids)}\n")
        f.write(f"Cluster_ex {i+1} (Size: {len(expanded_indices)}): {', '.join(expanded_cell_ids)}\n")

