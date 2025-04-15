import scanpy as sc
import numpy as np
import networkx as nx
from scipy.spatial import KDTree
import pandas as pd

# 读取h5ad文件
adata = sc.read_h5ad('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/D01972D1_cellbin.final.celltype.h5ad')
cellist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/01.data/at/D01972D1.input")
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/01.data/at/D01972D1.at",index_col=0)
data = adata[cellist["cell"],:]
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

# 将结果写入文件
output_file = '/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/01.data/output.txt'  # 替换为输出文件的路径和名称

with open(output_file, 'w') as f:
    for i, (cluster, size) in enumerate(zip(clusters, cluster_sizes)):
        # 获取cluster的obs index或第一列id
        cell_ids = adata.obs.index[list(cluster)]  # 替换为obs的index列或第一列id的名称
        f.write(f"Cluster {i+1} (Size: {size}): {', '.join(cell_ids)}\n")
