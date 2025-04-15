import pandas as pd
import matplotlib.pyplot as plt

# 聚类文件路径
cluster_file_path = "cluster.csv"

# 读取聚类文件
cluster_df = pd.read_csv(cluster_file_path)

# 创建一个空的字典，用于存储每个cluster的GC_UCell列数值
cluster_gc_ucell_values = {}

for file_path in file_paths:
    # 提取样品id
    sample_id = file_path.split('/')[-1].split('.')[0]
    # 读取csv文件
    df = pd.read_csv(file_path, index_col=0)
    # 筛选LM列为Lymphoid的细胞
    lymphoid_cells = df[(df['LM'] == 'Lymphoid') & df['celltype'].str.contains('Lymphoid_B')]
    
    # 匹配每个cell的聚类
    lymphoid_cells_with_cluster = lymphoid_cells.merge(cluster_df, left_index=True, right_on='cell_id')
    
    # 保存筛选结果到文件
    output_file_path = f'./lym.gc.at/{sample_id}_lymphoid.csv'
    lymphoid_cells_with_cluster.to_csv(output_file_path, index=True)
    
    # 构建UCell.score.csv文件路径
    ucell_file_path = f'/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/05.ucell/03.stat.ucell/gc.out/{sample_id}.UCell.score.csv'
    # 读取UCell.score.csv文件
    ucell_df = pd.read_csv(ucell_file_path, index_col=0)
    
    # 判断lymphoid_cells.index是否在ucell_df.index中
    valid_indices = lymphoid_cells_with_cluster.index.isin(ucell_df.index)
    
    # 筛选上一步每个样品筛选的cell之后，匹配上一步的cluster聚类，根据聚类提取GC_UCell列数值
    filtered_cells = ucell_df.loc[lymphoid_cells_with_cluster.index[valid_indices]]
    
    # 提取每个cluster的GC_UCell列数值
    for cluster_id, cluster_cells in filtered_cells.groupby('cluster'):
        cluster_gc_ucell_values.setdefault(cluster_id, []).extend(cluster_cells['GC_UCell'])
    
    # 绘制每个cluster的GC_UCell列数值分布图
    for cluster_id, gc_ucell_values in cluster_gc_ucell_values.items():
        plt.hist(gc_ucell_values, bins=500)
        plt.title(f'Distribution of GC_UCell for {sample_id}, Cluster {cluster_id}')
        plt.xlabel('GC_UCell')
        plt.ylabel('Count')
        plt.savefig(f'./stat.gc.cluster.png/{sample_id}_cluster{cluster_id}.png')
        plt.close()
