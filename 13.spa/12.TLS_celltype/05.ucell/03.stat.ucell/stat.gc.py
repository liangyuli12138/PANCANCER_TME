import pandas as pd
import matplotlib.pyplot as plt
import os

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/05.ucell/03.stat.ucell")

# 读取at.list文件路径
with open('at.list', 'r') as file:
    file_paths = [line.strip() for line in file]

# 存储所有样品的GC_UCell数值
gc_ucell_values = []

# 处理每个样品的csv文件
for file_path in file_paths:
    # 提取样品id
    sample_id = file_path.split('/')[-1].split('.')[0]
    # 读取csv文件
    df = pd.read_csv(file_path, index_col=0)
    # 筛选LM列为Lymphoid的细胞
    lymphoid_cells = df[(df['LM'] == 'Lymphoid') & df['celltype'].str.contains('Lymphoid_B')]
    # 保存筛选结果到文件
    output_file_path = f'./lym.gc.at/{sample_id}_lymphoid.csv'
    lymphoid_cells.to_csv(output_file_path, index=True)
    # 构建UCell.score.csv文件路径
    ucell_file_path = f'/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/05.ucell/03.stat.ucell/gc.out/{sample_id}.UCell.score.csv'
    # 读取UCell.score.csv文件
    ucell_df = pd.read_csv(ucell_file_path, index_col=0)
    # 判断lymphoid_cells.index是否在ucell_df.index中
    valid_indices = lymphoid_cells.index.isin(ucell_df.index)
    # 筛选上一步每个样品筛选的cell
    filtered_cells = ucell_df.loc[lymphoid_cells.index[valid_indices]]    
    # 提取GC_UCell列数值
    gc_ucell_values.extend(filtered_cells['GC_UCell'])
    # 绘制分布图
    plt.hist(filtered_cells['GC_UCell'], bins=20)
    plt.title(f'Distribution of GC_UCell for {sample_id}')
    plt.xlabel('GC_UCell')
    plt.ylabel('Count')
    plt.savefig(f'./stat.gc.png/{sample_id}.png')
    plt.close()

# 绘制所有样品的GC_UCell数值分布图
plt.hist(gc_ucell_values, bins=20)
plt.title('Distribution of GC_UCell for All Samples')
plt.xlabel('GC_UCell')
plt.ylabel('Count')
plt.ylim(0, 1500)  # 设置纵坐标的范围
plt.savefig('./stat.gc.png/all_samples.png')
plt.close()
