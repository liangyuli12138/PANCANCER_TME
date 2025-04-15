import pandas as pd
import matplotlib.pyplot as plt
import os

# 读取CSV文件
df = pd.read_csv('all.cell.dist.csv')

# 筛选距离为-2000到2000的数据
filtered_df = df[(df['distance'] >= -2000) & (df['distance'] <= 2000)]

# 获取celltype列下的所有细胞类型
cell_types = filtered_df['celltype'].unique()

# 设置窗口范围
start_range = -2000
end_range = 2000
window_size = 250
ranges = [(i, i + window_size) for i in range(start_range, end_range, window_size)]

# 统计每个细胞类型在各个窗口下的细胞数量占该窗口总细胞数量的百分比
def calculate_percentages(counts):
    total_count = sum(counts)
    percentages = [count / total_count * 100 for count in counts]
    return percentages

for cell_type in cell_types:
    cell_type_counts = []
    for start, end in ranges:
        cell_count = len(filtered_df[(filtered_df['celltype'] == cell_type) & (filtered_df['distance'] >= start) & (filtered_df['distance'] < end)])
        cell_type_counts.append(cell_count)
    cell_type_percentages = calculate_percentages(cell_type_counts)
    # 绘制折线图
    x = [f"{start}-{end}" for (start, end) in ranges]
    plt.plot(x, cell_type_percentages, marker='o')
    plt.xlabel('Distance Range')
    plt.ylabel('Percentage')
    plt.title(f'Cell Type: {cell_type}')
    plt.xticks(rotation=90)
    plt.tight_layout()
    # 在第8，9个window标签中间画红色虚线
    xticks = plt.gca().get_xticks()
    xtick_range = xticks[7:9]  # 取第8，9个标签
    xtick_mid = sum(xtick_range) / 2  # 计算中点位置
    plt.axvline(x=xtick_mid, linestyle='dashed', color='red')
    # 创建保存目录
    if not os.path.exists('png.out'):
        os.makedirs('png.out')
    # 保存图像为PNG文件
    plt.savefig(f'png.out/{cell_type}_percentages.png')
    plt.close()
