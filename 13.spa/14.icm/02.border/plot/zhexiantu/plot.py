import pandas as pd
import matplotlib.pyplot as plt

# 读取CSV文件
df = pd.read_csv('all.cell.dist.csv')

# 筛选距离为-1000到1000的数据
filtered_df = df[(df['distance'] >= -2000) & (df['distance'] <= 2000)]

# 筛选ICM1, ICM2, ICM3类细胞
icm1_cells = filtered_df[filtered_df['ICM'].str.contains('ICM1')]
icm2_cells = filtered_df[filtered_df['ICM'].str.contains('ICM2')]
icm3_cells = filtered_df[filtered_df['ICM'].str.contains('ICM3')]

# 统计每个区域的细胞数量
def count_cells_in_range(cells, start, end):
    return len(cells[(cells['distance'] >= start) & (cells['distance'] < end)])

ranges = [(-2000,-1750),(-1750,-1500),(-1500,-1250),(-1250,-1000),(-1000, -750), (-750, -500), (-500, -250), (-250, 0), (0, 250), (250, 500), (500, 750), (750, 1000),(1000,1250),(1250,1500),(1500,1750),(1750,2000)]

icm1_counts = [count_cells_in_range(icm1_cells, start, end) for start, end in ranges]
icm2_counts = [count_cells_in_range(icm2_cells, start, end) for start, end in ranges]
icm3_counts = [count_cells_in_range(icm3_cells, start, end) for start, end in ranges]

# 定义计算细胞数量的函数
def count_cells_in_range(cells, start, end):
    return len(cells[(cells['distance'] >= start) & (cells['distance'] <= end)])

# 计算每个区域的ICM1, ICM2, ICM3细胞的百分比
def calculate_percentages(counts):
    total_count = len(window_all_cells)
    percentages = [count / total_count * 100 for count in counts]
    return percentages

icm1_percentages = calculate_percentages(icm1_counts)
icm2_percentages = calculate_percentages(icm2_counts)
icm3_percentages = calculate_percentages(icm3_counts)

# 绘制折线图
x = [f"{start}-{end}" for (start, end) in ranges]
plt.plot(x, icm1_percentages, marker='o', label='ICM1')
plt.plot(x, icm2_percentages, marker='o', label='ICM2')
plt.plot(x, icm3_percentages, marker='o', label='ICM3')

plt.xlabel('Distance Range')
plt.ylabel('Percentage')
plt.title('Percentage of ICM1, ICM2, ICM3 Cells in Different Distance Ranges')
plt.legend()
plt.xticks(rotation=90)

xticks = plt.gca().get_xticks()
xtick_range = xticks[7:9]  # 取第8，9个标签
xtick_mid = sum(xtick_range) / 2  # 计算中点位置
plt.axvline(x=xtick_mid, linestyle='dashed', color='red')

plt.tight_layout()  # 调整图像布局，确保内容完整显示
plt.savefig('icm_percentages.png')
plt.savefig('icm_percentages.pdf')
plt.close()

