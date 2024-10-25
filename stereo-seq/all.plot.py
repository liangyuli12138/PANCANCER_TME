import pandas as pd
import matplotlib.pyplot as plt

# 读取CSV文件
df = pd.read_csv('all.cell.dist.csv')

# 删除数据大于4000的数据
#df = df[df['distance'] <= 4000]

# 筛选距离为-1000到1000的数据
filtered_df = df[(df['distance'] >= -1000) & (df['distance'] <= 1000)]

# 统计筛选后的数据的距离分布
distances = filtered_df['distance']

icm_cells = filtered_df[filtered_df['ICM'].str.contains('ICM')] 
#icm_counts = pd.cut(icm_cells['distance'], bins=range(-1000, 1001, 50), include_lowest=True).value_counts().sort_index()
#all_counts = pd.cut(distances, bins=range(-1000, 1001, 50), include_lowest=True).value_counts().sort_index()
icm_counts = icm_cells['distance'].value_counts().sort_index()
all_counts = distances.value_counts().sort_index()
#dist_counts_pos = icm_counts / all_counts * 100
dist_counts_pos = icm_counts[all_counts > 100] / all_counts[all_counts > 100] * 100

# 绘制分布图
plt.bar(dist_counts_pos.index[dist_counts_pos.index >= 0], dist_counts_pos.values[dist_counts_pos.index >= 0], color='#4472C4', alpha=0.7, label='Peritumoral') 
plt.bar(dist_counts_pos.index[dist_counts_pos.index < 0], dist_counts_pos.values[dist_counts_pos.index < 0], color='#ED7D31', alpha=0.7, label='Intratumoral')
#plt.bar(dist_counts_pos.head(20).index.astype(str), dist_counts_pos.head(20).values, color='#4472C4', alpha=0.7, label='Peritumoral')
#plt.bar(dist_counts_pos.tail(20).index.astype(str), dist_counts_pos.tail(20).values, color='#ED7D31', alpha=0.7, label='Intratumoral')

plt.xlabel('Distance (um)')
plt.ylabel('Percentage of immune cells')
plt.title('Distribution of Distance ')
plt.legend()
plt.savefig('filtered_dist_plot.percent.pdf')
plt.close()

# 绘制原始数据的分布图
distances = df['distance']

icm_cells = df[df['ICM'].str.contains('ICM')]
icm_counts = icm_cells['distance'].value_counts().sort_index()
all_counts = distances.value_counts().sort_index()
#dist_counts_pos = icm_counts / all_counts * 100
dist_counts_pos = icm_counts[all_counts > 100] / all_counts[all_counts > 100] * 100

# 绘制分布图
plt.bar(dist_counts_pos.index[dist_counts_pos.index >= 0], dist_counts_pos.values[dist_counts_pos.index >= 0], color='#4472C4', alpha=0.7, label='Peritumoral')
plt.bar(dist_counts_pos.index[dist_counts_pos.index < 0], dist_counts_pos.values[dist_counts_pos.index < 0], color='#ED7D31', alpha=0.7, label='Intratumoral')

plt.xlabel('Distance (um)')
plt.ylabel('Percentage of immune cells')
plt.title('Distribution of Distance')
plt.legend()
plt.savefig('all_dist_plot.percent.png')
plt.close()

