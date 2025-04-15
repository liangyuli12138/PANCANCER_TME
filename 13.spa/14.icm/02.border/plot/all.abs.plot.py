import pandas as pd
import matplotlib.pyplot as plt

# 读取CSV文件
df = pd.read_csv('all.immune.dist.csv')

# 删除数据大于4000的数据
#df = df[df['distance'] <= 4000]

# 筛选距离为-1000到1000的数据
filtered_df = df[(df['distance'] >= -1000) & (df['distance'] <= 1000)]

# 统计筛选后的数据的距离分布
distances = filtered_df['distance']
dist_counts_pos = distances[distances > 0].value_counts().sort_index()
dist_counts_neg = distances[distances < 0].value_counts().sort_index()

# 绘制分布图
plt.bar(dist_counts_pos.index, dist_counts_pos.values, color='#4472C4', alpha=0.7, label='Peritumoral')
plt.bar(dist_counts_neg.index, dist_counts_neg.values, color='#ED7D31', alpha=0.7, label='Intratumoral')
plt.xlabel('Distance (um)')
plt.ylabel('Count')
plt.title('Distribution of Distance (filtered)')
plt.legend()
plt.savefig('filtered_dist_plot.abs.png')
plt.close()

# 绘制原始数据的分布图
distances = df['distance']
dist_counts_pos = distances[distances > 0].value_counts().sort_index()
dist_counts_neg = distances[distances < 0].value_counts().sort_index()

# 绘制分布图
plt.bar(dist_counts_pos.index, dist_counts_pos.values, color='#4472C4', alpha=0.7, label='Peritumoral')
plt.bar(dist_counts_neg.index, dist_counts_neg.values, color='#ED7D31', alpha=0.7, label='Intratumoral')
plt.xlabel('Distance (um)')
plt.ylabel('Count')
plt.title('Distribution of Distance (original)')
plt.legend()
plt.savefig('original_dist_plot.abs.png')
plt.close()
