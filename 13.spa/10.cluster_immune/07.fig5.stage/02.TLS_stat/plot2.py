import pandas as pd
import matplotlib.pyplot as plt

# 读取CSV文件，指定表头为第0行
data = pd.read_csv('input.csv', header=0)

# 筛选'TLS_type'列中为'HE_STO'和'STO'的行
filtered_data = data[data['TLS_type'].isin(['HE_STO', 'STO'])]

# 对第一列和第二列进行统计
count = filtered_data.groupby(['groups', 'TLS_type']).size().unstack()

# 输出统计结果到文件
count.to_csv('output.csv')

# 绘制叠加柱状图
count.plot(kind='bar', stacked=True)

# 设置图表标题和坐标轴标签
plt.title('Count by groups and TLS_type')
plt.xlabel('TLS_type')
plt.ylabel('Count')

# 保存图表到文件
plt.savefig('stats.group.type.plot.png')

