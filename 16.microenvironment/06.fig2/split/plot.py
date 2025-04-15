import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 读取CSV文件
data = pd.read_csv('merge.REACTOME.filter.csv')

# 根据variable分组
grouped_data = data.groupby('variable')

# 创建画布和子图
fig, ax = plt.subplots()

# 遍历每个variable，画小提琴图
for name, group in grouped_data:
    # 绘制小提琴图
    sns.violinplot(x='sample', y='score', data=group, ax=ax)
    
# 设置x轴标签旋转角度
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')

# 保存为PDF文件
plt.savefig('violin_plot.pdf', format='pdf')

# 保存为PNG文件
plt.savefig('violin_plot.png', format='png')

# 显示图形
plt.show()
