import pandas as pd
import matplotlib.pyplot as plt

# 读取文件数据
data = pd.read_csv('merge.immue.csv', delimiter='\t')

# 提取横坐标和纵坐标数据
x = data['Cell_type']
y = data.drop('Cell_type', axis=1)  # 移除'Tissue'列

# 绘制堆积图
plt.figure(figsize=(10, 6))  # 设置图形大小

# 绘制柱状堆积图

bottom = [0] * len(x)  # 底部初始值
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive']  # 颜色列表

for i, col in enumerate(y.columns):
    plt.bar(x, y[col], bottom=bottom, color=colors[i % len(colors)], label=col)
    bottom += y[col]  # 更新底部值

plt.legend(loc='upper right')  # 添加图例
plt.xlabel('Cell_type')  # 设置横坐标标签
plt.ylabel('Cell_Quantity')  # 设置纵坐标标签

plt.xticks(rotation=45)  # 设置横坐标刻度的旋转角度
plt.tight_layout()

plt.savefig("merge.stat.pdf")
