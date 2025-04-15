import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# 读取csv文件
df = pd.read_csv('your_file.csv', header=None)

# 获取第三列数据
column_3 = df.iloc[:, 2]

# 计算分箱边界
bins = np.arange(min(column_3)-0.0025, max(column_3) + 0.0025, 0.005)

# 绘制分布图
fig, ax = plt.subplots()
ax.hist(column_3, bins=bins, alpha=0.5, color='blue')
ax.set_xlabel('Values')
ax.set_ylabel('Frequency')
ax.set_ylim(0, 50000)
ax.set_xlim(0, 0.2)

# 保存为png文件
plt.savefig('distribution.png')

# 关闭图像窗口
plt.close()
