import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/03.stats/pearson")

# 读取第一个文件
df1 = pd.read_csv('icar.group.celltype.stat.csv', index_col=0)

# 读取第二个文件
df2 = pd.read_csv('icm.group.stat.csv', index_col=0)

import numpy as np

# 计算df1每列与df2每列的相关系数矩阵
correlation_matrix = np.corrcoef(df1.T, df2.T)

# 创建相关系数矩阵的DataFrame
columns = list(df1.columns) + list(df2.columns)
correlation_df = pd.DataFrame(correlation_matrix, columns=columns, index=columns)

# 将相关系数矩阵保存为CSV文件
correlation_df.to_csv('correlation_matrix.csv', index=True)

