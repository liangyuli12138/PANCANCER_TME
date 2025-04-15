import pandas as pd
import numpy as np

# 读取csv文件
df = pd.read_csv('data.csv')

# 按celltype分组
grouped = df.groupby('celltype')

# 创建一个新的DataFrame保存结果
result = pd.DataFrame(columns=['celltype', 'min5_avg_distance'])

# 遍历每个celltype分组
for celltype, group in grouped:
    # 计算每个cell与其他cell的距离矩阵
    distances = np.zeros((len(group), len(group)))
    for i in range(len(group)):
        for j in range(len(group)):
            distances[i, j] = np.sqrt((group.iloc[i]['x'] - group.iloc[j]['x'])**2 + (group.iloc[i]['y'] - group.iloc[j]['y'])**2)
    
    # 计算每个cell与最近的5个cell的平均距离
    min5_avg_distances = []
    for i in range(len(group)):
        sorted_distances = np.sort(distances[i])
        min5_avg_distances.append(np.mean(sorted_distances[1:6]))  # 排除自身距离
        
    # 将结果添加到新的DataFrame中
    result = result.append({'celltype': celltype, 'min5_avg_distance': min5_avg_distances}, ignore_index=True)

# 输出结果到csv文件
result.to_csv('result.csv', index=False)
