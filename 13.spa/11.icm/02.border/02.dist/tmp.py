import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances

# 读取csv文件
df = pd.read_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/11.icm/01.plot/at/aaaa.all.at')

# 拆分id列获取x和y坐标
df['x'] = df['id'].apply(lambda x: int(x.split('_')[0]))
df['y'] = df['id'].apply(lambda x: int(x.split('_')[1]))

# 计算每个细胞距离最近10个other细胞的平均距离
df['dist_other'] = df.apply(lambda row: round(euclidean_distances(row[['x','y']].values.reshape(1, -1), df[df['celltype'] == 'other'][['x','y']])[0].mean()), axis=1)

# 计算每个细胞距离最近10个malignant细胞的平均距离
df['dist_malignant'] = df.apply(lambda row: round(euclidean_distances(row[['x','y']].values.reshape(1, -1), df[df['celltype'] == 'malignant'][['x','y']])[0].mean()), axis=1)

# 保存结果为新的csv文件
df[['id', 'celltype', 'dist_other', 'dist_malignant']].to_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/11.icm/02.border/02.dist/out/aaaa.dist.csv', index=False)
