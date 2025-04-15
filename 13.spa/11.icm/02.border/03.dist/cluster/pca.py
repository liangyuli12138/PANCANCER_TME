import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# 读取CSV文件
data = pd.read_csv('all.immune.dist.stat.csv', index_col=0)

# 创建PCA模型
pca = PCA(n_components=2)

# 对数据进行降维
data_pca = pca.fit_transform(data)

# 绘制PCA平面图
plt.scatter(data_pca[:, 0], data_pca[:, 1])
plt.xlabel('PCA Component 1')
plt.ylabel('PCA Component 2')
plt.title('PCA Plot')

# 在每个点上添加样品名标签
for i, sample in enumerate(data.index):
    plt.text(data_pca[i, 0], data_pca[i, 1], sample, fontsize=3)

plt.tight_layout()

# 保存为PNG文件
plt.savefig('pca_plot.png', dpi=300)
plt.close()

