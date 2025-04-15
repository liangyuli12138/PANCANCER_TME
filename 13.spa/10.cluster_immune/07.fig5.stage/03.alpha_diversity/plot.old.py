import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

data = pd.read_csv('l2m.stat.merge.csv.sn', index_col=0)

fig, ax = plt.subplots(figsize=(3,4))  # 修改图像大小为(3,4)
ax = sns.boxplot(data, x='groups', y='Shannon', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'], showfliers=False)
sns.stripplot(data, x='groups', y='Shannon', jitter=0.2, size=2, alpha=0.8, color='black', order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])
plt.xticks(rotation=90)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=8)
ax.set_ylabel('Shannon index of cell types', fontsize=8)
ax.yaxis.set_label_coords(-0.25, 0.5)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
plt.savefig('output.shannon.TLS.png', dpi=300)

