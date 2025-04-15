import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/05.ucell/04.ucell.plot")

data = pd.read_csv('merge.merge.all.ucell.list.at.input', index_col=0)

fig, ax = plt.subplots(figsize=(8,6))  # 修改图像大小为(8,6)
ax = sns.boxplot(data, x='groups', y='GC_ucell', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3'], showfliers=False)
sns.stripplot(data, x='groups', y='GC_ucell', jitter=0.2, size=2, alpha=0.8, color='black', order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3'])
plt.xticks(rotation=45)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=15)
ax.set_ylabel('GC score', fontsize=15)
ax.yaxis.set_label_coords(-0.08, 0.5)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
plt.savefig('output.GC.score.png', dpi=300)

fig, ax = plt.subplots(figsize=(8,6))  # 修改图像大小为(8,6)
ax = sns.boxplot(data, x='groups', y='FDC_ucell', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3'], showfliers=False)
sns.stripplot(data, x='groups', y='FDC_ucell', jitter=0.2, size=2, alpha=0.8, color='black', order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3'])
plt.xticks(rotation=45)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=15)
ax.set_ylabel('FDC score', fontsize=15)
ax.yaxis.set_label_coords(-0.08, 0.5)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
plt.savefig('output.FDC.score.png', dpi=300)

