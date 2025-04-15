import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/07.fig5.stage/02.TLS_stat")
data = pd.read_csv('l2m.stat.merge.csv.tls.gc', index_col=0)

fig, ax = plt.subplots(figsize=(3,4))  # 修改图像大小为(3,4)
ax = sns.boxplot(data, x='TLS_type', y='Lymphoid_percentage', palette='tab10', linewidth=0.8, order=['HE_STO','STO','Control'], showfliers=False)
sns.stripplot(data, x='TLS_type', y='Lymphoid_percentage', jitter=0.2, size=2, alpha=0.8, color='black', order=['HE_STO','STO','Control'])
plt.xticks(rotation=30)  # 将横轴的标签字体倾斜45度
ax.set_ylabel('Lymphoid_percentage', fontsize=10)
ax.yaxis.set_label_coords(-0.25, 0.5)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
plt.savefig('output.celltype.TLS.png', dpi=300)

fig, ax = plt.subplots(figsize=(3,4))  # 修改图像大小为(3,4)
ax = sns.boxplot(data, x='TLS_type', y='TLS_ucell', palette='tab10', linewidth=0.8, order=['HE_STO','STO','Control'], showfliers=False)
sns.stripplot(data, x='TLS_type', y='TLS_ucell', jitter=0.2, size=2, alpha=0.8, color='black', order=['HE_STO','STO','Control'])
plt.xticks(rotation=30)  # 将横轴的标签字体倾斜45度
ax.set_ylabel('TLS signature', fontsize=10)
ax.yaxis.set_label_coords(-0.25, 0.5)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
plt.savefig('output.ucell.TLS.png', dpi=300)

fig, ax = plt.subplots(figsize=(3,4))  # 修改图像大小为(3,4)
ax = sns.boxplot(data, x='TLS_type', y='GC_ucell', palette='tab10', linewidth=0.8, order=['HE_STO','STO','Control'], showfliers=False)
sns.stripplot(data, x='TLS_type', y='GC_ucell', jitter=0.2, size=2, alpha=0.8, color='black', order=['HE_STO','STO','Control'])
plt.xticks(rotation=30)  # 将横轴的标签字体倾斜45度
ax.set_ylabel('GC signature', fontsize=10)
ax.yaxis.set_label_coords(-0.25, 0.5)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
plt.savefig('output.ucell.GC.png', dpi=300)

