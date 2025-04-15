import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/09.icm")

data = pd.read_csv('l2m.stat.merge.csv.at.icm', index_col=0)

fig, ax = plt.subplots(figsize=(8,6))  # 修改图像大小为(8,6)
ax = sns.boxplot(data, x='groups', y='ICM1_percent', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'], showfliers=False)
ax = sns.boxplot(data, x='groups', y='ICM1_percent', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])
plt.xticks(rotation=45)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=15)
ax.set_ylabel('ICM1_percent', fontsize=15)
ax.yaxis.set_label_coords(-0.08, 0.5)
groups = data['groups'].unique()
data_per_group = [data[data['groups'] == group]['ICM1_percent'] for group in groups]
kw_result = stats.kruskal(*data_per_group)
kw_text = f"Kruskal-Wallis, p = {kw_result.pvalue:.4f}"
plt.text(0.5, 0.95, kw_text, ha='center', va='top', transform=ax.transAxes, fontsize=12)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
plt.savefig('output.ICM1_percent.score.pdf', dpi=300)

fig, ax = plt.subplots(figsize=(8,6))  # 修改图像大小为(8,6)
ax = sns.boxplot(data, x='groups', y='ICM2_percent', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'], showfliers=False)
ax = sns.boxplot(data, x='groups', y='ICM2_percent', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])
plt.xticks(rotation=45)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=15)
ax.set_ylabel('ICM2_percent', fontsize=15)
ax.yaxis.set_label_coords(-0.08, 0.5)
groups = data['groups'].unique()
data_per_group = [data[data['groups'] == group]['ICM2_percent'] for group in groups]
kw_result = stats.kruskal(*data_per_group)
kw_text = f"Kruskal-Wallis, p = {kw_result.pvalue:.4f}"
plt.text(0.5, 0.95, kw_text, ha='center', va='top', transform=ax.transAxes, fontsize=12)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
plt.savefig('output.ICM2_percent.score.pdf', dpi=300)

fig, ax = plt.subplots(figsize=(8,6))  # 修改图像大小为(8,6)
ax = sns.boxplot(data, x='groups', y='ICM3_percent', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'], showfliers=False)
ax = sns.boxplot(data, x='groups', y='ICM3_percent', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])
plt.xticks(rotation=45)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=15)
ax.set_ylabel('ICM3_percent', fontsize=15)
ax.yaxis.set_label_coords(-0.08, 0.5)
groups = data['groups'].unique()
data_per_group = [data[data['groups'] == group]['ICM3_percent'] for group in groups]
kw_result = stats.kruskal(*data_per_group)
kw_text = f"Kruskal-Wallis, p = {kw_result.pvalue:.4f}"
plt.text(0.5, 0.95, kw_text, ha='center', va='top', transform=ax.transAxes, fontsize=12)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
plt.savefig('output.ICM3_percent.score.pdf', dpi=300)


