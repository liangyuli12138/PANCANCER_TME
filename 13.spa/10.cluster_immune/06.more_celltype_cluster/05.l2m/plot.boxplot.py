import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/05.l2m/plot")
data = pd.read_csv('../l2m.stat.merge.csv', index_col=0)

fig, ax = plt.subplots(figsize=(8,6))  # 修改图像大小为(8,6)
ax = sns.boxplot(data, x='groups', y='Lymphoid_percentage', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'], showfliers=False)
sns.stripplot(data, x='groups', y='Lymphoid_percentage', jitter=0.2, size=2, alpha=0.8, color='black', order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])
plt.axhline(y=0.4, linestyle='dashed', color='red', zorder=0)   # 添加纵轴上的虚线
plt.xticks(rotation=45)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=15)
ax.set_ylabel('Lymphoid_percentage', fontsize=15)
ax.yaxis.set_label_coords(-0.08, 0.5)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
# 进行差异度t检验
group1 = data[data['groups'].isin(['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3'])]['Lymphoid_percentage']
group2 = data[data['groups'].isin(['Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])]['Lymphoid_percentage']
t_statistic, p_value = stats.ttest_ind(group1, group2)
# 在图表上添加p值
p_value_text = f"p value of t-test between Lymphoid and Myeloid groups: {p_value:.4f}"
plt.text(0.5, 1.03, p_value_text, transform=ax.transAxes, ha='center', fontsize=12)
plt.savefig('output.Lymphoid_percentage.pdf', dpi=300)

fig, ax = plt.subplots(figsize=(8,6))  # 修改图像大小为(8,6)
ax = sns.boxplot(data, x='groups', y='Myeloid_percentage', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'], showfliers=False)
sns.stripplot(data, x='groups', y='Myeloid_percentage', jitter=0.2, size=2, alpha=0.8, color='black', order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])
plt.axhline(y=0.6, linestyle='dashed', color='red', zorder=0)   # 添加纵轴上的虚线
plt.xticks(rotation=45)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=15)
ax.set_ylabel('Myeloid_percentage', fontsize=15)
ax.yaxis.set_label_coords(-0.08, 0.5)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
# 进行差异度t检验
group1 = data[data['groups'].isin(['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3'])]['Myeloid_percentage']
group2 = data[data['groups'].isin(['Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])]['Myeloid_percentage']
t_statistic, p_value = stats.ttest_ind(group1, group2)
# 在图表上添加p值
p_value_text = f"p value of t-test between Lymphoid and Myeloid groups: {p_value:.4f}"
plt.text(0.5, 1.03, p_value_text, transform=ax.transAxes, ha='center', fontsize=12)
plt.savefig('output.Myeloid_percentage.pdf', dpi=300)

fig, ax = plt.subplots(figsize=(8,6))  # 修改图像大小为(8,6)
ax = sns.boxplot(data, x='groups', y='log10_Lym2Mye_ratio', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'], showfliers=False)
sns.stripplot(data, x='groups', y='log10_Lym2Mye_ratio', jitter=0.2, size=2, alpha=0.8, color='black', order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])
plt.axhline(y=0, linestyle='dashed', color='red', zorder=0)   # 添加纵轴上的虚线
plt.xticks(rotation=45)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=15)
ax.set_ylabel('log10_Lym2Mye_ratio', fontsize=15)
ax.yaxis.set_label_coords(-0.08, 0.5)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
# 进行差异度t检验
group1 = data[data['groups'].isin(['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3'])]['log10_Lym2Mye_ratio']
group2 = data[data['groups'].isin(['Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])]['log10_Lym2Mye_ratio']
t_statistic, p_value = stats.ttest_ind(group1, group2)
# 在图表上添加p值
p_value_text = f"p value of t-test between Lymphoid and Myeloid groups: {p_value:.4f}"
plt.text(0.5, 1.03, p_value_text, transform=ax.transAxes, ha='center', fontsize=12)
plt.savefig('output.log10_Lym2Mye_ratio.pdf', dpi=300)

fig, ax = plt.subplots(figsize=(8,6))  # 修改图像大小为(8,6)
ax = sns.boxplot(data, x='groups', y='density', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'], showfliers=False)
sns.stripplot(data, x='groups', y='density', jitter=0.2, size=2, alpha=0.8, color='black', order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])
#plt.axhline(y=0.4, linestyle='dashed', color='red', zorder=0)   # 添加纵轴上的虚线
plt.xticks(rotation=45)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=15)
ax.set_ylabel('density', fontsize=15)
ax.yaxis.set_label_coords(-0.08, 0.5)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
# 进行差异度t检验
group1 = data[data['groups'].isin(['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3'])]['density']
group2 = data[data['groups'].isin(['Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])]['density']
t_statistic, p_value = stats.ttest_ind(group1, group2)
# 在图表上添加p值
p_value_text = f"p value of t-test between Lymphoid and Myeloid groups: {p_value:.4f}"
plt.text(0.5, 1.03, p_value_text, transform=ax.transAxes, ha='center', fontsize=12)
plt.savefig('output.density.pdf', dpi=300)

fig, ax = plt.subplots(figsize=(8,6))  # 修改图像大小为(8,6)
ax = sns.boxplot(data, x='groups', y='area', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'], showfliers=False)
sns.stripplot(data, x='groups', y='area', jitter=0.2, size=2, alpha=0.8, color='black', order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])
#plt.axhline(y=0.4, linestyle='dashed', color='red', zorder=0)   # 添加纵轴上的虚线
plt.xticks(rotation=45)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=15)
ax.set_ylabel('area', fontsize=15)
plt.ylim(0, 0.8)
ax.yaxis.set_label_coords(-0.08, 0.5)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
# 进行差异度t检验
group1 = data[data['groups'].isin(['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3'])]['area']
group2 = data[data['groups'].isin(['Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])]['area']
t_statistic, p_value = stats.ttest_ind(group1, group2)
# 在图表上添加p值
p_value_text = f"p value of t-test between Lymphoid and Myeloid groups: {p_value:.4f}"
plt.text(0.5, 1.03, p_value_text, transform=ax.transAxes, ha='center', fontsize=12)
plt.savefig('output.area.pdf', dpi=300)

fig, ax = plt.subplots(figsize=(8,6))  # 修改图像大小为(8,6)
ax = sns.boxplot(data, x='groups', y='elongation', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'], showfliers=False)
sns.stripplot(data, x='groups', y='elongation', jitter=0.2, size=2, alpha=0.8, color='black', order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])
#plt.axhline(y=0.4, linestyle='dashed', color='red', zorder=0)   # 添加纵轴上的虚线
plt.xticks(rotation=45)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=15)
ax.set_ylabel('elongation', fontsize=15)
ax.yaxis.set_label_coords(-0.08, 0.5)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
# 进行差异度t检验
group1 = data[data['groups'].isin(['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3'])]['elongation']
group2 = data[data['groups'].isin(['Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])]['elongation']
t_statistic, p_value = stats.ttest_ind(group1, group2)
# 在图表上添加p值
p_value_text = f"p value of t-test between Lymphoid and Myeloid groups: {p_value:.4f}"
plt.text(0.5, 1.03, p_value_text, transform=ax.transAxes, ha='center', fontsize=12)
plt.savefig('output.elongation.pdf', dpi=300)

fig, ax = plt.subplots(figsize=(8,6))  # 修改图像大小为(8,6)
ax = sns.boxplot(data, x='groups', y='distance', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'], showfliers=False)
sns.stripplot(data, x='groups', y='distance', jitter=0.2, size=2, alpha=0.8, color='black', order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])
#plt.axhline(y=0.4, linestyle='dashed', color='red', zorder=0)   # 添加纵轴上的虚线
plt.xticks(rotation=45)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=15)
ax.set_ylabel('distance', fontsize=15)
ax.yaxis.set_label_coords(-0.08, 0.5)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
# 进行差异度t检验
group1 = data[data['groups'].isin(['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3'])]['distance']
group2 = data[data['groups'].isin(['Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])]['distance']
t_statistic, p_value = stats.ttest_ind(group1, group2)
# 在图表上添加p值
p_value_text = f"p value of t-test between Lymphoid and Myeloid groups: {p_value:.4f}"
plt.text(0.5, 1.03, p_value_text, transform=ax.transAxes, ha='center', fontsize=12)
plt.savefig('output.distance.pdf', dpi=300)

