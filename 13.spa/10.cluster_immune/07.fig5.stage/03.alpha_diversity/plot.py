import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/07.fig5.stage/03.alpha_diversity")

data = pd.read_csv('l2m.stat.merge.csv.sn.filter', index_col=0)
group1 = data[data['groups'].isin(['Lymphoid1'])]['Shannon']
group2 = data[data['groups'] == 'Lymphoid2']['Shannon']
group3 = data[data['groups'] == 'Lymphoid3']['Shannon']

#t_result = stats.ttest_ind(group1, group2) 
statistic, p_value = stats.kruskal(group1,group2,group3)
t_text = f"Kruskal-Wallis, p = {p_value:.4f}"

fig, ax = plt.subplots(figsize=(3,5))  # 修改图像大小为(3,4)
ax = sns.boxplot(data, x='groups', y='Shannon', palette='tab10', linewidth=0.8, order=['Lymphoid1','Lymphoid2','Lymphoid3'], showfliers=False)
sns.stripplot(data, x='groups', y='Shannon', jitter=0.2, size=2, alpha=0.8, color='black', order=['Lymphoid1','Lymphoid2','Lymphoid3'])
plt.xticks(rotation=70)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=12)
ax.set_ylabel('Shannon index of cell types', fontsize=12)
ax.yaxis.set_label_coords(-0.25, 0.5)
plt.text(0.5, 1.1, t_text, ha='center', va='center', transform=ax.transAxes, fontsize=8, wrap=True)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
plt.savefig('output.shannon.LYM.TLS.pdf', dpi=300)

