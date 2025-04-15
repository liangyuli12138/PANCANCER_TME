import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/07.fig5.stage/07.GC/new_pvalue/")

data = pd.read_csv('input.csv')
data['log_pvalue'] = -np.log10(data['pvalue'])
group1 = data[data['groups'].isin(['Lymphoid0', 'Lymphoid1', 'Lymphoid2'])]['log_pvalue']
group2 = data[data['groups'] == 'Lymphoid3']['log_pvalue']
t_result = stats.ttest_ind(group1, group2) 
t_text = f"t-test of Lymphoid3 to others, p = {t_result.pvalue:.4f}"

fig, ax = plt.subplots(figsize=(4,5))  # 修改图像大小为(3,4)
ax = sns.boxplot(data, x='groups', y='log_pvalue', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3'], showfliers=False)
sns.stripplot(data, x='groups', y='log_pvalue', jitter=0.2, size=2, alpha=0.8, color='black', order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3'])
plt.xticks(rotation=70)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=12)
ax.set_ylabel('Pvalue of Co-localization', fontsize=12)
ax.yaxis.set_label_coords(-0.25, 0.5)
plt.text(0.5, 1.1, t_text, ha='center', va='center', transform=ax.transAxes, fontsize=8, wrap=True)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
plt.savefig('output.colocalization.LYM.TLS.png', dpi=300)

