import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/07.fig5.stage/06.distance")
data = pd.read_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/07.fig5.stage/06.distance/l2m.stat.merge.csv.sn.stage', index_col=0)

fig, ax = plt.subplots(figsize=(3,4))  # 修改图像大小为(3,4)
ax = sns.boxplot(data, x='stage', y='distance', palette='tab10', linewidth=0.8, order=['early','late'], showfliers=False)
sns.stripplot(data, x='stage', y='distance', jitter=0.2, size=2, alpha=0.8, color='black', order=['early','late'])
plt.axhline(y=0, linestyle='dashed', color='red', zorder=0)   # 添加纵轴上的虚线
#plt.xticks(rotation=30)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('stage', fontsize=10)
ax.set_ylabel('distance', fontsize=10)
ax.yaxis.set_label_coords(-0.25, 0.5)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)

early_data = data[data['stage'] == 'early']['distance']
late_data = data[data['stage'] == 'late']['distance']
t_statistic, p_value = stats.ttest_ind(early_data, late_data)
ax.annotate(f'p = {p_value:.4f}', xy=(0.5, 1.05), xycoords='axes fraction', ha='center', fontsize=10)

plt.tight_layout()  # 调整图像布局
plt.savefig('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/07.fig5.stage/06.distance/stat.distance.png', dpi=300)

