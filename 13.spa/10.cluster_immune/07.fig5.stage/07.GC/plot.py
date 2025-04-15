import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns


df = pd.read_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/07.fig5.stage/07.GC/input.csv')

# 创建不同stage对应的颜色字典
color_dict = {'early': 'red', 'late': 'blue'}

# 根据stage绘制散点图
for index, row in df.iterrows():
    if index == 0:
        continue  # 跳过表头行
    stage = row['stage']
    zscore = row['Zscore']
    pvalue = row['pvalue']
    # 对pvalue取-log10处理
    pvalue = -np.log10(pvalue)
    plt.scatter(pvalue, zscore, color=color_dict[stage])

# 设置标题和标签
plt.title('GC cluster')
plt.xlabel('-log10(pvalue)')
plt.ylabel('Zscore')

plt.savefig('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/07.fig5.stage/07.GC/gc.cluster.png', dpi=300)

data=df
early_data = data[data['stage'] == 'early']['log_pvalue']
late_data = data[data['stage'] == 'late']['log_pvalue']

t_statistic, p_value = stats.ttest_ind(early_data, late_data)

data['log_pvalue'] = -np.log10(data['pvalue'])
fig, ax = plt.subplots(figsize=(3,4))  # 修改图像大小为(3,4)
ax = sns.boxplot(data, x='stage', y='log_pvalue', palette='tab10', linewidth=0.8, order=['early','late'], showfliers=False)
sns.stripplot(data, x='stage', y='log_pvalue', jitter=0.2, size=2, alpha=0.8, color='black', order=['early','late'])
#plt.axhline(y=0, linestyle='dashed', color='red', zorder=0)   # 添加纵轴上的虚线
#plt.xticks(rotation=30)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('stage', fontsize=10)
ax.set_ylabel('-log10(pvalue)', fontsize=10)
ax.yaxis.set_label_coords(-0.25, 0.5)
plt.subplots_adjust(top=0.9)
ax.set_ylim([0, 110]) 
plt.text(0.5, 100, f"p-value: {p_value:.4f}", ha='center', va='center', fontsize=8)

plt.tight_layout()  # 调整图像布局
plt.savefig('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/07.fig5.stage/07.GC/gc.pvalue.png', dpi=300)
plt.close()
#GC_ucell
data=pd.read_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/07.fig5.stage/07.GC/input2.csv') 
early_data = data[data['stage'] == 'early']['GC_ucell']
late_data = data[data['stage'] == 'late']['GC_ucell']
t_statistic, p_value = stats.ttest_ind(early_data, late_data)

fig, ax = plt.subplots(figsize=(3,4))  # 修改图像大小为(3,4)
ax = sns.boxplot(data, x='stage', y='GC_ucell', palette='tab10', linewidth=0.8, order=['early','late'], showfliers=False)
sns.stripplot(data, x='stage', y='GC_ucell', jitter=0.2, size=2, alpha=0.8, color='black', order=['early','late'])
#plt.axhline(y=0, linestyle='dashed', color='red', zorder=0)   # 添加纵轴上的虚线
#plt.xticks(rotation=30)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('stage', fontsize=10)
ax.set_ylabel('GC_ucell', fontsize=10)
#ax.yaxis.set_label_coords(-0.25, 0.5)
#plt.subplots_adjust(top=0.9)
plt.text(0.5, 80, f"p-value: {p_value:.4f}", ha='center', va='center', fontsize=8)

plt.tight_layout()  # 调整图像布局
plt.savefig('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/07.fig5.stage/07.GC/gc.ucell.png', dpi=300)
plt.close()


