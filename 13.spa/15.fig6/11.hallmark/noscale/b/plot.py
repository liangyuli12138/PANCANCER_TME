import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import pandas as pd
import os


# 读取文件数据
data = pd.read_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/11.hallmark/noscale/merge.gene.signature.score.filter.csv')

# 获取celltype的列表
cell_types = data.columns[1:]

# 分组标签
groups = ['Lymphoid1', 'Lymphoid2', 'Lymphoid3']
group_colors = ['#306F9F', '#306F9F', '#DF7D26']

# 遍历每个celltype
for cell_type in cell_types:
    # 创建箱线图
    fig, ax = plt.subplots(figsize=(4, 4))
    # 存储每个分组的数值列表
    values = []
    # 根据分组获取对应的数值
    for group in groups:
        group_values = data[data['groups'] == group][cell_type]
        values.append(group_values.tolist())
    # 绘制箱线图
    bp = ax.boxplot(values, labels=groups, showfliers=False, patch_artist=True)
    # 设置箱线图颜色
    for patch, color in zip(bp['boxes'], group_colors):
        patch.set_facecolor(color)
    # 进行t检验
    ttest_result = ttest_ind(values[1], values[2])
    p_value = ttest_result.pvalue
    # 在图的正上方添加t检验结果
    ax.set_title(f'Boxplot for {cell_type}\nt-test p-value: {p_value:.4f}')
    ax.set_xlabel('Groups')
    ax.set_ylabel('Gene signatuer score')
    # 保存图表为pdf和png文件
    plt.savefig(f'boxplot_{cell_type}.b.png')
    plt.close('all')
