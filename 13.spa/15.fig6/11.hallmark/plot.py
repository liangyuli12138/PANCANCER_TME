import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import pandas as pd

# 读取文件数据
data = pd.read_csv('merge.gene.signature.score.filter.csv')

# 获取celltype的列表
cell_types = data.columns[1:]

# 分组标签
groups = ['Lymphoid1', 'Lymphoid2', 'Lymphoid3']
group_colors = ['blue', 'blue', 'orange']

# 遍历每个celltype
for cell_type in cell_types:
    # 创建箱线图
    fig, ax = plt.subplots()
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
    ttest_result = ttest_ind(values[0], values[1])
    p_value = ttest_result.pvalue
    # 在图的正上方添加t检验结果
    ax.set_title(f'Boxplot for {cell_type}\nt-test p-value: {p_value:.4f}')
    ax.set_xlabel('Groups')
    ax.set_ylabel('Cell Type Proportion')
    # 保存图表为pdf和png文件
    plt.savefig(f'boxplot_{cell_type}.png')
