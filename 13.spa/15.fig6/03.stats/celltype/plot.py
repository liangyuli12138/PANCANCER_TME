import pandas as pd
import matplotlib.pyplot as plt

# 读取文件数据
data = pd.read_csv('icar.group.celltype.stat.csv')

# 获取celltype的列表
cell_types = data.columns[1:]

# 去重得到分组标签
groups = data['groups'].unique()

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
#    ax = sns.violinplot(data=values, inner='quartile', palette='Set3',  showmedians=True, showextrema=False)
    ax.boxplot(values, labels=groups, showfliers=False)
#    sns.stripplot(values, jitter=0.2, size=3, alpha=0.8, color='black',outlier_prop=0.95)
    # 进行Kruskal-Wallis检验
    kruskal_result = kruskal(*values)
    p_value = kruskal_result.pvalue
    # 在图的正上方添加Kruskal-Wallis检验结果
    ax.set_title(f'Boxplot for {cell_type}\nKruskal-Wallis p-value: {p_value:.4f}')
    ax.set_xlabel('Groups')
    ax.set_ylabel('Cell Type Proportion')
    # 保存图表为pdf和png文件
    plt.savefig(f'boxplot_{cell_type}.png')

