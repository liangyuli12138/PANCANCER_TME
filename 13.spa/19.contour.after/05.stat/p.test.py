import pandas as pd
from scipy import stats

# 读取表格数据
data = pd.read_csv('your_file.csv', sep='\t')
# 将每行数据拆分为多列
#data = data[0].str.split('\t', expand=True)
# 提取第一组和第二组的数据
group1 = data.iloc[0:4, 1:]
group2 = data.iloc[4:, 1:]
# 初始化结果列表
results = []
# 遍历每个细胞类型
for cell_type in data.columns[1:]:
    # 提取第一组和第二组的数据
    group1_data = group1[cell_type].tolist()
    group2_data = group2[cell_type].tolist()
    # 进行t检验
    t_statistic, p_value = stats.ttest_ind(group1_data, group2_data)
    # 添加结果到列表
    results.append([cell_type, p_value])
# 创建结果DataFrame
results_df = pd.DataFrame(results, columns=['Cell Type', 'p-value'])
# 输出结果到csv文件
results_df.to_csv('results.csv', index=False)
