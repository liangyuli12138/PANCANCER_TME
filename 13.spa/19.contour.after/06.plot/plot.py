import pandas as pd
import matplotlib.pyplot as plt

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/19.contour.after/06.plot")

# 读取csv文件
data = pd.read_csv('pancancer.icar.contour.cell.stat.xls',sep='\t')

# 提取细胞类型的名称
cell_types = data.columns[2:]

# 创建空列表来存储每个细胞类型的数据
all_classes = []
all_lymphoid3_counts = []
all_lymphoid1_2_counts = []

# 循环处理每个细胞类型
for cell_type in cell_types:
    # 提取Lymphoid3和Lymphoid_1_2的数量数据
    lymphoid3_data = data[data['type'] == 'Lymphoid3'][cell_type]
    lymphoid1_2_data = data[data['type'] == 'Lymphoid_1_2'][cell_type]
    # 提取class和数量数据
    classes = lymphoid3_data.index.tolist()
    lymphoid3_counts = lymphoid3_data.tolist()
    lymphoid1_2_counts = lymphoid1_2_data.tolist()
    # 添加到所有细胞类型的数据列表中
    all_classes.append(classes)
    all_lymphoid3_counts.append(lymphoid3_counts)
    all_lymphoid1_2_counts.append(lymphoid1_2_counts)

# 绘制曲线图
plt.figure(figsize=(10, 6))

for i in range(len(cell_types)):
    plt.plot(all_classes[i], all_lymphoid3_counts[i], label='Lymphoid3')
    plt.plot(all_classes[i], all_lymphoid1_2_counts[i], label='Lymphoid_1_2')
    plt.xlabel('class')
    plt.ylabel('count')
    plt.title(cell_types[i])
    plt.legend()
    # 保存为PNG文件
    plt.savefig(cell_types[i] + '.png')
    # 清除当前图形以便绘制下一个图形
    # 设置横坐标刻度和标签
    #plt.xticks([0, 100, 300, 500], ['0', '0-100', '100-300', '300-500'])
    plt.xticks(classes, ['0', '0-100', '100-300', '300-500'])
    plt.clf()
