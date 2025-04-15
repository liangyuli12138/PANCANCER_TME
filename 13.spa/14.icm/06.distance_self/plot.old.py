import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('Lymphoid123.distance.stat.csv')


# 提取Lymphoid列的唯一值
lymphoid_values = ["Lymphoid1", "Lymphoid2", "Lymphoid3"]

# 创建整合图
fig, ax = plt.subplots()

# 遍历每个lymphoid值
for lymphoid in lymphoid_values:
    # 筛选相应的lymphoid列数据
    lymphoid_data = data[data['lym'] == lymphoid]
    # 遍历每个type1值
    for type1 in lymphoid_data['type1'].unique():
        # 提取对应type1和lymphoid的dis平均值数据
        dis_avg = lymphoid_data[lymphoid_data['type1'] == type1]['dis'].mean()
        # 绘制曲线图
        ax.plot([lymphoid], [dis_avg], marker='o', label=type1)  # 每个type1对应一条曲线，标签为type1的值
# 设置横坐标标签
ax.set_xticklabels(lymphoid_values)
# 添加图例
ax.legend()
# 设置纵坐标标签
ax.set_ylabel('dis')
# 保存为PDF文件
plt.savefig('Lym123.celltype.self.dis.output.pdf')
