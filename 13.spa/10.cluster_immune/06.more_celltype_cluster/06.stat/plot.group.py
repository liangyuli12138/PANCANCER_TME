import pandas as pd
import matplotlib.pyplot as plt

# 清除之前的图像信息
plt.clf()

# 读取csv文件
data = pd.read_csv("group.stat.csv", delimiter="\t")

# 提取数据列
tissues = data["Tissue"]
lymphoid0 = data["Lymphoid0"]
lymphoid1 = data["Lymphoid1"]
lymphoid2 = data["Lymphoid2"]
lymphoid3 = data["Lymphoid3"]
myeloid4 = data["Myeloid4"]
myeloid5 = data["Myeloid5"]
myeloid6 = data["Myeloid6"]
myeloid7 = data["Myeloid7"]
myeloid8 = data["Myeloid8"]

# 调整图形大小
fig, ax = plt.subplots(figsize=(12, 6))

# 设置颜色
colors = plt.get_cmap('tab10')

# 绘制堆积柱状图，使用颜色循环
ax.bar(tissues, lymphoid0, label="Lymphoid0", color=colors(0))
ax.bar(tissues, lymphoid1, bottom=lymphoid0, label="Lymphoid1", color=colors(1))
ax.bar(tissues, lymphoid2, bottom=lymphoid0+lymphoid1, label="Lymphoid2", color=colors(2))
ax.bar(tissues, lymphoid3, bottom=lymphoid0+lymphoid1+lymphoid2, label="Lymphoid3", color=colors(3))
ax.bar(tissues, myeloid4, bottom=lymphoid0+lymphoid1+lymphoid2+lymphoid3, label="Myeloid4", color=colors(4))
ax.bar(tissues, myeloid5, bottom=lymphoid0+lymphoid1+lymphoid2+lymphoid3+myeloid4, label="Myeloid5", color=colors(5))
ax.bar(tissues, myeloid6, bottom=lymphoid0+lymphoid1+lymphoid2+lymphoid3+myeloid4+myeloid5, label="Myeloid6", color=colors(6))
ax.bar(tissues, myeloid7, bottom=lymphoid0+lymphoid1+lymphoid2+lymphoid3+myeloid4+myeloid5+myeloid6, label="Myeloid7", color=colors(7))
ax.bar(tissues, myeloid8, bottom=lymphoid0+lymphoid1+lymphoid2+lymphoid3+myeloid4+myeloid5+myeloid6+myeloid7, label="Myeloid8", color=colors(8))

# 添加图例和标签
ax.legend()
ax.set_xlabel("Tissue")
ax.set_ylabel("ICAR Quantity")

# 调整x轴标签角度
plt.xticks(rotation=45)

# 调整布局以确保横坐标标签完全显示
plt.tight_layout()

# 保存为pdf文件
plt.savefig("group.stat.pdf")

