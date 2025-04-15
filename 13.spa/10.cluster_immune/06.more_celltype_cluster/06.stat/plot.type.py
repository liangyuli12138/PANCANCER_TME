import pandas as pd
import matplotlib.pyplot as plt
import os

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/06.stat")

# 读取csv文件
data = pd.read_csv("type.stat.csv", delimiter="\t")
plt.clf()

# 提取数据列
tissues = data["Tissue"]
lymphoid = data["Lymphoid"]
myeloid = data["Myeloid"]

# 绘制堆积柱状图
plt.bar(tissues, lymphoid, label="Lymphoid")
plt.bar(tissues, myeloid, bottom=lymphoid, label="Myeloid")

# 添加图例和标签
plt.legend()
plt.xlabel("Tissue")
plt.ylabel("ICAR Quantity")

# 调整x轴标签角度
plt.xticks(rotation=45)
plt.tight_layout()

# 保存为pdf文件
plt.savefig("type.stat.pdf")
