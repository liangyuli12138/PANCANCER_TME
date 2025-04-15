import pandas as pd
import matplotlib.pyplot as plt

# 读取csv文件
df = pd.read_csv("cell.merge.at.group")

# 根据cell列进行过滤
df_lymphoid = df[df["cell"].str.contains("Lymphoid")]

# 计算每个基因下的平均值和方差
gene_expr_mean = df_lymphoid.groupby("cell").mean()

# 选择需要绘制的基因列表
genes = ["Lymphoid_B_memory", "Lymphoid_B_naive", "Lymphoid_T", "Lymphoid_T_inhibitory", "Myeloid_cDC", "Myeloid_Marco_C1QC", "Myeloid_Marco_LYVE1", "Myeloid_Marco_SPP1"]

# 绘制折线图和方差图
mean_plot = gene_expr_mean.plot(kind="line", y=genes, marker="o")

mean_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.tight_layout()  # 调整图像布局

plt.savefig('output.celltype.line.TLS.png', dpi=300)

