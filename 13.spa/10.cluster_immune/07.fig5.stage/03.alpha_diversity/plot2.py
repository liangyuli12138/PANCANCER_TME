import pandas as pd
import matplotlib.pyplot as plt

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/07.fig5.stage/03.alpha_diversity")
# 读取csv文件
df = pd.read_csv("cell.merge.at.group.new")

# 根据cell列进行过滤
df_lymphoid = df[df["cell"].str.contains("Lymphoid")]

# 计算每个基因下的平均值和方差
gene_expr_mean = df_lymphoid.groupby("cell").mean()

# 选择需要绘制的基因列表
genes = ["B_memory","B_naive","CD4_Tfh","CD4_Treg","T_naïve","DC","Marco"]

# 绘制折线图和方差图
mean_plot = gene_expr_mean.plot(kind="line", y=genes, marker="o")

mean_plot.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.tight_layout()  # 调整图像布局

plt.savefig('output.celltype.line.TLS.new.pdf', dpi=300)

