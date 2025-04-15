import os 
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns 
import scipy.stats as stats

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/16.microenvironment/06.fig2/select")

data = pd.read_csv('select.REACTOME.csv', index_col=0)

data = data.sort_values(by='variable') # 按照variable列进行排序

unique_variables = data['variable'].unique() # 获取variable列的唯一值

order = sorted(sub_data['sample'].unique())

for variable in unique_variables: 
    fig, ax = plt.subplots(figsize=(8,6)) # 修改图像大小为(8,6) 
    sub_data = data[data['variable'] == variable] # 获取当前variable的子数据集
    ax = sns.boxplot(data=sub_data, x='sample', y='score', palette='tab10', linewidth=0.8, showfliers=False, order=order)
    plt.xticks(rotation=90)  # 将横轴的标签字体倾斜45度
    ax.set_xlabel('sample', fontsize=15)
    ax.set_ylabel('UCell score', fontsize=15)
    ax.yaxis.set_label_coords(-0.08, 0.5)
    ax.set_title(variable)
    kw_result = stats.kruskal(*[sub_data[sub_data['sample'] == sample]['score'] for sample in sub_data['sample'].unique()])
    kw_text = f"Kruskal-Wallis, p = {kw_result.pvalue:.4f}"
    plt.text(0.5, 0.95, kw_text, ha='center', va='top', transform=ax.transAxes, fontsize=12)
    plt.tight_layout()  # 调整图像布局
    plt.subplots_adjust(top=0.9)
    plt.savefig(f'REACTOME_ucell/output_{variable}.png', dpi=300)  # 根据variable的值生成对应的文件名，保存图像
