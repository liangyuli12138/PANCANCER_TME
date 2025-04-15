import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/10.cluster_immune/06.more_celltype_cluster/07.gene_20240701/plot.gc")

data = pd.read_csv('merge.all.ucell.list.csv', index_col=0)

def plot_box_strip(data, lymphoid_col):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax = sns.boxplot(data=data, x='groups', y=lymphoid_col, palette='tab10', linewidth=0.8, order=['Lymphoid0', 'Lymphoid1', 'Lymphoid2', 'Lymphoid3', 'Myeloid4', 'Myeloid5', 'Myeloid6', 'Myeloid7', 'Myeloid8'], showfliers=False)
    sns.stripplot(data=data, x='groups', y=lymphoid_col, jitter=0.2, size=2, alpha=0.8, color='black', order=['Lymphoid0', 'Lymphoid1', 'Lymphoid2', 'Lymphoid3', 'Myeloid4', 'Myeloid5', 'Myeloid6', 'Myeloid7', 'Myeloid8'])
    plt.xticks(rotation=45)
    ax.set_xlabel(lymphoid_col, fontsize=15)
    ax.set_ylabel('GC signature socre', fontsize=15)
    ax.yaxis.set_label_coords(-0.08, 0.5)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    group1 = data[data['groups'].isin(['Lymphoid0', 'Lymphoid1', 'Lymphoid2', 'Lymphoid3'])][lymphoid_col]
    group2 = data[data['groups'].isin(['Myeloid4', 'Myeloid5', 'Myeloid6', 'Myeloid7', 'Myeloid8'])][lymphoid_col]
    t_statistic, p_value = stats.ttest_ind(group1, group2)
    p_value_text = f"p value of t-test between Lymphoid and Myeloid groups: {p_value:.4f}"
    plt.text(0.5, 1.03, p_value_text, transform=ax.transAxes, ha='center', fontsize=12)
    plt.savefig('output.'+lymphoid_col+'.pdf', dpi=300)

lymphoid_cols = ['GC_ucell']
for lymphoid_col in lymphoid_cols:
    plot_box_strip(data, lymphoid_col)

