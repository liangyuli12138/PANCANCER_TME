import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/14.icm/02.border/plot/allcell/celltype_merge_zz")
data = pd.read_csv('all.cell.dist.csv.type.ICM').query("-2000 < distance < 2000")
packs = [(int(i), int(i+250)) for i in np.linspace(-2000, 2000-250, num=16)]
packs_labels = [f'{int(start)}-{int(end)}' for start, end in packs]

def _assign_label(df):
    for start, end in packs:
        if start<df['distance']<=end:
            label = f'{int(start)}-{int(end)}'
    return label

data['label'] = data.apply(_assign_label, axis=1)
data = data[data.celltype!='Unknown']
counts = (data.groupby(['label'])['celltype']
          .value_counts(normalize=True)
          .to_frame(name='ratio')
          .reset_index()
         )

mapping = data.groupby(['celltype_merge']).celltype.apply(lambda x: sorted(list(set(x)))).to_dict()
counts.label = pd.Categorical(counts.label, categories=packs_labels)

#color_df = pd.read_csv('celltype.color.csv')

for k, v in mapping.items():
    _,ax = plt.subplots(figsize=(6,4))
    df_to_plot = counts[counts.celltype.isin(v)]
    df_to_plot = df_to_plot.merge(color_df, left_on='celltype', right_on='Sub_cluster', how='left')
    ax=sns.lineplot(df_to_plot, x='label', y='ratio', hue='celltype', ax=ax)
    # Add dashed line between "-250--0" and "0-250"
    ax.axvline(x=7.5, linestyle='dashed', color='red')
    plt.xlabel('Distance Range')
    plt.ylabel('Percentage')
    plt.title(f'Cell Type: {k}')
    # sns.scatterplot(df_to_plot, x='label', y='ratio',ax=ax, color='black')
    for t in ax.get_xticklabels():
        t.set(rotation=90)
    handles, labels = ax.get_legend_handles_labels()  # 获取图例的标签和句柄
    sorted_labels, sorted_handles = zip(*sorted(zip(labels, handles)))  # 按字母顺序对标签进行排序
    plt.legend(sorted_handles, sorted_labels, loc='center left', bbox_to_anchor=(1.02, 0.5), fontsize=6)  # 将图例放置在右边外部
    plt.tight_layout()
#    plot_legend(ax=ax)
    # Save figure as png and pdf
    plt.savefig(f'png.out/{k}.png', dpi=300)
    plt.savefig(f'png.out/{k}.pdf')
    plt.close("all")
