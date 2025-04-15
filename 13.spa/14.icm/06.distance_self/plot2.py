for cell_type in cell_types:
    # 筛选出特定细胞类型的数据
    cell_data = data[data['type1'] == cell_type]
    # 创建一个绘图窗口
    fig, ax = plt.subplots()
    # 绘制箱线图
    sns.boxplot(x='lym', y='dis', data=cell_data, ax=ax)
    # 设置图形标题
    ax.set_title(f'{cell_type} Cell Type')
    # 保存图形为PDF文件
    plt.savefig(f'{cell_type}_output.pdf')
    plt.close(fig)
