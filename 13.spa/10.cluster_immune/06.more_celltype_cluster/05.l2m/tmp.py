fig, ax = plt.subplots(figsize=(8,6))  # 修改图像大小为(8,6)
ax = sns.boxplot(data, x='groups', y='aaaa', palette='tab10', linewidth=0.8, order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'], showfliers=False)
sns.stripplot(data, x='groups', y='aaaa', jitter=0.2, size=2, alpha=0.8, color='black', order=['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3','Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])
plt.axhline(y=0.4, linestyle='dashed', color='red', zorder=0)   # 添加纵轴上的虚线
plt.xticks(rotation=45)  # 将横轴的标签字体倾斜45度
ax.set_xlabel('Groups', fontsize=15)
ax.set_ylabel('aaaa', fontsize=15)
ax.yaxis.set_label_coords(-0.08, 0.5)
plt.tight_layout()  # 调整图像布局
plt.subplots_adjust(top=0.9)
# 进行差异度t检验
group1 = data[data['groups'].isin(['Lymphoid0','Lymphoid1','Lymphoid2','Lymphoid3'])]['aaaa']
group2 = data[data['groups'].isin(['Myeloid4','Myeloid5','Myeloid6','Myeloid7','Myeloid8'])]['aaaa']
t_statistic, p_value = stats.ttest_ind(group1, group2)
# 在图表上添加p值
p_value_text = f"p value of t-test between Lymphoid and Myeloid groups: {p_value:.4f}"
plt.text(0.5, 1.03, p_value_text, transform=ax.transAxes, ha='center', fontsize=12)
plt.savefig('output.aaaa.png', dpi=300)

