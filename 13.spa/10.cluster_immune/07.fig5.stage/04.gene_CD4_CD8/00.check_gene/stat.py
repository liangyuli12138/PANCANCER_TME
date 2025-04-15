# 根据stage和celltype进行分组，并统计每个组中的cell数目
cell_counts = adata.obs[adata.obs['celltype'].isin(['Lymphoid_CD4_Treg', 'Lymphoid_CD8_Tex'])].groupby(['stage', 'celltype']).size()

# 统计每个stage的细胞总数
stage_counts = adata.obs.groupby('stage').size()

# 计算每个stage下，celltype为Lymphoid_CD4_Treg和Lymphoid_CD8_Tex的细胞数量占比
percentages = cell_counts / stage_counts * 100

# 打印数量和百分比结果
for index, count in cell_counts.items():
    stage = index[0]
    celltype = index[1]
    percentage = percentages[index]
    total_count = stage_counts[stage]
    print(f"Stage: {stage}, Celltype: {celltype}, Count: {count}, Percentage: {percentage:.2f}%, Total Count: {total_count}")

