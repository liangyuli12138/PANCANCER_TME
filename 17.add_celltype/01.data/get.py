import scanpy as sc

# 读取h5ad文件为adata
adata = sc.read_h5ad('pancancer.ref.umap.0723.h5ad')

# 提取B细胞的数据
adata_b = adata[adata.obs['groups_ori'] == 'B'].copy()
adata_b.write_h5ad('B_cells.h5ad')
adata_b.obs.to_csv('B_cells.obs.csv')
adata_b.var.to_csv('B_cells.var.csv')

# 提取EC细胞的数据
adata_ec = adata[adata.obs['groups_ori'] == 'EC'].copy()
adata_ec.write_h5ad('EC_cells.h5ad')
adata_ec.obs.to_csv('EC_cells.obs.csv')

# 提取Myeloid细胞的数据
adata_myeloid = adata[adata.obs['groups_ori'] == 'Myeloid'].copy()
adata_myeloid.write_h5ad('Myeloid_cells.h5ad')
adata_myeloid.obs.to_csv('Myeloid_cells.obs.csv')
