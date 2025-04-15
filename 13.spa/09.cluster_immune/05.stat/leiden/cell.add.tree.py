adata_concat.obs['cell_type_sum'] = np.sum(adata_concat[:, ['Lymphoid_CD4_Treg', 'Lymphoid_CD8_Tex', 'Myeloid_pDC']].X, axis=1)
color_map = {'Lymphoid_CD4_Treg': 'green', 'Lymphoid_CD8_Tex': 'red', 'Myeloid_pDC': 'blue'}
sc.pl.umap(adata_concat, color=['cell_type_sum'], palette=color_map, frameon=False, na_color="grey", save=".combined.cell.png")

