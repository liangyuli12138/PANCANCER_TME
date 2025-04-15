import scanpy as sc
import matplotlib.pyplot as plt

os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/16.diff_all/plot")

adata_concat = sc.read_h5ad("pancancer.icar.all.cell.filter.h5ad")

marker_genes_dict = {
"Lymphoid1":["FDCSP","FCER2","CXCL13",],
"Lymphoid2":["C7","CXCL16","CXCL12",],
"Lymphoid3":["S100A8","CD2","S100A9","C3","IL10","HLA-A","HLA-B","HLA-C","CXCL9","CXCL1","CXCL11","CXCL10","CXCL14"],
}

adata_concat = adata_concat[adata_concat.obs['new_groups'].isin(['Lymphoid1', 'Lymphoid2','Lymphoid3']), :]
sc.pp.scale(adata_concat, max_value=10)



# 指定要绘制小提琴图的基因列表
genes = ["FDCSP", "FCER2", "C7", "CXCL13", "CXCL9", "CXCL5", "CXCL1", "CXCL2", "CXCL3", "CXCL8", "CXCL11", "CXCL17", "CXCL16", "CXCL10", "CXCL12", "CXCL14", "S100A8", "CD2", "S100A9", "C3", "IL10", "HLA-A", "HLA-B", "HLA-C"]

# 创建子图布局
fig, axes = plt.subplots(nrows=len(genes)//4 + (len(genes)%4 > 0), ncols=4, figsize=(12, 3*(len(genes)//4 + (len(genes)%4 > 0))))

# 绘制小提琴图
for i, gene in enumerate(genes):
    ax = axes[i//4, i%4]
#    adata_filtered = adata_concat[adata_concat.X[:, gene] != 0]
    sc.pl.violin(adata_filtered, gene, groupby="new_groups",  stripplot=False, nner="box",  ax=ax)
    ax.set_title(gene)

# 调整子图之间的间距
plt.tight_layout()

# 保存为PDF文件
plt.savefig("violin_plots.png")
