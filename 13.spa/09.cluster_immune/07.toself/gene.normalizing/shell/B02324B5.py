import scanpy as sc
import pandas as pd

adata = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/B02324B5_cellbin.final.celltype.h5ad"
adata = sc.read_h5ad(adata)
sc.pp.filter_cells(adata, min_genes=100)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt < 20, :]
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

cellist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/gene/at/B02324B5.input")
atlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/gene/at/B02324B5.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)


# 读取基因列表文件
with open('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/gene/gene.list', 'r') as f:
    gene_list = [line.strip() for line in f.readlines()]

# 获取每个tlsclu下所有基因的平均表达量
group_info = pd.DataFrame(index=adata.obs['tlsclu'].unique(), columns=gene_list)
for tlsclu in adata.obs['tlsclu'].unique():
    for gene in gene_list:
        if gene in adata.var_names:
            mean_expression = adata[adata.obs['tlsclu'] == tlsclu].obs_vector(gene).mean()
        else:
            mean_expression = 0
        group_info.loc[tlsclu, gene] = mean_expression

# 保存为geneAmean.csv文件
group_info.to_csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/gene.normalizing/out/B02324B5.geneAmean.csv')
