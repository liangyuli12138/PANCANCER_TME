import scanpy as sc
import pandas as pd
import os
import anndata
import harmonypy as hm


os.chdir("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/03.cluster_split/01.B")

# 读取h5ad文件和注释的csv文件，并进行处理
def process_sample(h5ad_file, csv_file, sample_id):
    # 读取h5ad文件
    adata = sc.read(h5ad_file)
    print("Sample ID:", sample_id, end="\n")    
    # 读取注释的csv文件
    annotation = pd.read_csv(csv_file)
    
    # 筛选出符合条件的细胞
    filtered_cells = annotation[(annotation['celltype'].str.contains('Lymphoid_B')) & (annotation['LM'] == 'Lymphoid')]
    if len(filtered_cells) == 0:
        print("Filtered cells is empty, skipping sample:", sample_id, end="\n")
        return None
    
    annotation = filtered_cells.copy()
 
    # 提取符合条件的cell的id
    cell_ids = filtered_cells['cell'].tolist()
    
    # 根据cell id进行子集选取
    adata_subset = adata[adata.obs_names.isin(cell_ids)].copy()
    
    # 过滤掉gene id包含RPL或RPS开头的核糖体基因
    filtered_genes = [gene for gene in adata_subset.var_names if not gene.startswith('RPL') and not gene.startswith('RPS')]
    adata_subset = adata_subset[:, filtered_genes]

    filtered_cells_number = len(adata_subset.obs_names)
    print("Filtered Cells:", filtered_cells_number, end="\n")

    # 给这些cell重新命名为：sample id+原cell id的格式
    adata_subset.obs_names = sample_id + '_' + adata_subset.obs_names
    
    # 在adata.obs文件里加一列为”sample“，记录样品id
    adata_subset.obs['sample'] = sample_id
    
    for col in annotation.columns:
        adata_subset.obs[col] = annotation[col].values

    # 删除处理好的adata对象
    del adata
    
    return adata_subset

# 处理所有样品并合并
def process_samples(h5ad_directory, csv_directory):
    # 获取所有样品的文件路径
    h5ad_files = [os.path.join(h5ad_directory, f) for f in os.listdir(h5ad_directory) if f.endswith('_cellbin.final.h5ad')]
    csv_files = [os.path.join(csv_directory, f) for f in os.listdir(csv_directory) if f.endswith('.tls.at')]

    # 合并所有处理好的adata对象
    merged_adata = None

    for h5ad_file in h5ad_files:
        # 获取样品id
        sample_id = os.path.splitext(os.path.basename(h5ad_file))[0].split('_')[0]

        # 获取对应的csv文件
        csv_file = None
        for file in csv_files:
            if sample_id in file:
                csv_file = file
                break
#import anndata as ad
#ad.concat([])
        if csv_file is not None:
            # 处理单个样品
            adata_subset = process_sample(h5ad_file, csv_file, sample_id)
	    
            # 合并adata对象
            if adata_subset is not None:
                if merged_adata is None:
                    merged_adata = adata_subset
                else:
                    merged_adata = merged_adata.concatenate(adata_subset,batch_key=None,index_unique=False)

            # 删除处理好的adata对象
            del adata_subset

    return merged_adata

# 设置不同的leiden参数进行聚类
def run_leiden_clustering(adata, leiden_resolutions):
    for resolution in leiden_resolutions:
        sc.tl.leiden(adata, key_added='tlscell', resolution=resolution)
        sc.pl.umap(adata, color='tlscell', save=f'leiden_resolution_{resolution}.png')
        adata.write(f'leiden_resolution_{resolution}.h5ad')
        adata.obs.to_csv(f'leiden_resolution_{resolution}_adata_obs.csv')

# 设置输出log文件
log_file = open("log.txt", "w")

# 设置log输出到文件
sc.settings.verbosity = 3
sc.settings.log_file = log_file

# 读取数据并进行处理
h5ad_directory = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/01.cluster/h5ad"
csv_directory = "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/01.cluster/at"
merged_adata = process_samples(h5ad_directory, csv_directory)

# 数据处理
sc.pp.normalize_total(merged_adata)
sc.pp.log1p(merged_adata)
sc.pp.highly_variable_genes(merged_adata, n_top_genes=2000)
sc.pp.scale(merged_adata, max_value=10)
sc.pp.pca(merged_adata)

meta_data = merged_adata.obs
data_mat = merged_adata.obsm["X_pca"]
ho = hm.run_harmony(data_mat, meta_data, "sample", theta=3)
merged_adata.obsm["X_harmony"] = ho.Z_corr.T

sc.pp.neighbors(merged_adata, use_rep = "X_harmony")
sc.tl.umap(merged_adata)

# 设置不同的leiden参数进行聚类并输出结果
leiden_resolutions = [0.3, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0]
c.pl.umap(adata_concat, color=genes,
 s=2, frameon=False, ncols=4, vmax='p99.9', na_color='grey' ,save="st.tls.umap.markergene.target.all.png")
un_leiden_clustering(merged_adata, leiden_resolutions)

log_file.close()

