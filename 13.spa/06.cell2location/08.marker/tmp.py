adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/aaaa_cellbin.final.celltype.h5ad")
cellist = pd.read_csv("aaaa_cellbin.final.celltype.obs.csv.input")
atlist = pd.read_csv("aaaa_cellbin.final.celltype.obs.csv.at",index_col=0)
adata = adata[cellist["cell"],:]
adata.obs = adata.obs.join(atlist)
adata_concat = adata_concat.concatenate(adata)

