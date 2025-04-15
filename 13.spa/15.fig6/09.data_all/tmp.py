adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/aaaa_cellbin.final.celltype.h5ad")
atlist = pd.read_csv("at100.group/aaaa.ex.at",index_col=0)
adata.obs = adata.obs.join(atlist)
adata.obs.index = "aaaa_" + adata.obs.index
adata_concat = adata_concat.concatenate(adata, join="outer")
del adata

