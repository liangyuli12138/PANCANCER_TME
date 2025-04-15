adata = sc.read_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/aaaa/aaaa_cellbin.tmp.h5ad")
idlist = pd.read_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/aaaa/aaaa_cellbin.tmp.obs.id",index_col=0)
adata.obs.index = idlist.center_xy.values.astype(str)
adata.write_h5ad("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/aaaa/aaaa_cellbin.final.h5ad")
adata.obs.to_csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/aaaa/aaaa_cellbin.final.obs")

