library(Seurat)
library(SeuratDisk)

fileh="/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/01.data/h5ad2rds/pancancer.icar.all.cell.h5ad"

Convert(fileh, dest = "h5seurat", overwrite = FALSE)
pbm <- LoadH5Seurat("pancancer.icar.all.cell.h5seurat")
saveRDS(pbm,file="pancancer.icar.all.cell.rds")

