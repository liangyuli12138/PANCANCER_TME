library(Seurat)
library(SeuratDisk)

fileh="/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/03.cluster_primary/anndata_0.7/h5ad2rds/pancancer.ref.0723.final.h5ad"

Convert(fileh, dest = "h5seurat", overwrite = FALSE)
pbm <- LoadH5Seurat("pancancer.ref.0723.final.h5seurat")
saveRDS(pbm,file="pancancer.ref.0723.final.rds")

