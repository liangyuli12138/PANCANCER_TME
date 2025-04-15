library(Seurat)
library(SeuratDisk)

fileh="/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/01.all.cluster/oridata_allcancer/h5ad2rds/pancancer.ref.raw.0917.h5ad"

Convert(fileh, dest = "h5seurat", overwrite = FALSE)
pbm <- LoadH5Seurat("pancancer.ref.raw.0917.h5seurat")
saveRDS(pbm,file="pancancer.ref.raw.0917.rds")

