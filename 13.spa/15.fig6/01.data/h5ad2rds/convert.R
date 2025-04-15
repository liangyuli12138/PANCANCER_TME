library(sceasy)
library(SeuratDisk)
library(Seurat)
#Convert

args <- commandArgs(TRUE)
input_file <- args[1] 
from <- args[2]  ## anndata, seurat, sce, loom
to <- args[3]    ## anndata, seurat, sce, loom
outpfx <- args[4]

if (from == "anndata" & to=="seurat"){
    sceasy::convertFormat(input_file, from=from, to=to,
                       outFile=paste0(outpfx, '.rds'))

}

if (from=="seurat" & to=="anndata") {
    obj <- readRDS(input_file)
    #DefaultAssay(obj) <- "RNA"
    DefaultAssay(obj) <- "Spatial"
    new.obj <- CreateSeuratObject(obj@assays$Spatial@counts)
#    new.obj <- RenameCells(object = new.obj, add.cell.id = outpfx)

    #new.obj@meta.data <- obj@meta.data
    for(i in colnames(new.obj@meta.data)){
        new.obj@meta.data[[i]] <- as.character(new.obj@meta.data[[i]])
    }
    sceasy::convertFormat(new.obj, from=from, to=to,
                       outFile=paste0(outpfx, '.h5ad'))
}
