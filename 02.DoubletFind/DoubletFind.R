#使用R422环境
.libPaths("/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R411/lib/R/library/")
library(Seurat)
library(DoubletFinder)
args <- commandArgs(T)
input_rds <- args[1]
cancer <- args[2]
output <- args[3]
cat('Input_rds:',input_rds,'\n')
cat('Cancer:',cancer,'\n')
cat('outdir:',output,'\n')
data <- readRDS(input_rds)
data <- subset(data,subset=Tissue==cancer)

data_tmp <- SplitObject(data,split.by='batch')

data_result <- sapply(data_tmp,function(data_tmp1) {
    data_tmp1 <- NormalizeData(data_tmp1)
    data_tmp1 <- FindVariableFeatures(data_tmp1, selection.method = "vst", nfeatures = 2000)
    data_tmp1 <- ScaleData(data_tmp1)
    data_tmp1 <- RunPCA(data_tmp1)
    data_tmp1 <- RunUMAP(data_tmp1, dims = 1:30)
    sweep.res.list_skin <- paramSweep_v3(data_tmp1, PCs = 1:30, sct = FALSE) #sct=TRUE representing SCTransform was used
    sweep.stats_skin <- summarizeSweep(sweep.res.list_skin, GT = FALSE)
    bcmvn_skin <- find.pK(sweep.stats_skin)
    mpK<-as.numeric(as.vector(bcmvn_skin$pK[which.max(bcmvn_skin$BCmetric)]))
    nExp_poi <- round(0.075*nrow(data_tmp1@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
    data_tmp1<- doubletFinder_v3(data_tmp1, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    data_tmp1 <- data_tmp1@meta.data
    tmp_colname <- colnames(data_tmp1)
    tmp_colnames <- c(tmp_colname[1:(length(tmp_colname)-2)],'pANN','DF.classifications')
    colnames(data_tmp1) <- tmp_colnames
    return(list(data_tmp1))
})

data_result <- Reduce('rbind',data_result)
write.csv(data_result,paste0(output,'/',cancer,'_DoubletFind_result.csv'),quote = FALSE)
