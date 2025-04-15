.libPaths("/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R411/lib/R/library/")
#suppressMessages(library(SingleCellExperiment,lib.loc="/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R4/lib/R/library"))
#suppressMessages(library(monocle3))
suppressMessages(library(Seurat))
#suppressMessages(library(SeuratWrappers))
suppressMessages(library(patchwork))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(spacexr))
library(viridis)
library(ggdark)

args<-commandArgs(T)

index <- "bbbb"

cellbin_data_file <- paste0('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/12.spa/05.cell_correct/result/',index,'/',index,'_cellbin.final.rds')

output_path <- paste0('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/12.spa/05.cell_correct/result/',index,'/',index,'_cellbin.final.qc.csv')

cellbin_data <- readRDS(cellbin_data_file)

###数据过滤
cellbin_data[["percent.mt"]] <- PercentageFeatureSet(cellbin_data, pattern = "^MT-")#线粒体基因比例
max_feature <- max(cellbin_data$nFeature_Spatial)*0.9
data <- subset(cellbin_data, subset = nFeature_Spatial>100 & percent.mt < 15 & nFeature_Spatial<max_feature & area>50 & area<2500)

write.csv(data@meta.data, output_path)

