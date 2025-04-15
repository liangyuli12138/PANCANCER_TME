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

index <- "D01972B6"
outfix_sample <- "Gastric.D01972B6"

#inpute
RCTD_file <- paste0('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/03.TLS/umap_all/at_cell/',index,'.cell')#cell2location结果文件

#cellbin数据文件
cellbin_data_file <- paste0('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/',index,'/',index,'_cellbin.final.rds')
#芯片号
#outfix_sample <- 'D01972B6'
#输出文件夹
output_path <- '/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/03.TLS/umap_all/plot_cell/'


####58个亚类颜色设置#####
#color_use <- read.csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/03.TLS/umap/colour.list')
color_use <- read.csv('/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/liuhui/project/pancancer/CellBin_FastMid/RCTD/results/plot/颜色设置marker.csv')
cellType_cols <- color_use$Colour
names(cellType_cols) <- color_use$Sub_cluster

########################
output_path <- paste0(output_path,'/',outfix_sample)
dir.create(output_path,showWarnings = F,recursive = T)

RCTD_result <- read.csv(RCTD_file,header=TRUE,row.names=2)
cellbin_data <- readRDS(cellbin_data_file)
###数据过滤
cellbin_data[["percent.mt"]] <- PercentageFeatureSet(cellbin_data, pattern = "^MT-")#线粒体基因比例
#max_feature <- max(cellbin_data$nFeature_Spatial)*0.95
data <- subset(cellbin_data, subset = nFeature_Spatial>100 & percent.mt < 15 & nFeature_Spatial<3000 & area>50 & area<2500)

#results_df_tmp <- RCTD_result$celltype
results_df_tmp <- RCTD_result[colnames(data),]
identical(rownames(results_df_tmp),colnames(data))
data@meta.data <- cbind(data@meta.data,results_df_tmp)
sf <- sf::st_as_sf(data@meta.data, wkt = c("contour"))

#first type的空间分布
options(repr.plot.width =16,repr.plot.height =8)
plot2 <- ggplot() + geom_sf(data = sf,aes(fill=celltype),color=NA) + 
             theme_void()+theme(plot.title=element_text(hjust=0.5,vjust=0.5))+
scale_fill_manual(values = cellType_cols)

pdf(paste0(output_path,'/',outfix_sample,'_RCTD_firstType_sf.pdf'),w=16,h=8)
plot2
dev.off()

plot2 <- ggplot() + geom_sf(data = sf,aes(fill=celltype),color=NA) +
             theme_void()+theme(plot.title=element_text(hjust=0.5,vjust=0.5))+ ggtitle(outfix_sample)+theme(legend.position="none")
scale_fill_manual(values = cellType_cols)

png(paste0(output_path,'/',outfix_sample,'_RCTD_firstType_sf.png'),res=300,w=8,h=8,units='in')
plot2
dev.off()

#plot2

#每种细胞类型的空间分布
SpatialPlot_SF_each <- function(plot_data,celltype_infor){
    plot_list <- list()
    for(cell_name in unique(celltype_infor)){
        index <- celltype_infor %in% cell_name
        plot_list[[cell_name]] <- ggplot() + geom_sf(data = sf[!index,],fill='grey',color=NA) +
             geom_sf(data = sf[index,],fill='red',color=NA) + ggtitle(cell_name)+
             theme_void()+theme(plot.title=element_text(hjust=0.5,vjust=0.5))
    }
   return(plot_list)
}
plot_list <- SpatialPlot_SF_each(plot_data=sf,celltype_infor=sf$celltype)
for(celltype in names(plot_list)){
  png(paste0(output_path,'/',outfix_sample,'_',celltype,'_Spatial.png'),res=300,w=8,h=8,units='in')
  print(plot_list[[celltype]])
  dev.off()
}
