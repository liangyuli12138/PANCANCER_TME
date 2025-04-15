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
library(sf)

args<-commandArgs(T)

index <- "D01872D1"
outfix_sample <- "Colorectal.D01872D1"

#inpute
RCTD_file <- paste0('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/01.Hexcell/region_merge/D01872D1_region.csv')#cell2location结果文件

#cellbin数据文件
cellbin_data_file <- paste0('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/',index,'/',index,'_cellbin.final.rds')
#芯片号
#outfix_sample <- 'D01972B6'
#输出文件夹
output_path <- '/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/01.Hexcell/plot/plot'


####58个亚类颜色设置#####
#color_use <- read.csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/02.data/plot/colour.list')
color_use <- read.csv('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/01.Hexcell/plot/color.list')
cellType_cols <- color_use$Colour
names(cellType_cols) <- color_use$Sub_cluster

########################
output_path <- paste0(output_path,'/',outfix_sample)
dir.create(output_path,showWarnings = F,recursive = T)

RCTD_result <- read.csv(RCTD_file,header=TRUE,row.names=1)
cellbin_data <- readRDS(cellbin_data_file)
###数据过滤
cell_filter <- rownames(RCTD_result)
data <- subset(cellbin_data, cells=cell_filter)

#results_df_tmp <- RCTD_result$celltype
results_df_tmp <- RCTD_result[colnames(data),]
identical(rownames(results_df_tmp),colnames(data))
data@meta.data <- cbind(data@meta.data,results_df_tmp)
sf <- sf::st_as_sf(data@meta.data, wkt = c("contour"))

#polygon_sf <- st_read("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/03.net/out/D01872D1.original_polygons.shp")
#first ty.pe的空间分布
options(repr.plot.width =16,repr.plot.height =8)
plot2 <- ggplot() + geom_sf(data = sf,aes(fill=region),color=NA) + 
#geom_sf(data = polygon_sf, color = "red", fill = NA, size = 1, linetype = "solid", alpha = 0.7) +
scale_fill_manual(values = cellType_cols)+
#  scale_fill_manual(values = c("other" = "gray")) +
#  scale_fill_identity() +
scale_linetype_manual(values = c(polygon = "dashed")) +
theme_void()+theme(plot.title=element_text(hjust=0.5,vjust=0.5, size = 20),legend.position = "right")

pdf(paste0(output_path,'/',outfix_sample,'.ori_sf.pdf'),w=16,h=8)
plot2
dev.off()

plot2 <- ggplot() + geom_sf(data = sf,aes(fill=region),color=NA) +
#geom_sf(data = polygon_sf, color = "red", fill = NA, size = 1, linetype = "solid", alpha = 0.7) +
scale_fill_manual(values = cellType_cols)+
#  scale_fill_manual(values = c("other" = "gray")) +
#  scale_fill_identity() +
scale_linetype_manual(values = c(polygon = "dashed")) +
             theme_void()+theme(plot.title=element_text(hjust=0.5,vjust=0.5, size = 20))+ ggtitle(outfix_sample)+theme(legend.position="right")

png(paste0(output_path,'/',outfix_sample,'.ori.png'),res=300,w=8,h=8,units='in')
plot2
dev.off()

