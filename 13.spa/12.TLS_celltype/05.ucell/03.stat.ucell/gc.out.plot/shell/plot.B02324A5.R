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

index <- "B02324A5"
outfix_sample <- "Kidney.B02324A5"

#inpute
RCTD_file <- paste0('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/05.ucell/03.stat.ucell/gc.out/B02324A5.UCell.score.csv.filter')
cellbin_data_file <- paste0('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/',index,'/',index,'_cellbin.final.rds')
output_path <- '/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/05.ucell/03.stat.ucell/gc.out.plot/tls.gc.plot'

output_path <- paste0(output_path,'/',outfix_sample)
dir.create(output_path,showWarnings = F,recursive = T)

RCTD_result <- read.csv(RCTD_file,header=TRUE,row.names=1)
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

polygon_sf <- st_read("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/03.net/out/B02324A5.original_polygons.shp")
#first ty.pe的空间分布
options(repr.plot.width =16,repr.plot.height =8)

plot_function <- function(sf, polygon_sf, aaaa, bbbb, cccc, dddd, eeee, output_file){
   max_value <- aaaa * max(sf$GC_UCell, na.rm = TRUE)
   min_value <- bbbb * max(sf$GC_UCell, na.rm = TRUE)
   plot2 <- ggplot() + 
     geom_sf(data = sf, aes(fill = GC_UCell), color = NA) +
     scale_fill_gradient(low = cccc, high = dddd, na.value = eeee, limits = c(min_value, max_value)) +
     labs(fill = "GC score") +
     geom_sf(data = polygon_sf, color = "red", fill = NA, size = 1, linetype = "solid", alpha = 0.7) +
     scale_linetype_manual(values = c(polygon = "dashed")) +
     theme_void() +
     theme(plot.title = element_text(hjust = 0.5, vjust = 0.5))
   png(output_file, res = 300, w = 8, h = 8, units = 'in')
   print(plot2)
   dev.off()
}

output_file <- "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/05.ucell/03.stat.ucell/gc.out.plot/tls.gc.plot/B02324A5.png"

plot <- plot_function(sf, polygon_sf,0.8,0.3,"lightblue","#800080","lightgrey",output_file)

