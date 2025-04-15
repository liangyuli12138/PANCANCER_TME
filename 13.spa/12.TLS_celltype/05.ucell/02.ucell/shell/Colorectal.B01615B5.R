library(Seurat)
library(ggplot2)
library(UCell)
library(pheatmap)
.libPaths("/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R411/lib/R/library/")
suppressMessages(library(patchwork))
suppressMessages(library(ggplot2))
suppressMessages(library(spacexr))
library(viridis)
library(ggdark)
library(sf)


setwd("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/05.ucell/02.ucell/out/")
gene.set <- "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/05.ucell/02.ucell/tls.list"
gene.set <- read.csv(gene.set,header=FALSE)
gene.set.name <- unique(gene.set[,1])
markers <- list()
for (ii in 1:length(gene.set.name))
{       
        markers[[ii]] <- gene.set[gene.set[,1]==gene.set.name[ii],2]
}
names(markers) <- gene.set.name

B01615B5 <- readRDS("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/B01615B5/B01615B5_cellbin.final.rds")
data.name <- "Colorectal.B01615B5_cellbin"
cell.subset <- B01615B5
cell.subset[["percent.mt"]] <- PercentageFeatureSet(cell.subset, pattern = "^MT-")
cell.subset <- subset(cell.subset, subset = nFeature_Spatial > 200 & percent.mt < 15)
cell.subset <- NormalizeData(cell.subset)

counts <- cell.subset@assays$Spatial@counts
slide.seq <- CreateSeuratObject(counts = counts, assay = 'Spatial')
cell.loc <- colnames(cell.subset)
cell.loc.2 <- strsplit(cell.loc,split="_")
cell.len <- length(cell.loc.2)
cell.axis <- matrix(0,nrow=cell.len, ncol=2)
for (j in 1:cell.len)
{
	cell.axis[j,1]=as.numeric(cell.loc.2[[j]][1])
	cell.axis[j,2]=as.numeric(cell.loc.2[[j]][2])
}

coord.df = data.frame(x=-cell.axis[,2], y=cell.axis[,1], stringsAsFactors=FALSE) # (stringsAsFactors only if also have a separate barcodes column)
rownames(coord.df) <- colnames(cell.subset)

slide.seq@images$image =  new(
    Class = 'SlideSeq',
    assay = "Spatial",
    key = "image_",
    coordinates = coord.df
)
slide.seq <- NormalizeData(slide.seq)
slide.seq <- AddModuleScore_UCell(slide.seq, features = markers,ncores=4,force.gc=T, slot= "data")

polygon_sf <- st_read("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/03.net/out/B01615B5.original_polygons.shp")

marker.name <- colnames(slide.seq@meta.data)
index <- grep("_UCell",marker.name)
marker.name <- marker.name[index]
for (ii in 1:length(marker.name))
{
	# Get the sequence depth plot for input data
        plot2 <- SpatialFeaturePlot(slide.seq, features = marker.name[ii],stroke=NA,pt.size.factor=0.5)+
                theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
        outpdf <- paste(data.name,".",marker.name[ii],".SpatialFeaturePlot.pdf",sep="")
        pdf(outpdf, w=6, h=8)
        print(plot2)
        while (!is.null(dev.list()))  dev.off()
        plot3 <- SpatialFeaturePlot(slide.seq, features = marker.name[ii],stroke=NA,pt.size.factor=0.5)+
                theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
                theme(plot.title = element_text(size = 200))+
                theme(legend.key.size = unit(2, "cm"))
        outpng <- paste(data.name,".",marker.name[ii],".SpatialFeaturePlot.png",sep="")
        png(outpng,width = 1200, height = 1600)
        print(plot3)
        while (!is.null(dev.list()))  dev.off()
        plot4 <- SpatialFeaturePlot(slide.seq, features = marker.name[ii],stroke=NA,pt.size.factor=0.5)+
                theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
                theme(plot.title = element_text(size = 200))+
                theme(legend.key.size = unit(2, "cm"))
        polygon_plot <- geom_sf(data = polygon_sf, color = "red", fill = NA, size = 1, linetype = "solid", alpha = 0.7)
        combined_plot <- plot4 + polygon_plot
        outpng <- paste(data.name,".",marker.name[ii],".sf.SpatialFeaturePlot.png",sep="")
        png(outpng,width = 1200, height = 1600)
        print(combined_plot)
        while (!is.null(dev.list()))  dev.off()
}

csv.name <- "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/05.ucell/02.ucell/out/Colorectal.B01615B5.UCell.score.csv"
write.csv(slide.seq@meta.data,csv.name)
# End
