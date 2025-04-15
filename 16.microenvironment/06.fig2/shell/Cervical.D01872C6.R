library(Seurat)
library(ggplot2)
library(UCell)
library(pheatmap)

setwd("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/16.microenvironment/06.fig2/out/")
gene.set <- "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/16.microenvironment/06.fig2/REACTOME.csv.input"
gene.set <- read.csv(gene.set,header=FALSE)
gene.set.name <- unique(gene.set[,1])
markers <- list()
for (ii in 1:length(gene.set.name))
{       
        markers[[ii]] <- gene.set[gene.set[,1]==gene.set.name[ii],2]
}
names(markers) <- gene.set.name

D01872C6 <- readRDS("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/D01872C6/D01872C6_cellbin.final.rds")
data.name <- "Cervical.D01872C6_cellbin"
cell.subset <- D01872C6
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
}

csv.name <- "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/16.microenvironment/06.fig2/out/Cervical.D01872C6.UCell.score.csv"
write.csv(slide.seq@meta.data,csv.name)
# End
