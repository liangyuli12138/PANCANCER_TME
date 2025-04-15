library(Seurat)
library(ggplot2)
library(UCell)
library(pheatmap)

setwd("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/16.microenvironment/04.hallmarker.all/scale/out/")

B02324A5 <- readRDS("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/B02324A5/B02324A5_cellbin.final.rds")
data.name <- "Kidney.B02324A5_cellbin"
cell.subset <- B02324A5
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
#slide.seq <- AddModuleScore_UCell(slide.seq, features = markers,ncores=4,force.gc=T, slot= "data")
checklist <- read.csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/16.microenvironment/04.hallmarker.all/scale/split/Kidney.B02324A5.UCell.score.csv",header=TRUE)
rownames(checklist) <- checklist$X
checklist <- checklist[, -1]
slide.seq@meta.data <- checklist

marker.name <- colnames(slide.seq@meta.data)
index <- grep("_UCell",marker.name)
marker.name <- marker.name[index]
for (ii in 1:length(marker.name))
{
	# Get the sequence depth plot for input data
	plot2 <- SpatialFeaturePlot(slide.seq, features = marker.name[ii],stroke=NA,pt.size.factor=0.5, min.cutoff=-1, max.cutoff=1)+  
		theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
	outpdf <- paste(data.name,".",marker.name[ii],".SpatialFeaturePlot.pdf",sep="")
	pdf(outpdf, w=6, h=8)
	print(plot2)
	while (!is.null(dev.list()))  dev.off()
        plot3 <- SpatialFeaturePlot(slide.seq, features = marker.name[ii],stroke=NA,pt.size.factor=0.5, min.cutoff=-1, max.cutoff=1)+
                theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
                theme(plot.title = element_text(size = 100))+
                theme(legend.key.size = unit(2, "cm"))
        outpng <- paste(data.name,".",marker.name[ii],".SpatialFeaturePlot.png",sep="")
        png(outpng,width = 1200, height = 1600)
        print(plot3)
        while (!is.null(dev.list()))  dev.off()
}

#csv.name <- "/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/16.microenvironment/04.hallmarker.all/out/Kidney.B02324A5.UCell.score.csv"
#write.csv(slide.seq@meta.data,csv.name)
# End
