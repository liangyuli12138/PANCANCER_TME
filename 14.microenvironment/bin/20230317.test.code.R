################################################################################
# revised at 20230316

################################################################################
# pathway visualization
setwd("C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion")
library(Seurat)
library(ggplot2)
library(UCell)
library(pheatmap)

# Input gene list file
#gene.set <- "C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230302.test.set.yangyi\\AD_meta_pathway_geneset.csv"
#gene.set <- "C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\20221207.geneset.collection.ovary.csv"
gene.set <- "C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\20230317.test.pathway.ovary.csv"
gene.set <- read.csv(gene.set,header=FALSE)

gene.set.name <- unique(gene.set[,1])
markers <- list()
for (ii in 1:length(gene.set.name))
{
	markers[[ii]] <- gene.set[gene.set[,1]==gene.set.name[ii],2]
}
names(markers) <- gene.set.name


# Input 10x format stereo-seq data
#name = "C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230302.test.set.yangyi\\PTB469_1_A"
name = "C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\11y1"

# Output prefix
data.name <- "11y1_bin200"

# Generate spatial transcriptome Seurat object 
#obj <- readRDS("C:\\Users\\superman\\Desktop\\ShareData\\20221216.breast.cancer\\20221222.stereoseq.bin200\\20221221.bin200.30.stereoseq.combined.RDS")
#cell.subset <- obj[[11]]
data <- Read10X(data.dir = name)
cell.subset <- CreateSeuratObject(counts = data, project = name)
cell.subset[["percent.mt"]] <- PercentageFeatureSet(cell.subset, pattern = "^MT-")
cell.subset <- subset(cell.subset, subset = nFeature_RNA > 200 & percent.mt < 20)
cell.subset <- NormalizeData(cell.subset)
#cell.subset <- ScaleData(cell.subset)

counts <- cell.subset@assays$RNA@counts
slide.seq <- CreateSeuratObject(counts = counts, assay = 'Spatial')

cell.loc <- colnames(cell.subset)
cell.loc.2 <- strsplit(cell.loc,split="-")
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
	plot2 <- SpatialFeaturePlot(slide.seq, features = marker.name[ii],stroke=NA)+
		theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
	outpdf <- paste(data.name,".",marker.name[ii],".SpatialFeaturePlot.pdf",sep="")
	pdf(outpdf, w=6, h=8)
	print(plot2)
	while (!is.null(dev.list()))  dev.off()
}

csv.name <- paste(data.name,".UCell.score.csv",sep="")
write.csv(slide.seq@meta.data,csv.name)
# End
################################################################################


################################################################################
# test the pathway correlations
	mat <- read.csv(csv.name,header=TRUE)
	zz <- t(mat[,5:dim(mat)[2]])
	test.cell <- round(dim(zz)[2]*0.05)

	exp <- as.matrix(zz)	
	genes.1 <- rownames(exp)
	top.cell.list <- list()
	for (j in 1:length(genes.1))
	{
		tmp.1 <- exp[j,]
		tmp.3 <- sort(tmp.1,decreasing = TRUE)
	if(tmp.3[test.cell]>0)
	{
			thre <- tmp.3[test.cell]
		} else {
			thre <- 0	
		}	
		m <- which(tmp.1>thre)
		top.cell.list[[j]] <- m
	}
						
	a <- length(top.cell.list)
	marker.cor <- matrix(0, nrow =a , ncol = a)
	for (i in 1:a)
	{
		m.1 <- top.cell.list[[i]]			
		for (j in i:a)
		{
			if (i!=j)
			{
				m.2 <- top.cell.list[[j]]
				m.3 <- unique(intersect(m.1, m.2))
				m.4 <- unique(c(m.1, m.2))
				marker.cor[i,j] <- length(m.3)/max(1,length(m.4))
				marker.cor[j,i] <- length(m.3)/max(1,length(m.4))						
			}	
		}
	}
	rownames(marker.cor) <- genes.1
	colnames(marker.cor) <- genes.1

	marker.cor[marker.cor>quantile(marker.cor,0.99)] <- quantile(marker.cor,0.99)
	outpdf <- paste(data.name,".test.pathway.cor.heatmap.pdf",sep="")
	pdf(outpdf, w=30, h=30)
	pheatmap(marker.cor)
	while (!is.null(dev.list()))  dev.off()
################################################################################



















################################################################################
# bulk-RNA-seq analysis
set.A <- c(
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\2m",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\11y",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\11y1",
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\13w",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\17y",
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\18w",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\18y",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\18y4",
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\22w",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\24y",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\27y"
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\26w",
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\30w"
)

set.B <- c(
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\32y",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\33y",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\36y1",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\40y",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\42y",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\42y2",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\47y",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\49y",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\52y",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\56y",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\59y",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin200\\71y"
)

set <- c(set.A,set.B)
data <- Read10X(data.dir = set[1])
tmp.1 <- names(rowSums(data))

for (ii in 2:length(set))
{
	print(ii)
	data <- Read10X(data.dir = set[ii])
	tmp <- names(rowSums(data))
	tmp.1 <- c(tmp.1,tmp)
}
tmp.1 <- unique(tmp.1)

exp <- matrix(0,ncol=length(set),nrow=length(tmp.1))
rownames(exp) <- tmp.1
#gene <- tmp.1
#exp <- c()
for (ii in 1:length(set))
{
	print(ii)
	data <- Read10X(data.dir = set[ii])
	tmp <- rowSums(data)
	gene<- names(tmp)
	exp[gene,ii] <- tmp 
	#exp <- cbind(exp,tmp)
}

exp.norm <- exp/colSums(exp)*1000000
DEGs <- exp.norm[,c(1,2)]
for (ii in 1:dim(exp.norm)[1])
{
	print(ii)
	tmp <- wilcox.test(exp.norm[ii,1:7],exp.norm[ii,8:19])
	FC <- max(median(exp.norm[ii,1:7]),1)/max(median(exp.norm[ii,8:19]),1)
	DEGs[ii,1] <- FC
	DEGs[ii,2] <- tmp$p.value
}
write.csv(DEGs,"DEGs.1.csv")
genes.1 <- rownames(DEGs[(DEGs[,2]<0.05)&(DEGs[,1]>5),])



##############################
# stereo-seq gene module detection

 get test bin50 datasets
set.50 <- c(
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin50\\2m",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin50\\11y1",
"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin50\\17y"
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin50\\71y",
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin50\\49y",
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin50\\32y"
)

#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin50\\13w",
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin50\\18y",
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin50\\24y",
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin50\\33y",
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin50\\36y1",
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin50\\40y"
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin50\\47y",
#"C:\\Users\\superman\\Desktop\\20221006.Immu.Endo.Metabolic.GeneSet\\20230227.microenvironment.discussion\\20230303.test.set.yanghuanjie\\bin50\\59y"
#)
data <- Read10X(data.dir = set.50[1])
obj <- CreateSeuratObject(data)
for (ii in 2:length(set.50))
{
	print(ii)
	data <- Read10X(data.dir = set.50[ii])
	tmp <- CreateSeuratObject(data)
	obj <- merge(obj,tmp)
}
##############################
genes.1 <- read.csv("test.csv",header=FALSE)
genes.1 <- genes.1[,2]
	obj <- NormalizeData(obj)
# test the pathway correlations
	test.cell <- round(dim(obj)[2]*0.05)
	genes.1 <- intersect(rownames(obj),genes.1)
	exp <- as.matrix(obj@assays$RNA@data[genes.1,])

	top.cell.list <- list()
	for (j in 1:length(genes.1))
	{
		print(j)
		tmp.1 <- exp[j,]
		tmp.3 <- sort(tmp.1,decreasing = TRUE)
	if(tmp.3[test.cell]>0)
	{
			thre <- tmp.3[test.cell]
		} else {
			thre <- 0	
		}	
		m <- which(tmp.1>thre)
		top.cell.list[[j]] <- m
	}
						
	a <- length(top.cell.list)

if (FALSE)
{
##############################	
library(snowfall)
top.cell.list.len <- length(top.cell.list)
marker_cor_calculate <- function(xxx)
{
	res.vect <- c()
	m.1 <- top.cell.list[[xxx]]
	for (j in 1:top.cell.list.len)
	{
		m.2 <- top.cell.list[[j]]
		m.3 <- unique(intersect(m.1,m.2))
		m.4 <- unique(c(m.1, m.2))
		res.vect <- c(res.vect,length(m.3)/max(1,length(m.4)))
	}
	res.set <- list(xxx, res.vect)
}
sfInit(parallel = TRUE, cpus = 6) #初始化
top.cell.list.len <- length(top.cell.list)
	marker.cor <- matrix(0, nrow =top.cell.list.len , ncol = top.cell.list.len)
	sfExport("top.cell.list","marker.cor","top.cell.list.len")         # 载入依赖的对象
	result <- sfLapply(1:top.cell.list.len, marker_cor_calculate)
	
	for (i in 1:length(result))
	{
		m.1 <- as.matrix(result[[i]][[2]])
		m.1[result[[i]][[1]]] <- 0
		marker.cor[result[[i]][[1]],] <- m.1			
	}						
	rownames(marker.cor) <- genes.1
	colnames(marker.cor) <- genes.1		
	
sfStop()		
##############################	
}
	marker.cor <- matrix(0, nrow =a , ncol = a)
	for (i in 1:a)
	{
		print(i)
		m.1 <- top.cell.list[[i]]			
		for (j in i:a)
		{
			if (i!=j)
			{
				m.2 <- top.cell.list[[j]]
				m.3 <- unique(intersect(m.1, m.2))
				m.4 <- unique(c(m.1, m.2))
				marker.cor[i,j] <- length(m.3)/max(1,length(m.4))
				marker.cor[j,i] <- length(m.3)/max(1,length(m.4))						
			}	
		}
	}
	rownames(marker.cor) <- genes.1
	colnames(marker.cor) <- genes.1
##############################


	marker.cor[marker.cor>quantile(marker.cor,0.99)] <- quantile(marker.cor,0.99)
	data.name <- "DownDEGs.3"
	outpdf <- paste(data.name,".test.pathway.cor.heatmap.pdf",sep="")
	pdf(outpdf, w=24, h=24)
	#p <- pheatmap(marker.cor,clustering_method="ward.D2",cutree_rows=8)
	pheatmap(marker.cor,clustering_method="ward.D2")
	while (!is.null(dev.list()))  dev.off()
################################################################################

if (FALSE)
{
## 对基因进行分类
row_cluster <- cutree(p$tree_row,k=8)
newOrder <- marker.cor[p$tree_row$order,]
newOrder=cbind(newOrder,row_cluster[match(rownames(newOrder),names(row_cluster))])
colnames(newOrder)[ncol(newOrder)]="Cluster"

## 输出结果
write.table(newOrder,"cluster.list.out",quote=F,row.names=T,col.names=T,sep="\t")
}