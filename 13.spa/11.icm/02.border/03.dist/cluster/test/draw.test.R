library(pheatmap)
mycolors <- colorRampPalette(c("#5E4EA1","#5267AD","#3F7EB8","#4596B7","#5AAEAD","#6CC4A4","#8ACFA4","#A5DAA4","#BFE4A0","#D7EE9B","#EAF69E","#F4FAAE","#FFFFBF","#FFF2A9","#FEE593","#FED380","#FEBE6E","#FCA85E","#F98E51","#F47245","#E95E47","#DC4A4C","#CB364C","#B42147","#9E0041"))(25)
data = read.table("./pearson.merge.csv",header = T,row.names=1,sep = "\t");

pheatmap(data,scale="none",fontsize_col=10,fontsize_row=10,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="average",file ='./test.png/output.pearson.average.png',width=8,height=10,col=mycolors)
pheatmap(data,scale="none",fontsize_col=10,fontsize_row=10,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="single",file ='./test.png/output.pearson.single.png',width=8,height=10,col=mycolors)
pheatmap(data,scale="none",fontsize_col=10,fontsize_row=10,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="complete",file ='./test.png/output.pearson.complete.png',width=8,height=10,col=mycolors)
pheatmap(data,scale="none",fontsize_col=10,fontsize_row=10,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="mcquitty",file ='./test.png/output.pearson.mcquitty.png',width=8,height=10,col=mycolors)
pheatmap(data,scale="none",fontsize_col=10,fontsize_row=10,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="median",file ='./test.png/output.pearson.median.png',width=8,height=10,col=mycolors)
pheatmap(data,scale="none",fontsize_col=10,fontsize_row=10,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="centroid",file ='./test.png/output.pearson.centroid.png',width=8,height=10,col=mycolors)
pheatmap(data,scale="none",fontsize_col=10,fontsize_row=10,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="ward.D",file ='./test.png/output.pearson.ward.D.png',width=8,height=10,col=mycolors)
pheatmap(data,scale="none",fontsize_col=10,fontsize_row=10,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="ward.D2",file ='./test.png/output.pearson.ward.D2.png',width=8,height=10,col=mycolors)

