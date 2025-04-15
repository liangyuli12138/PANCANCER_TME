library(pheatmap)
mycolors <- colorRampPalette(c("#5E4EA1","#5267AD","#3F7EB8","#4596B7","#5AAEAD","#6CC4A4","#8ACFA4","#A5DAA4","#BFE4A0","#D7EE9B","#EAF69E","#F4FAAE","#FFFFBF","#FFF2A9","#FEE593","#FED380","#FEBE6E","#FCA85E","#F98E51","#F47245","#E95E47","#DC4A4C","#CB364C","#B42147","#9E0041"))(25)
data = read.table("./input.pearson.xls",header = T,row.names=1,sep = "\t");
pheatmap(data,scale="none",fontsize_col=10,fontsize_row=10,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="average",file ='./output.pearson.average.pdf',width=8,height=10,col=mycolors)
pheatmap(data,scale="none",fontsize_col=10,fontsize_row=10,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="single",file ='./output.pearson.single.pdf',width=8,height=10,col=mycolors)
pheatmap(data,scale="none",fontsize_col=10,fontsize_row=10,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="complete",file ='./output.pearson.complete.pdf',width=8,height=10,col=mycolors)
pheatmap(data,scale="none",fontsize_col=10,fontsize_row=10,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="mcquitty",file ='./output.pearson.mcquitty.pdf',width=8,height=10,col=mycolors)
pheatmap(data,scale="none",fontsize_col=10,fontsize_row=10,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="median",file ='./output.pearson.median.pdf',width=8,height=10,col=mycolors)
pheatmap(data,scale="none",fontsize_col=10,fontsize_row=10,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="centroid",file ='./output.pearson.centroid.pdf',width=8,height=10,col=mycolors)
