library(pheatmap)
mycolors <- colorRampPalette(c("#5E4EA1","#5267AD","#3F7EB8","#4596B7","#5AAEAD","#6CC4A4","#8ACFA4","#A5DAA4","#BFE4A0","#D7EE9B","#EAF69E","#F4FAAE","#FFFFBF","#FFF2A9","#FEE593","#FED380","#FEBE6E","#FCA85E","#F98E51","#F47245","#E95E47","#DC4A4C","#CB364C","#B42147","#9E0041"))(25)
data = read.table("./B.monify.xls",header = T,row.names=1,sep = "\t");
pheatmap(data,scale="none",fontsize_col=5,fontsize_row=5,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="ward.D2",file ='./pancancer.ref.0905.final.obs.csv.from.windows.zscore.B.1.pdf',width=8,height=10,col=mycolors)
pheatmap(data,scale="none",fontsize_col=5,fontsize_row=5,border_color = "grey",cluster_rows=T,cluster_cols=T,file ='./pancancer.ref.0905.final.obs.csv.from.windows.zscore.B.2.pdf',width=8,height=10,col=mycolors)
