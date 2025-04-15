library(pheatmap)
mycolors <- colorRampPalette(c("#5E4EA1","#5267AD","#3F7EB8","#4596B7","#5AAEAD","#6CC4A4","#8ACFA4","#A5DAA4","#BFE4A0","#D7EE9B","#EAF69E","#F4FAAE","#FFFFBF","#FFF2A9","#FEE593","#FED380","#FEBE6E","#FCA85E","#F98E51","#F47245","#E95E47","#DC4A4C","#CB364C","#B42147","#9E0041"))(25)
data = read.csv("./correlation_df.xls",header = T,row.names=1);
pheatmap(data,scale="none",fontsize_col=4,fontsize_row=4,border_color = "grey",cluster_rows=T,cluster_cols=T,clustering_method="ward.D2",file ='./all.tls.pdf',width=24,height=30,col=mycolors)
row_labels_df <- data.frame(row_labels)
write.csv(row_labels_df, file = "cluster.sort.csv", row.names = FALSE)
