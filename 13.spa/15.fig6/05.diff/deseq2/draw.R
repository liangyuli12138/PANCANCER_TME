library(pheatmap)
mycolors <- colorRampPalette(c("#5E4EA1","#5267AD","#3F7EB8","#4596B7","#5AAEAD","#6CC4A4","#8ACFA4","#A5DAA4","#BFE4A0","#D7EE9B","#EAF69E","#F4FAAE","#FFFFBF","#FFF2A9","#FEE593","#FED380","#FEBE6E","#FCA85E","#F98E51","#F47245","#E95E47","#DC4A4C","#CB364C","#B42147","#9E0041"))(25)
data = read.csv("./cluster_expr.csv",header = T,row.names=1);
#pheatmap(data,scale="column",fontsize_col=5,fontsize_row=5,border_color = "grey",cluster_rows=F,cluster_cols=F,file ='./test.pdf',width=8,height=10,col=mycolors,breaks=c(-2, -1, 0, 1, 2))

pheatmap(data,
         show_rownames = F,
         show_colnames = T,
         cluster_cols = F,
         cluster_rows=F,
         height=10,  
         breaks = c(seq(-2,-0.01,by=0.01),seq(0,2,by=0.01)),legend_breaks = -2:2,
         scale =  "row",
         frontsize = 10,
         angle_col=90, 
         # color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100),
         color =colorRampPalette(c("navy", "white", "firebrick3"))(400),
         clustering_method = 'single',
         file ='./test.png'
) 
