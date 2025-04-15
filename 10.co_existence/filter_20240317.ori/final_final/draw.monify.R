library(pheatmap)

data <- read.table("./input.pearson.xls", header = T, row.names = 1, sep = "\t")

#mycolors <- colorRampPalette(c("#5E4EA1","#5267AD","#3F7EB8","#4596B7","#5AAEAD","#6CC4A4","#8ACFA4","#A5DAA4","#BFE4A0","#D7EE9B","#EAF69E","#F4FAAE","#FFFFBF","#FFF2A9","#FEE593","#FED380","#FEBE6E","#FCA85E","#F98E51","#F47245","#E95E47","#DC4A4C","#CB364C","#B42147","#9E0041"))(100)

pheatmap(data, scale = "none", fontsize_col = 10, fontsize_row = 10, border_color = "grey", cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "complete", file ="./output.pearson.complete.pdf", width = 9, height = 10,  breaks = seq(-0.8, 0.8, 0.02), treeheight_col = 50, treeheight_row = 0, colorbar_legend = list(at = seq(-0.8, 0.8, 0.2), labels = c("<= -0.8", "-0.4", "-0.2", "0", "0.2", "0.4", "0.8")))
pheatmap(data, scale = "none", fontsize_col = 10, fontsize_row = 10, border_color = "grey", cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "complete", file ="./output.pearson.complete.png", width = 9, height = 10,  breaks = seq(-0.8, 0.8, 0.02), treeheight_col = 50, treeheight_row = 0, colorbar_legend = list(at = seq(-0.8, 0.8, 0.2), labels = c("<= -0.8", "-0.4", "-0.2", "0", "0.2", "0.4", "0.8")))

