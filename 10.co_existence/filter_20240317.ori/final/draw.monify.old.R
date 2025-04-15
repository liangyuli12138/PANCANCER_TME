library(pheatmap)

data <- read.table("./input.pearson.xls", header = T, row.names = 1, sep = "\t")

highlight_mat <- matrix(NA, nrow = length(rownames(data)), ncol = 1)
rownames(highlight_mat) <- rownames(data)
colnames(highlight_mat) <- "Highlight"

highlight_mat[rownames(data) %in% c("Lymphoid_CD4_Tn","Lymphoid_CD8_Tn","Lymphoid_B_naive","Lymphoid_B_memory","Lymphoid_CD8_Tm","Myeloid_Marco_LYVE1","Lymphoid_CD4_Tcm","Lymphoid_NKT","Myeloid_cDC2","Lymphoid_CD8_Tisg","Myeloid_Mono"), ] <- 1
highlight_mat[rownames(data) %in% c("Myeloid_Marco_SPP1","Lymphoid_CD4_Treg","Lymphoid_CD8_Tex","Lymphoid_NK_CD56+","Lymphoid_CD4_Tfh","Myeloid_cDC3"), ] <- 2
highlight_mat[rownames(data) %in% c("Lymphoid_CD4_CTL","Lymphoid_NK_CD16+","Lymphoid_CD4_Th17","Lymphoid_ILC","Myeloid_pDC","Myeloid_cDC1","Myeloid_Marco_C1QC","Lymphoid_Plamsa_IGLC","Lymphoid_Plamsa_IGKC","Lymphoid_CD4_Tstr","Lymphoid_CD8_Tstr","Lymphoid_CD8_Teff","Lymphoid_MAIT"), ] <- 3

annotation_colors = list(highlight_mat = setNames(rep("#800080", nrow(highlight_mat)), rownames(highlight_mat)))

mycolors <- colorRampPalette(c("#5E4EA1","#5267AD","#3F7EB8","#4596B7","#5AAEAD","#6CC4A4","#8ACFA4","#A5DAA4","#BFE4A0","#D7EE9B","#EAF69E","#F4FAAE","#FFFFBF","#FFF2A9","#FEE593","#FED380","#FEBE6E","#FCA85E","#F98E51","#F47245","#E95E47","#DC4A4C","#CB364C","#B42147","#9E0041"))(100)

pheatmap(data, scale = "none", fontsize_col = 10, fontsize_row = 10, border_color = "grey", cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "complete", file ="./output.pearson.complete.pdf", width = 10, height = 10, col = mycolors, breaks = seq(-0.8, 0.8, 0.02), treeheight_col = 50, treeheight_row = 0, annotation_row = highlight_mat, annotation_colors = list(highlight_mat = c("1" = "#800080", "2" = "#800080", "3" = "#800080")), annotation_names_row = TRUE)

pheatmap(data, scale = "none", fontsize_col = 10, fontsize_row = 10, border_color = "grey", cluster_rows = TRUE, cluster_cols = TRUE, clustering_method = "complete", file ="./output.pearson.complete.png", width = 10, height = 10, col = mycolors, breaks = seq(-0.8, 0.8, 0.02), treeheight_col = 50, treeheight_row = 0, annotation_row = highlight_mat, annotation_colors = list(highlight_mat = c("1" = "#800080", "2" = "#800080", "3" = "#800080")), annotation_names_row = TRUE)

