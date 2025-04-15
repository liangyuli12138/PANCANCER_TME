########寻找moudle#########
library(igraph)
adj <- read.table('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\filter.list.yes.gene.stat.filter.matrix')#基因邻接矩阵
adj <- as.matrix(adj)
g_adj <- graph_from_adjacency_matrix(adj, diag = FALSE, mode = 'undirected', weighted = TRUE)
plot(g_adj,layout=layout_in_circle)
imc <- cluster_infomap(g_adj,nb.trials = 100,modularity=F) #寻找模块
modules <- communities(imc)
plot(imc,g_adj)
names(modules) <- paste0('module-',names(modules))
saveRDS(modules,'D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\modules_clutser_infomap.rds')

########模块GO富集注释#######
##模块过滤后富集注释
library(org.Hs.eg.db)
library(stringr)
library(clusterProfiler)
library(ggplot2)

module_gene_number <- sapply(modules,function(x) length(x))
module_filter <- modules[module_gene_number>9]
module_filter_gene <- sapply(names(module_filter),function(x){
  tmp <- data.frame(module_name=x,gene=paste0(module_filter[[x]],collapse=','))
  list(tmp)
})
module_filter_gene <- Reduce('rbind',module_filter_gene)
write.csv(module_filter_gene,file='D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\modules_filter_gene.csv',quote = F,row.names = F)
modules_filter_enrich <- sapply(module_filter,function(x) list(x))
module_GO <- compareCluster(geneClusters = modules_filter_enrich,fun = 'enrichGO',ont='ALL',keyType='SYMBOL',OrgDb=org.Hs.eg.db,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
#write.csv(module_GO@compareClusterResult,file='D:\\project\\pancancer\\module_analysis_6_15\\modules_filter_GO_enrich.csv',quote = F,row.names = F)
write.table(module_GO@compareClusterResult,file='D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\modules_filter_GO_enrich.txt',quote = F,row.names = F,sep='\t')

png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\modules_filter_GO_enrich.png',w=13, h=17,units='in',res=300)
dotplot(module_GO)+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+xlab('')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=100))
dev.off()

##模块过滤后GO富集注释(分MF、CC、BP)
#BP
module_GO_BP <- compareCluster(geneClusters = modules_filter_enrich,fun = 'enrichGO',ont='BP',keyType='SYMBOL',OrgDb=org.Hs.eg.db,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
write.table(module_GO_BP@compareClusterResult,file='D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\modules_filter_GO_BP_enrich.txt',quote = F,row.names = F,sep='\t')
png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\modules_filter_GO_BP_enrich.png',w=15, h=17,units='in',res=300)
dotplot(module_GO_BP)+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+xlab('')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))
dev.off()

#CC
module_GO_CC <- compareCluster(geneClusters = modules_filter_enrich,fun = 'enrichGO',ont='CC',keyType='SYMBOL',OrgDb=org.Hs.eg.db,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
write.table(module_GO_CC@compareClusterResult,file='D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\modules_filter_GO_CC_enrich.txt',quote = F,row.names = F,sep='\t')
png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\modules_filter_GO_CC_enrich.png',w=15, h=17,units='in',res=300)
dotplot(module_GO_CC)+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+xlab('')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))
dev.off()

#MF
module_GO_MF <- compareCluster(geneClusters = modules_filter_enrich,fun = 'enrichGO',ont='MF',keyType='SYMBOL',OrgDb=org.Hs.eg.db,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
write.table(module_GO_MF@compareClusterResult,file='D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\modules_filter_GO_MF_enrich.txt',quote = F,row.names = F,sep='\t')
png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\modules_filter_GO_MF_enrich.png',w=14, h=16,units='in',res=300)
dotplot(module_GO_MF)+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+xlab('')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))
dev.off()

#KEGG
gene_used <- sapply(modules_filter_enrich,function(x){
  tmp <- bitr(x,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
  list(tmp$ENTREZID)
})

module_KEGG <- compareCluster(geneClusters = gene_used,fun = 'enrichKEGG',organism = 'hsa',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
write.table(module_KEGG@compareClusterResult,file='D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\modules_filter_KEGG_enrich.txt',quote = F,row.names = F,sep='\t')
png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\modules_filter_KEGG_enrich.png',w=13, h=14,units='in',res=300)
dotplot(module_KEGG)+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+xlab('')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))
dev.off()

#########heatmap############
#模块热图
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(pheatmap))

data <- read.table('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\final.0616.malig.gene.modules.txt',header=T,row.names = 1)
data <- t(data)

ann_col<-data.frame(Cancer=colnames(data))
rownames(ann_col) <- ann_col$Cancer
ann_col$Cancer <-data.table::tstrsplit(ann_col$Cancer,'_')[[1]]
cols_set <- c('#F68282','#B95FBB','#78BE94','#ff9a36','#2FF18B','#B84D64','#faf4cf','#CCB1F1','#25aff5','#A4DFF2','#7CA878','#AC8F14',
              '#35A132','#8DD3C7','#CFECBB','#F4F4B9','#AF98B5','#8952A0','#F4867C','#C0979F','#86B1CD')
ann_col$col <- cols_set[as.numeric(as.factor(ann_col$Cancer))]
ann_col1 <- unique(ann_col)
ann_col2 <- ann_col1$col
names(ann_col2) <- ann_col1$Cancer
ha <- HeatmapAnnotation(Cancer=ann_col$Cancer,
                        col = list(Cancer=ann_col2))

png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\module_heatmap.png',res=300,units='in',h=8,w=8)
Heatmap(data, name="Overlap",column_names_rot=45,cluster_rows=F,cluster_columns = F,top_annotation = ha,show_column_names=F,
        row_names_gp=grid::gpar(fontsize=10),show_column_dend = F,show_row_dend = F,row_names_side = "left",
        col=c('#1760A9','#F2F6F5','#B6252F'))
dev.off()

cols_set1 <- rev(heat.colors(20))[c(-1,-2,-3,-4)]
#library("viridis")  
#viridis(10)
#cols_set2 <- c("#440154FF","#482878FF","#3E4A89FF","#31688EFF","#26828EFF","#1F9E89FF", "#35B779FF","#6DCD59FF", "#B4DE2CFF" ,"#FDE725FF")
cols_set2 <- c("#3E4A89FF","#31688EFF","#26828EFF","#1F9E89FF", "#35B779FF","#6DCD59FF", "#B4DE2CFF" ,"#FDE725FF")

pdf('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\module_heatmap.pdf',h=6,w=8)
Heatmap(data, name="Overlap",column_names_rot=45,cluster_rows=F,cluster_columns = F,top_annotation = ha,show_column_names=F,
        row_names_gp=grid::gpar(fontsize=10),show_column_dend = F,show_row_dend = F,row_names_side = "left",
        col=c(cols_set2,cols_set1),height = nrow(data)*unit(10, "mm"))
dev.off()

png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\module_heatmap.png',h=6,w=8,res=300,units='in')
Heatmap(data, name="Overlap",column_names_rot=45,cluster_rows=F,cluster_columns = F,top_annotation = ha,show_column_names=F,
        row_names_gp=grid::gpar(fontsize=10),show_column_dend = F,show_row_dend = F,row_names_side = "left",
        col=c(cols_set2,cols_set1),height = nrow(data)*unit(10, "mm"))
dev.off()

png('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\module_heatmap.png',h=5,w=8,res=300,units='in')
Heatmap(data, name="Overlap",column_names_rot=45,cluster_rows=F,cluster_columns = F,top_annotation = ha,show_column_names=F,
        row_names_gp=grid::gpar(fontsize=10),show_column_dend = F,show_row_dend = F,row_names_side = "left",
        col=c(cols_set2,cols_set1),height = nrow(data)*unit(8, "mm"))
dev.off()

pdf('D:\\project\\pancancer\\module-analysis\\module_analysis_6_16\\module_heatmap.pdf',h=5,w=8)
Heatmap(data, name="Overlap",column_names_rot=45,cluster_rows=F,cluster_columns = F,top_annotation = ha,show_column_names=F,
        row_names_gp=grid::gpar(fontsize=10),show_column_dend = F,show_row_dend = F,row_names_side = "left",
        col=c(cols_set2,cols_set1),height = nrow(data)*unit(8, "mm"))
dev.off()
