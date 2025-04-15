data <- data.table::fread('D:\\project\\pancancer\\8-15\\fib.diff.list1')
library(org.Hs.eg.db)
library(stringr)
library(clusterProfiler)
library(ggplot2)

gene_list <- split(data$gene,data$cluster)
gene_list$cluster<-NULL

#BP
GO_BP <- compareCluster(geneClusters = gene_list,fun = 'enrichGO',ont='BP',keyType='SYMBOL',OrgDb=org.Hs.eg.db,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
write.table(GO_BP@compareClusterResult,file='D:\\project\\pancancer\\8-15\\fib_GO_BP_enrich_P0.05.txt',quote = F,row.names = F,sep='\t')
png('D:\\project\\pancancer\\8-15\\fib_GO_BP_enrich_P0.05.png',w=12, h=12,units='in',res=300)
dotplot(GO_BP)+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+xlab('')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))
dev.off()

#CC
GO_CC <- compareCluster(geneClusters = gene_list,fun = 'enrichGO',ont='CC',keyType='SYMBOL',OrgDb=org.Hs.eg.db,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
write.table(GO_CC@compareClusterResult,file='D:\\project\\pancancer\\8-15\\fib_GO_CC_enrich_P0.05.txt',quote = F,row.names = F,sep='\t')
png('D:\\project\\pancancer\\8-15\\fib_GO_CC_enrich_P0.05.png',w=12, h=12,units='in',res=300)
dotplot(GO_CC)+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+xlab('')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))
dev.off()

#MF
GO_MF <- compareCluster(geneClusters = gene_list,fun = 'enrichGO',ont='MF',keyType='SYMBOL',OrgDb=org.Hs.eg.db,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
write.table(GO_MF@compareClusterResult,file='D:\\project\\pancancer\\8-15\\fib_GO_MF_enrich_P0.05.txt',quote = F,row.names = F,sep='\t')
png('D:\\project\\pancancer\\8-15\\fib_GO_MF_enrich_P0.05.png',w=12, h=12,units='in',res=300)
dotplot(GO_MF)+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+xlab('')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=100))
dev.off()

#ALL GO
GO_ALL <- compareCluster(geneClusters = gene_list,fun = 'enrichGO',ont='ALL',keyType='SYMBOL',OrgDb=org.Hs.eg.db,pvalueCutoff = 0.05,qvalueCutoff = 0.05)
#write.csv(module_GO@compareClusterResult,file='D:\\project\\pancancer\\module_analysis_6_15\\modules_filter_GO_enrich_P0.2.csv',quote = F,row.names = F)
write.table(GO_ALL@compareClusterResult,file='D:\\project\\pancancer\\8-15\\fib_GO_enrich_P0.05.txt',quote = F,row.names = F,sep='\t')

png('D:\\project\\pancancer\\8-15\\fib_GO_enrich_P0.05.png',w=12, h=12,units='in',res=300)
dotplot(GO_ALL)+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+xlab('')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=100))
dev.off()

#KEGG
gene_used <- sapply(gene_list,function(x){
  tmp <- bitr(x,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
  list(tmp$ENTREZID)
})

module_KEGG <- compareCluster(geneClusters = gene_used,fun = 'enrichKEGG',organism = 'hsa',pvalueCutoff = 0.05,qvalueCutoff = 0.05)
write.table(module_KEGG@compareClusterResult,file='D:\\project\\pancancer\\8-15\\fib_KEGG_enrich_P0.05.txt',quote = F,row.names = F,sep='\t')
png('D:\\project\\pancancer\\8-15\\fib_KEGG_enrich_P0.05.png',w=10, h=9,units='in',res=300)
dotplot(module_KEGG)+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+xlab('')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))
dev.off()
