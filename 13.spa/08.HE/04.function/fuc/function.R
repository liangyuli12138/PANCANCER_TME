data <- data.table::fread('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/04.function/out/fib.diff.list')
library(org.Hs.eg.db)
library(stringr)
library(clusterProfiler)
library(ggplot2)

gene_list <- split(data$gene,data$cluster)
gene_list$cluster<-NULL

#BP
GO_BP <- compareCluster(geneClusters = gene_list,fun = 'enrichGO',ont='BP',keyType='SYMBOL',OrgDb=org.Hs.eg.db,pvalueCutoff = 0.1,qvalueCutoff = 0.1)
write.table(GO_BP@compareClusterResult,file='/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/04.function/fuc/TLS.GO_BP_enrich_P0.1.txt',quote = F,row.names = F,sep='\t')
png('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/04.function/fuc/TLS.GO_BP_enrich_P0.1.png',w=12, h=12,units='in',res=300)
dotplot(GO_BP)+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+xlab('')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))
dev.off()

#CC
GO_CC <- compareCluster(geneClusters = gene_list,fun = 'enrichGO',ont='CC',keyType='SYMBOL',OrgDb=org.Hs.eg.db,pvalueCutoff = 0.1,qvalueCutoff = 0.1)
write.table(GO_CC@compareClusterResult,file='/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/04.function/fuc/TLS.GO_CC_enrich_P0.1.txt',quote = F,row.names = F,sep='\t')
png('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/04.function/fuc/TLS.GO_CC_enrich_P0.1.png',w=12, h=12,units='in',res=300)
dotplot(GO_CC)+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+xlab('')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))
dev.off()

#MF
GO_MF <- compareCluster(geneClusters = gene_list,fun = 'enrichGO',ont='MF',keyType='SYMBOL',OrgDb=org.Hs.eg.db,pvalueCutoff = 0.1,qvalueCutoff = 0.1)
write.table(GO_MF@compareClusterResult,file='/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/04.function/fuc/TLS.GO_MF_enrich_P0.1.txt',quote = F,row.names = F,sep='\t')
png('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/04.function/fuc/TLS.GO_MF_enrich_P0.1.png',w=12, h=12,units='in',res=300)
dotplot(GO_MF)+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+xlab('')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))
dev.off()

#KEGG
gene_used <- sapply(gene_list,function(x){
  tmp <- bitr(x,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
  list(tmp$ENTREZID)
})

module_KEGG <- compareCluster(geneClusters = gene_used,fun = 'enrichKEGG',organism = 'hsa',pvalueCutoff = 0.1,qvalueCutoff = 0.1)
write.table(module_KEGG@compareClusterResult,file='/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/04.function/fuc/TLS.KEGG_enrich_P0.1.txt',quote = F,row.names = F,sep='\t')
png('/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/04.function/fuc/TLS.KEGG_MF_enrich_P0.1.png',w=12, h=12,units='in',res=300)
dotplot(module_KEGG)+theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1))+xlab('')+
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))
dev.off()

