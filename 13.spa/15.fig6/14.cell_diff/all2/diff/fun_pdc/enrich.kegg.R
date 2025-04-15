library(org.Hs.eg.db)
library(stringr)
library(clusterProfiler)
library(ggplot2)

setwd("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/14.cell_diff/all2/diff/fun_pdc")

diff <- read.table(file="all.for.enrich.gene.list",sep="\t",header=T)
gene_used <- split(diff$gene,diff$cluster)

gene_used <- sapply(gene_used,function(x){
  tmp <- bitr(x,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = org.Hs.eg.db)
  list(tmp$ENTREZID)
})

result_KEGG <- compareCluster(geneClusters = gene_used,fun = 'enrichKEGG',organism = 'hsa',pvalueCutoff = 0.2,qvalueCutoff = 0.2)
write.table(result_KEGG@compareClusterResult,file="result.enrich.KEGG.csv",quote = F,row.names = F,sep='\t')

