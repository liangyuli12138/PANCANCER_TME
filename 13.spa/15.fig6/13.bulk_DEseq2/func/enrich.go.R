library(org.Hs.eg.db)
library(stringr)
library(clusterProfiler)
library(ggplot2)

setwd("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/13.bulk_DEseq2/func")

diff <- read.table(file="all.for.enrich.gene.list",sep="\t",header=T)
gene_used <- split(diff$gene,diff$cluster)

result_BP <- compareCluster(geneClusters = gene_used,fun = "enrichGO",ont="BP",keyType="SYMBOL",OrgDb=org.Hs.eg.db,pvalueCutoff = 0.1,qvalueCutoff = 0.1)
write.table(result_BP@compareClusterResult,file="result.enrich.GO.BP.csv",quote = F,row.names = F,sep='\t')

result_CC <- compareCluster(geneClusters = gene_used,fun = "enrichGO",ont="CC",keyType="SYMBOL",OrgDb=org.Hs.eg.db,pvalueCutoff = 0.1,qvalueCutoff = 0.1)
write.table(result_CC@compareClusterResult,file="result.enrich.GO.CC.csv",quote = F,row.names = F,sep='\t')

result_MF <- compareCluster(geneClusters = gene_used,fun = "enrichGO",ont="MF",keyType="SYMBOL",OrgDb=org.Hs.eg.db,pvalueCutoff = 0.1,qvalueCutoff = 0.1)
write.table(result_MF@compareClusterResult,file="result.enrich.GO.MF.csv",quote = F,row.names = F,sep='\t')

