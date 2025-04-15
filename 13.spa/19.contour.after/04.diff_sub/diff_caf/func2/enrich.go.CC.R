library(org.Hs.eg.db)
library(stringr)
library(clusterProfiler)
library(ggplot2)

setwd("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/19.contour.after/04.diff_sub/diff_caf/func2")

diff <- read.table(file="all.for.enrich.gene.list",sep="\t",header=T)
gene_used <- split(diff$gene,diff$cluster)
result_BP <- compareCluster(geneClusters = gene_used,fun = "enrichGO",ont="CC",keyType="SYMBOL",OrgDb=org.Hs.eg.db,pvalueCutoff = 0.1,qvalueCutoff = 0.1)
write.table(result_BP@compareClusterResult,file="result.enrich.GO.cc.csv",quote = F,row.names = F,sep='\t')

