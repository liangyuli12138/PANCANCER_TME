library(org.Hs.eg.db)
library(stringr)
library(clusterProfiler)
library(ggplot2)

setwd("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/14.cell_diff/all2/diff/fun/enrich")

diff <- read.table(file="all.for.enrich.gene.list",sep="\t",header=T)
gene_used <- split(diff$gene,diff$cluster)
result_BP <- compareCluster(geneClusters = gene_used,fun = "enrichGO",ont="BP",keyType="SYMBOL",OrgDb=org.Hs.eg.db,pvalueCutoff = 0.1,qvalueCutoff = 0.1)
write.table(result_BP@compareClusterResult,file="result.enrich.GO.csv",quote = F,row.names = F,sep='\t')

