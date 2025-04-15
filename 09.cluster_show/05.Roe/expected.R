library(epitools)
dat <- read.csv("celltype.tissue.groups_thi.csv",header = TRUE, row.names=1, sep = "\t")
oi <- as.matrix(dat)
ei <- expected(oi)
write.csv(ei,"Roe.out")

