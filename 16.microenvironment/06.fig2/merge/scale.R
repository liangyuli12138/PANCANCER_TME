library(Seurat)
checklist <- read.csv("merge.ucell.score.ori.csv",header=TRUE)
rownames(checklist) <- checklist$X
checklist <- checklist[, -1]
score <- checklist

path <- colnames(score)
for (i in 1:length(path)){
score[, i] <- scale(score[, i])
}

write.csv(score, file = "merge.ucell.score.scale.csv")

