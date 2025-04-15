library(Seurat)
checklist <- read.csv("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/16.microenvironment/02.subtypes.all/scale/merge/hallmark.30.score.ori.csv",header=TRUE)
rownames(checklist) <- checklist$X
checklist <- checklist[, -1]
score <- checklist

path <- colnames(score)
for (i in 1:length(path)){
score[, i] <- scale(score[, i])
}

write.csv(score, file = "hallmark.30.score.scale.csv")

