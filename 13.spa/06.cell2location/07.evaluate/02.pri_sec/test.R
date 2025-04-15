library(ggplot2)


library(gghalves)

library(tidyverse)

data <- read.csv("data.csv",header=T)

从中随机选取10个基因来画图。

data <- data[sample(1:nrow(data), 10),]

还是以前的做法，长矩阵转化为短矩阵。

data_new <- melt(data,id="gene")

换了一下列名字

colnames(data_new) <- c("Genes","Samples","Values")

利用样品的ID来添加了样品的分组。

data_new$group <- str_split(data_new$Samples, "_",simplify = T)[,4]




ggplot(data_new, aes(x = Genes, y = Values, fill = group)) +

geom_violin(position = position_dodge(width = 1), scale = 'width') +

theme(axis.text.x = element_text(angle = 45, hjust = 1),

        legend.position = "top",

        legend.justification = "right")


pdf("/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/07.evaluate/02.pri_sec/test.pdf",w=16,h=8)

dev.off()

