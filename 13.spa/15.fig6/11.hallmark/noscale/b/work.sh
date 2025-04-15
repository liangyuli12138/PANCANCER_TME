perl get.score.pl pancancer.icar.background.cell.obs hallmark.30.score.csv.filter hallmark.30.score.ori.csv > merge.gene.signature.score.csv &

cut -d "," -f 1,3-13 merge.gene.signature.score.csv > merge.gene.signature.score.filter.csv

