perl -e 'while(<>){chomp;@a=split;print "perl merge.score.group.cell.pl  /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/08.ucell/out/$a[1].$a[0].UCell.score.csv /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/04.plot/at.pla/$a[0].group.at /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/04.plot/at.pla/$a[0].ori.at > at/$a[0].merge.at\n"}' sn.list |sh &

perl get.mean.group.pl sn.list immune.cluster.r0.5.obs > tls.development.gene.at

perl get.mean.group.B.pl sn.list immune.cluster.r0.5.obs > tls.development.gene.at.B

