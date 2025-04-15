cut -f 1 sn.list |while read i;do echo perl get.at.pl ../03.net/out/$i.output.txt /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/$i""_cellbin.final.celltype.obs.csv \> at/$i.immune.at;done|sh

cut -f 1 sn.list |while read i;do echo perl get.at2.pl ../03.net/out/$i.output.ex.txt /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/$i""_cellbin.final.celltype.obs.csv $i;done|sh

for i in  ../03.net/out/*.output.txt;do echo perl get.stat.pl $i;done|sh > stat.num.xls

perl get.group.pl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/01.Hexcell/region_merge/D01972D1_region.csv ../05.stat/leiden/immune.cluster.r0.5.obs at/D01972D1.immune.at |les

perl -e 'while(<>){chomp;@a=split;print "perl get.group.pl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/01.Hexcell/region_merge/$a[0]","_region.csv ../05.stat/leiden/immune.cluster.r0.5.obs at/$a[0].immune.at > at/$a[0].group.at\n"}' sn.list |sh &

perl -e 'while(<>){chomp;@a=split;print "perl get.group.pl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/01.Hexcell/region_merge/$a[0]","_region.csv ../05.stat/leiden2/immune.cluster.r0.5.obs at/$a[0].immune.at > at.pla/$a[0].group.at\n"}' sn.list |sh &

perl -e 'while(<>){chomp;@a=split;print "perl get.group.pl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/01.Hexcell/region_merge/$a[0]","_region.csv ../05.stat/leiden3/immune.cluster.r1.obs at/$a[0].immune.at > at.group.zscore.r1/$a[0].group.at\n"}' sn.list |sh &

perl -e 'while(<>){chomp;@a=split;print "perl get.ori.add.pl at.pla/$a[0].group.at at/$a[0].ori.at > at.pla/$a[0].ori.at\n"}' sn.list |sh

