perl -e 'while(<>){chomp;@a=split;print "cp /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/16.microenvironment/05.TLS/out/$a[1].$a[0].UCell.score.csv gc.out/$a[0].UCell.score.csv\ncp /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/05.ucell/02.ucell/out/$a[1].$a[0].UCell.score.csv fdc.out/$a[0].UCell.score.csv\n"}' sn.list |sh

perl -e 'while(<>){chomp;@a=split;print "perl filter.pl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/04.cluster_pla/at/$a[0].tls.at gc.out/$a[0].UCell.score.csv > gc.out/$a[0].UCell.score.csv.filter\n"}' sn.list|sh

perl -e 'while(<>){chomp;@a=split;print "perl filter.pl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/04.cluster_pla/at/$a[0].tls.at  fdc.out/$a[0].UCell.score.csv > fdc.out/$a[0].UCell.score.csv.filter\n"}' sn.list|sh

