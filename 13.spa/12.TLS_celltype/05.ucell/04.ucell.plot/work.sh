perl -e 'while(<>){chomp;@a=split;print "perl get.ucell.pl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/04.plot/at/$a[0].immune.at /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/12.TLS_celltype/05.ucell/03.stat.ucell/fdc.out/$a[0].UCell.score.csv $a[0]\n"}' sn.list |sh > merge.all.ucell.list

perl -e 'open IN0,$ARGV[0];<IN0>;while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$a[-2]};open IN1,$ARGV[1];$t=<IN1>;print "cell,FDC_ucell\n";while(<IN1>){chomp;@a=split(/,/);print "$a[0],$ha{$a[0]}\n"}' merge.all.ucell.list immune.cluster.11.r1.5.new.obs > merge.all.ucell.list.at

perl -e 'open IN0,$ARGV[0];<IN0>;while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);if($ha{$a[0]}=~/Lym/){print "$_,$ha{$a[0]}\n"}}' l2m.stat.merge.csv merge.merge.all.ucell.list.at > merge.merge.all.ucell.list.at.input

