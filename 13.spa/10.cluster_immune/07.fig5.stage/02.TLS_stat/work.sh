#perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$_};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/\t/);if(exists $ha{$a[0]}){print "$ha{$a[0]},$a[1]\n"}};for $i(keys %ha){@a=split(/,/,$ha{$i});if($a[1]=~/My/){print "$ha{$i},Control\n"}}' l2m.stat.merge.csv tls.list > l2m.stat.merge.csv.tls
perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$_};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/\t/);if(exists $ha{$a[0]}){print "$ha{$a[0]},$a[1]\n"}};for $i(keys %ha){@a=split(/,/,$ha{$i});if($a[1]=~/Myeloid7/ || $a[1]=~/Myeloid8/){print "$ha{$i},Control\n"}}' l2m.stat.merge.csv tls.list > l2m.stat.merge.csv.tls


perl -e 'while(<>){chomp;@a=split;print "perl get.ucell.pl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/04.plot/at/$a[0].immune.at /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/16.microenvironment/05.TLS/out/$a[1].$a[0].UCell.score.csv $a[0]\n"}' sn.list |sh > merge.all.ucell.list

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}="$a[1],$a[2]"};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/\,/);if(exists $ha{$a[0]}){print "$_,$ha{$a[0]}\n"}}' merge.all.ucell.list l2m.stat.merge.csv.tls > l2m.stat.merge.csv.tls.gc

cut -d "," -f 2,12 l2m.stat.merge.csv.tls.gc > input.csv

