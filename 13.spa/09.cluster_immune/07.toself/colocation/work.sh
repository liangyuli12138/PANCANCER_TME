perl -e 'while(<>){chomp;@a=split;$x=-log($a[2])/log(10);$x=$x>50?50:$x;print "$x\n"}' all.toself.stat.csv > all.toself.stat.csv.log

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;s/\_Cluster/Cluster/;@a=split(/,/);$ha{$a[0]}=$a[1]};print "cell,log\n";open IN1,$ARGV[1];<IN1>;while(<IN1>){chomp;@a=split(/,/);print "$a[0],$ha{$a[0]}\n"}' all.toself.stat.csv.log immune.cluster.r0.5.obs > all.toself.stat.csv.log.at

