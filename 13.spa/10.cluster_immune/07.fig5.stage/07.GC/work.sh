perl -e 'open IN0,$ARGV[0];while(<IN0>){s/_Cluster/Cluster/;chomp;@a=split;$ha{$a[0]}="$a[1],$a[2]"};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);print "$a[-1],$ha{$a[0]}\n"}' all.toself.stat.csv l2m.stat.merge.csv.sn.stage > input.csv

perl -e 'open IN0,$ARGV[0];while(<IN0>){s/_Cluster/Cluster/;chomp;@a=split(/,/);$ha{$a[0]}="$a[-1]"};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);print "$a[-1],$ha{$a[0]}\n"}' l2m.stat.merge.csv.tls.gc l2m.stat.merge.csv.sn.stage > input2.csv 

perl -e 'open IN0,$ARGV[0];<IN0>;while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$a[-1]};open IN1,$ARGV[1];$t=<IN1>;print "cell,GC_ucell\n";while(<IN1>){chomp;@a=split(/,/);print "$a[0],$ha{$a[0]}\n"}' merge.all.ucell.list immune.cluster.11.r1.5.new.obs > merge.all.ucell.list.at
