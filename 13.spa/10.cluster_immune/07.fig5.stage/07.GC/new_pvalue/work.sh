perl -e 'open IN0,$ARGV[0];while(<IN0>){s/_Cluster/Cluster/;chomp;@a=split;$ha{$a[0]}="$a[1],$a[2]"};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);print "$a[1],$ha{$a[0]}\n"}' all.toself.stat.csv l2m.stat.merge.csv.sn.stage > input.csv

perl -e 'open IN0,$ARGV[0];while(<IN0>){s/_Cluster/Cluster/;chomp;@a=split;$ha{$a[0]}="$a[1],$a[2]"};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);print "$a[1],$ha{$a[0]}\n"}' all.toself.stat.csv l2m.stat.merge.csv.sn.stage|awk '/Lymphoid/' > input.csv

