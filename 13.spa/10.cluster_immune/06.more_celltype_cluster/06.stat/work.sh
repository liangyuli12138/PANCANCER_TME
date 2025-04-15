perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];<IN1>;while(<IN1>){chomp;@a=split(/,/);print "$a[0]\t$a[1]\t$ha{$a[1]}\n"}' sn.list cluster.immune.cell.csv.marker > sn.list.xls

