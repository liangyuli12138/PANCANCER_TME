perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[2]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);print "$ha{$a[0]},$_\n"}' sn.list.xls l2m.stat.merge.csv.xls > l2m.stat.merge.csv.tissue

