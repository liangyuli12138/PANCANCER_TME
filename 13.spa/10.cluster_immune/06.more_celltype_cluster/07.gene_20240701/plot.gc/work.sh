perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);print "$_,$ha{$a[0]}\n"}' l2m.stat.merge.csv.xls merge.all.ucell.list > merge.all.ucell.list.csv

