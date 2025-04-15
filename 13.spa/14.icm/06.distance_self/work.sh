perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$a[3]};open IN1,$ARGV[1];<IN1>;while(<IN1>){chomp;@a=split(/,/);if($a[2] eq $a[3]){print "$_,$ha{$a[0]}\n"}}' immune.cluster.11.r1.5.new.obs merge.csv > Lymphoid123.distance.stat.csv

