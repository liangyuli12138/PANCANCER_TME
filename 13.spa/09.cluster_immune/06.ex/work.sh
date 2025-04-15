perl merge.stat.pl immune.list score.list ori.list > merge.stat.ex.csv

perl -e '$t=<>;while(<>){chomp;@a=split;if($a[3] ne "Unknown" && $a[4] eq "Immune"){print "$_\n"}}' merge.stat.ex.csv > merge.stat.ex.csv.filter

perl -e 'open IN0,$ARGV[0];<IN0>;while(<IN0>){chomp;@a=split;@a=split;$id=$a[0]."_".$a[1];$ha{$id}=1};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split;@a=split;$id=$a[0]."_".$a[1];if(!exists $ha{$id}){print "$_\n"}}' ../05.stat/merge.stat.csv merge.stat.ex.csv > merge.stat.add.csv

