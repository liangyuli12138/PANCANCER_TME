perl get.stat.pl obs.list bin100.merge.sub.percent.csv bin100.merge.sub.obs.csv bin100.merge.big.percent.csv bin100.merge.big.obs.csv

perl -e '<>;while(<>){chomp;@a=split(/,/);$ha{$a[1]}++};for $i(sort {$a <=> $b} keys %ha){print "$i\t$ha{$i}\n"}' bin100.merge.sub.obs.csv > bin100.merge.sub.obs.csv.stat

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);if($a[1]>=5){$ha{$a[0]}=1}};open IN1,$ARGV[1];<IN1>;while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[0]}){print "$_\n"}}' bin100.merge.sub.obs.csv bin100.merge.sub.percent.csv > bin100.merge.sub.percent.csv.filter

