perl -e '<>;while(<>){chomp;@a=split(/,/);$ha{$a[1]}++};for $i(sort {$a <=> $b} keys %ha){print "$i\t$ha{$i}\n"}' bin200.r0.3.obs |les

perl -e 'open IN0,$ARGV[0];<IN0>;while(<IN0>){chomp;@a=split(/,/);$ha{$a[1]}++;$hb{$a[0]}=$a[1]};open IN1,$ARGV[1];$t=<IN1>;print "$t";while(<IN1>){chomp;@a=split(/,/);if($ha{$hb{$a[0]}}>100){print "$_\n"}}' bin200.r0.3.obs bin100.merge.sub.percent.csv > bin100.merge.sub.percent.csv.filter

