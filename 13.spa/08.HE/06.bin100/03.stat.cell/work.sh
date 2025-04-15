perl -e 'open IN0,$ARGV[0];<IN0>;while(<IN0>){chomp;@a=split(/,/);$ha{$a[1]}++;$hb{$a[0]}=$a[1]};open IN1,$ARGV[1];$t=<IN1>;print "$t";while(<IN1>){chomp;@a=split(/,/);if($ha{$hb{$a[0]}}>100){print "$_\n"}}' bin200.r0.1.obs bin100.merge.sub.percent.csv.filter > bin100.merge.sub.percent.csv.filter.filter

perl get.cell.pl bin200.r0.1.obs bin100.merge.sub.percent.csv.filter.filter > bin100.merge.sub.percent.csv.filter.filter.matrix

perl get.excel.pl bin200.r0.1.obs bin100.merge.sub.percent.csv.filter.filter > bin100.merge.sub.percent.csv.filter.filter.xls

