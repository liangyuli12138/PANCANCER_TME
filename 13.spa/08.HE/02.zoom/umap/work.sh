perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$a[-3]};print "\"\",\"id\",\"celltype\"\n";open IN1,$ARGV[1];<IN1>;while(<IN1>){chomp;$n++;@a=split(/,/);@z=split(/\_/,$a[0]);if($z[0]>16500 && $z[0]<20300 && $z[1]>9100 && $z[1]<13200 && /immune/){print "\"$n\",\"$a[0]\",\"$ha{$a[0]}\"\n"}elsif($z[0]>15917 && $z[0]<20883 && $z[1]>8850 && $z[1]<13636){print "\"$n\",\"$a[0]\",\"others\"\n"}}' D01972D1_cellbin.final.celltype.obs.csv D01972D1_region.csv > D01972D1_region.csv.at

perl -e '<>;while(<>){chomp;@a=split(/,/);$ha{$a[-1]}++};for $i(keys %ha){print "$i\t$ha{$i}\n"}' D01972D1_region.csv.at |sort -k 2,2nr > D01972D1_region.csv.at.stat

#perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/\t/);$ha{$a[0]}=1};print "\"\",\"id\",\"celltype\"\n";open IN1,$ARGV[1];<IN1>;while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[2]}){print "$_\n"}}' D01972D1_region.csv.at.stat.filter D01972D1_region.csv.at > D01972D1_region.csv.at.filter

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/\t/);$ha{$a[0]}=$a[2]};print "\"\",\"id\",\"celltype\"\n";open IN1,$ARGV[1];<IN1>;while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[2]}){print "$a[0],$a[1],\"$ha{$a[2]}\"\n"}}' D01972D1_region.csv.at.stat.filter D01972D1_region.csv.at > D01972D1_region.csv.at.filter

