perl get.pl icm.list all_cell_meta.csv > all_cell_meta.csv.icm

perl -e '@c=(ICM1,ICM2,ICM3);open IN0,$ARGV[0];<IN0>;while(<IN0>){chomp;@a=split(/,/);if($a[-1] !~ /ICM/){next};$x=$a[1].$a[3];$ha{$x}++;$hb{$x}{$a[-1]}++};open IN1,$ARGV[1];$t=<IN1>;chomp $t;print "$t,ICM1_percent,ICM2_percent,ICM3_percent\n";while(<IN1>){chomp;@a=split(/,/);print "$_";for ($i=0;$i<3;$i++){if($ha{$a[0]}{$c[$i]}>0){$p=$hb{$c[$i]}/$ha{$a[0]}{$c[$i]}}else{$p=0};print ",$p"};print "\n"}' all_cell_meta.csv.icm l2m.stat.merge.csv.at |les

perl -e '@c=(ICM1,ICM2,ICM3);open IN0,$ARGV[0];<IN0>;while(<IN0>){chomp;@a=split(/,/);if($a[-1] !~ /ICM/){next};$x=$a[1].$a[3];$ha{$x}++;$hb{$x}{$a[-1]}++};open IN1,$ARGV[1];$t=<IN1>;chomp $t;print "$t,ICM1_percent,ICM2_percent,ICM3_percent\n";while(<IN1>){chomp;@a=split(/,/);print "$_";for ($i=0;$i<3;$i++){if($ha{$a[0]}>0){$p=$hb{$a[0]}{$c[$i]}/$ha{$a[0]}}else{$p=0};print ",$p"};print "\n"}' all_cell_meta.csv.icm l2m.stat.merge.csv.at > l2m.stat.merge.csv.at.icm

