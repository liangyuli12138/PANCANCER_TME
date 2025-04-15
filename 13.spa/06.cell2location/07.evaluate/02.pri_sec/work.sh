cat list/*csv.score > all.csv.score &

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;s/\"//g;@a=split;if(!exists $ha{$a[0]}){next};$hb{$ha{$a[0]}}=1;$hx{$ha{$a[0]}}++;$hy{$hx{$ha{$a[0]}}}{$ha{$a[0]}}=$a[1];$hz{$hx{$ha{$a[0]}}}{$ha{$a[0]}}=$a[2]};for $i(sort {$a cmp $b} keys %hb){print "\t$i","_pri\t$i","_sec"};print "\n";for $j(sort {$a<=>$b} keys %hy){for $i(sort {$a cmp $b} keys %hb){print "\t$hy{$j}{$i}\t$hz{$j}{$i}"};print "\n"}' cluster.list all.csv.score > all.cell2location.score.csv.xls

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};print "cell,value,group\n";open IN1,$ARGV[1];while(<IN1>){chomp;s/\"//g;@a=split;if(!exists $ha{$a[0]}){next};$x=$a[1]/($a[1]+$a[2]);$y=$a[2]/($a[1]+$a[2]);print "$a[0],$x,primary\n$a[0],$y,secondary\n"}' cluster.list all.csv.score > all.cell2location.score.csv.xls &

