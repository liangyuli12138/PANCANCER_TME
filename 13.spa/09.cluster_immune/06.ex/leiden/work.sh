perl -e 'while(<>){chomp;@a=split;$id=$a[0]."_".$a[2];$ha{$id}{$a[3]}++;$hb{$a[3]}=1};for $i(sort {$a cmp $b} keys %hb){print ",$i"};print "\n";for $j(sort {$a cmp $b} keys %ha){print "$j";for $i(sort {$a cmp $b} keys %hb){$ha{$j}{$i}=$ha{$j}{$i}?$ha{$j}{$i}:0;print ",$ha{$j}{$i}"};print "\n"}' ../merge.stat.add.csv > merge.stat.add.csv.stat.xls

