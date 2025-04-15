perl -e '<>;while(<>){chomp;@a=split(/,/);if($a[-1]<-250){$c="Intratumoral"}elsif($a[-1]>250){$c="Peritumoral"}else{$c="Boundary"};$m="$c"."-"."$a[3]";$ha{$m}=1;$z="$a[0]_$a[1]";$hb{$z}{$m}++};for $i(sort {$a cmp $b} keys %ha){print ",$i"};print "\n";for $j(sort {$a cmp $b} keys %hb){print "$j";for $i(sort {$a cmp $b} keys %ha){$hb{$j}{$i}=$hb{$j}{$i}?$hb{$j}{$i}:0;print ",$hb{$j}{$i}"};print "\n"}' all.immune.dist.csv > all.immune.dist.stat.csv &

perl -e '$t=<>;chomp $t;print "$t,immune_type\n";while(<>){chomp;@a=split(/,/);if($a[-1]<-250){$c="Intratumoral"}elsif($a[-1]>250){$c="Peritumoral"}else{$c="Boundary"};print "$_,$c\n"}' all.all.dist.csv > all.all.dist.type.csv

