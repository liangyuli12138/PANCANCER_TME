
perl -e '$t=<>;chomp $t;print "$t,immune_type\n";while(<>){chomp;@a=split(/,/);if($a[-3]<-250){$c="Intratumoral"}elsif($a[-3]>250){$c="Peritumoral"}else{$c="Boundary"};print "$_,$c\n"}' all.cell.dist.csv.type > all.cell.dist.csv.type.type

perl -e '<>;while(<>){chomp;@a=split(/,/);$ha{$a[1]}{$a[7]}++;if($a[3]=~/ICM/){$hb{$a[1]}{$a[7]}{$a[3]}++;$hc{$a[7]}{$a[3]}=1}};for $i(sort {$a cmp $b} keys %hc){for $j(sort {$a cmp $b} keys %{$hc{$i}}){print ",$i-$j"}};print "\n";for $k(sort {$a cmp $b} keys %hb){print "$k";for $i(sort {$a cmp $b} keys %hc){for $j(sort {$a cmp $b} keys %{$hc{$i}}){$o=$hb{$k}{$i}{$j}/$ha{$k}{$i};$o=sprintf "%.4f",$o;print ",$o"}};print "\n"}' all.cell.dist.csv.type.type > all.cell.dist.csv.type.type.icm.stat.xls

