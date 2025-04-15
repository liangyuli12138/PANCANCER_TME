perl -e '<>;while(<>){chomp;@a=split(/,/);$ha{$a[-3]}{$a[-4]}{$a[-12]}++;$hb{$a[-12]}=1};for $i(sort {$a cmp $b} keys %hb){print ",$i"};print "\n";for $i(sort {$a cmp $b} keys %ha){for $j(sort {$a cmp $b} keys %{$ha{$i}}){print "$i\t$j";for $k(sort {$a cmp $b} keys %hb){print "\t$ha{$i}{$j}{$k}"};print "\n"}}' pancancer.icar.contour.cell.obs > pancancer.icar.contour.cell.stat.xls

