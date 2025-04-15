perl -e 'print ",marker\n";<>;while(<>){chomp;@a=split(/,/);$d=$a[0];$a[0]=~s/Cluster\d+//;print "$d,$a[0]\n"}' ../cluster.immune.cell.csv > cluster.immune.cell.csv.marker

perl -e '<>;while(<>){chomp;@a=split(/,/);$a[0]=~s/Cluster\d+//;$ha{$a[0]}{$a[1]}++;$hb{$a[1]}=1};for $i(sort {$a <=> $b} keys %hb){print ",group$i"};print "\n";for $j(sort {$a cmp $b} keys %ha){print "$j";for $i(sort {$a <=> $b} keys %hb){$ha{$j}{$i}=$ha{$j}{$i}?$ha{$j}{$i}:0;print ",$ha{$j}{$i}"};print "\n"}' ../immune.cluster.r0.5.obs > immune.cluster.r0.5.obs.stat.xls

