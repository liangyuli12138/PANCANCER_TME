perl -e '<>;while(<>){chomp;@a=split(/,/);if($a[8]=~/\_N/){next};$ha{$a[8]}{$a[-6]}++;$ht{$a[8]}++;$hr{$a[-6]}=1};for $i(keys %ha){for $j(keys %{$ha{$i}}){$hb{$j}{$i}=$ha{$i}{$j}}};for $i(sort {$a cmp $b} keys %ht){print "\t$i"};print "\n";for $j(sort {$a cmp $b} keys %hr){print "$j";for $i(sort {$a cmp $b} keys %ht){$hb{$j}{$i}=$hb{$j}{$i}?$hb{$j}{$i}:0;print "\t$hb{$j}{$i}"};print "\n"}' pancancer.ref.0905.final.obs.csv > pancancer.ref.0905.final.obs.csv.to.windows.xls

/hwfssz4/BC_PUB/Software/03.Soft_ALL/R-3.5.1/bin/Rscript draw.R

/hwfssz4/BC_PUB/Software/03.Soft_ALL/R-3.5.1/bin/Rscript draw.origin.R 
