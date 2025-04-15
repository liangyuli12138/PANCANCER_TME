perl -e '<>;while(<>){chomp;@a=split(/,/);$ha{$a[9]}=1;$hb{$a[17]}=1;$hc{$a[9]}{$a[17]}++};for $i(sort {$a cmp $b} keys %hb){print "\t$i"};print "\n";for $j(sort {$a cmp $b} keys %ha){print "$j";for $i(sort {$a cmp $b} keys %hb){print "\t$hc{$j}{$i}"};print "\n"}' pancancer.ref.0905.final.obs.csv > celltype.tissue.groups_thi.csv

/hwfssz4/BC_PUB/Software/07.User-defined/03.Animal_Plant/wubin/mamba/bin/Rscript expected.R
