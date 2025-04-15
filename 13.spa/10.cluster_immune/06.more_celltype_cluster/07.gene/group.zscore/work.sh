perl -e '<>;print "cell,new_groups\n";while(<>){chomp;@a=split(/,/);if(!exists $ha{$a[6]}){$ha{$a[6]}=$a[4]};};for $i(sort {$a cmp $b} keys %ha){print "$i,$ha{$i}\n"}' merge_25chip_immune_area.nor.obs > cluster2group.csv


