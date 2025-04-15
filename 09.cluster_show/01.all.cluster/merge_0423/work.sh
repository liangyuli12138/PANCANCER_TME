perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];<IN1>;print "cell,celltype\n";while(<IN1>){chomp;@a=split(/,/);if($ha{$a[15]} eq "delete"){next};if(exists $ha{$a[15]}){print "$a[0],$ha{$a[15]}\n"}}' list pancancer.ref.0905.final.obs.csv > merge.celltype.20240423.csv

