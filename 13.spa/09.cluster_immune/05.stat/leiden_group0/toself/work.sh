perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];<IN1>;print "cell,-log10\n";while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[0]}){print "$a[0],$ha{$a[0]}\n"}}' all.toself.stat.csv.log.at immune.group0.r0.5.obs > all.toself.stat.csv.log.group0.at

