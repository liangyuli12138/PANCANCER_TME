perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=1};open IN1,$ARGV[1];<IN1>;print "cell\n";while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[-3]}){print "$a[0]\n"}}' filter.list pancancer.icar.all.cell.obs > pancancer.icar.all.cell.obs.input


