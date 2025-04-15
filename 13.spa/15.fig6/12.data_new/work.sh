perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[1]}=1};open IN1,$ARGV[1];while(<IN1>){chomp;if(/^\,/){print "$_\n"};@a=split(/,/);if(exists $ha{$a[-3]}){print "$_\n"}}' filter.all.cluster.list pancancer.icar.all.cell.filter.obs > pancancer.icar.all.cell.filter.obs.filter.at

