perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=1};open IN1,$ARGV[1];<IN1>;while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[-3]}){print "$_\n"}}' filter.list pancancer.icar.all.cell.obs.ori > pancancer.icar.all.cell.obs.ori.filter

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[1]}=1};open IN1,$ARGV[1];print "cell,new_cluster,merge_groups,merge_celltype\n";while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[-1]}){print "$a[0],$a[-3],$a[-2],$a[-1]\n"}}' list pancancer.icar.all.cell.obs.ori.filter > pancancer.icar.all.cell.obs.immune.at


