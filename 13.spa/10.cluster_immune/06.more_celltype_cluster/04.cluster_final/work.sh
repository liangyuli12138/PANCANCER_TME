perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[1]}=$a[0]};open IN1,$ARGV[1];print "cell,new_groups\n";<IN1>;while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[0]}){print "$a[0],$ha{$a[0]}\n"}}' new.group.list immune.cluster.11.r1.5.obs > new.group.list.at

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split;$a[1]=~s/Cluster\d+//;print "$_\t$ha{$a[1]}\n"}' sn.list new.group.list > new.group.list.tissue

