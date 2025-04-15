perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);if($a[-1] eq "Doublet"){$ha{$a[0]}=1}};open IN1,$ARGV[1];print "cell,at\n";<IN1>;while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[0]}){print "$a[0],Doublet\n"}else{print "$a[0],Singlet\n"}}' merge.doublet.csv sub_cluster_EC.2000.30.obs.csv > doublet.at

