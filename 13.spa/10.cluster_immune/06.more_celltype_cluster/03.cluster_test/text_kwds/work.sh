perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];print "cell,TLS,colour\n";<IN1>;while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[0]}){print "$a[0],$a[0]-TLS,$ha{$a[0]}\n"}else{print "$a[0],NA,black\n"}}' list immune.cluster.11.r1.5.obs > list.at 

