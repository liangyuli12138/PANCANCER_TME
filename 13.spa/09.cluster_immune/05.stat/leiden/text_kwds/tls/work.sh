perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$x=$a[0]."Cluster".$a[1];$ha{$x}=1};open IN1,$ARGV[1];print "cell,TLS\n";<IN1>;while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[0]}){print "$a[0],$a[0]-TLS\n"}else{print "$a[0],NA\n"}}' list immune.cluster.r0.5.obs > list.at

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];print "cell,TLS,colour\n";<IN1>;while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[0]}){print "$a[0],$a[0]-TLS,$ha{$a[0]}\n"}else{print "$a[0],NA,black\n"}}' list immune.cluster.r0.5.obs > list.at

