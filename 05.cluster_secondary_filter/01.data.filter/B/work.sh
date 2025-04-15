perl -e 'print "cell\n";open IN0,$ARGV[0];<IN0>;while(<IN0>){chomp;@a=split(/,/);if($a[-1]==5 || $a[-1]==6 || $a[-1]==9 || $a[-1]==10){next}else{print "$a[0]\n"}}' pancancer.ref.0723.final.B.umap.obs.txt > B.input

