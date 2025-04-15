perl -e 'print "cell\n";open IN0,$ARGV[0];<IN0>;while(<IN0>){chomp;@a=split(/,/);if($a[-1]==4 || $a[-1]==9){print "$a[0]\n"}}' all.obs > NK.input

