perl -e 'print "cell\n";open IN0,$ARGV[0];<IN0>;while(<IN0>){chomp;@a=split(/,/);if($a[-1]==10 || $a[-1]==11 || $a[-1]==15){next}else{print "$a[0]\n"}}' my.obs > Myeloid.input

