perl -e 'print "cell\n";open IN0,$ARGV[0];<IN0>;while(<IN0>){chomp;@a=split(/,/);if($a[-1]==3){next}else{print "$a[0]\n"}}open IN1,$ARGV[1];<IN1>;while(<IN1>){chomp;@a=split(/,/);if($a[-1]==3 || $a[-1]==4){print "$a[0]\n"}}' CD8.obs C1.obs > CD8.input

