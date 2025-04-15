perl -e 'print "cell\n";open IN0,$ARGV[0];<IN0>;while(<IN0>){chomp;@a=split(/,/);if($a[-1] eq "0" || $a[-1]==1 || $a[-1]==2 || $a[-1]==3 || $a[-1]==4 || $a[-1]==5 || $a[-1]==6 || $a[-1]==7 || $a[-1]==8 || $a[-1]==11 ){print "$a[0]\n"}}' ec.obs > EC.input

