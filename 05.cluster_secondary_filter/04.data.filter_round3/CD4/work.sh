perl -e 'print "cell\n";while(<>){chomp;@a=split(/,/);if($a[-1] eq "0" || $a[-1] eq "1" || $a[-1] eq "2" || $a[-1] eq "3" || $a[-1] eq "4"){print "$a[0]\n"}}' pancancer.ref.0728.final.CD4.umap.obs.txt > CD4.input

