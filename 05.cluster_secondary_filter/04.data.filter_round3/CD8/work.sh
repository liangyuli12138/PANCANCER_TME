perl -e 'print "cell\n";while(<>){chomp;@a=split(/,/);if($a[-1] eq "6" || $a[-1] eq "1" || $a[-1] eq "2" || $a[-1] eq "3" || $a[-1] eq "5"){print "$a[0]\n"}}' pancancer.ref.0728.final.CD8.umap.obs.txt > CD8.input

