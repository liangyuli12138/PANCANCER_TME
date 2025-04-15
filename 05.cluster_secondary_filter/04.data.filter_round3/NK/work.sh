perl -e 'print "cell\n";while(<>){chomp;@a=split(/,/);if($a[-1] eq "5"){print "$a[0]\n"}}' pancancer.ref.0728.final.CD4.umap.obs.txt > NK.mait.input
perl -e 'while(<>){chomp;@a=split(/,/);if($a[-1] eq "4"){print "$a[0]\n"}}' pancancer.ref.0728.final.CD8.umap.obs.txt >> NK.mait.input
perl -e '<>;while(<>){chomp;@a=split(/,/);print "$a[0]\n"}' pancancer.ref.0728.final.NK.umap.obs.txt >> NK.mait.input 


