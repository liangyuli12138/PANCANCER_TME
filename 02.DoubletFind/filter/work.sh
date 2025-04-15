perl -e 'print "cell\n";while(<>){chomp;if(/^,/ || /T_like/ || /Doublet/){next};@a=split(/,/);print "$a[0]\n"}' merge.doublet.csv > merge.doublet.csv.input

perl -e 'print "cell,group\n";while(<>){chomp;if(/^,/ || /T_like/ || /Doublet/){next};@a=split(/,/);print "$a[0],$a[-5]\n"}' merge.doublet.csv > merge.doublet.csv.at

sed -i 's/SMC/Mural_cell/' merge.doublet.csv.at

sed -i 's/mastocyte/Mastocyte/' merge.doublet.csv.at

