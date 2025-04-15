perl -e 'print "cell,difftype\n";<>;while(<>){chomp;@a=split(/\,/);if($a[-1]=~/Lymphoid_B/){if($a[-5]=~/Lymphoid1/ || $a[-5]=~/Lymphoid2/ || $a[-5]=~/Lymphoid3/){print "$a[0],Lymphoid_B.$a[-5]\n"}elsif($a[-5]=~/Mye/){print "$a[0],Lymphoid_B.Myeloid\n"}else{print "$a[0],other\n"}}else{print "$a[0],other\n"}}' pancancer.icar.all.cell.obs > pancancer.icar.all.cell.obs.at

