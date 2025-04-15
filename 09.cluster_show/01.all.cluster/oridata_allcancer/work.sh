perl -e '$t=<>;chomp $t;print "$t,groups_eig\n";while(<>){chomp;@a=split(/,/);if($a[1]=~/Epithelium_Malig/ || $a[1]=~/Lymphoid_Plamsa/){print "$_,$a[1]\n"}else{print "$_,$a[-1]\n"}}' merge.all.all.level.sev.at > merge.all.all.level.eig.at

