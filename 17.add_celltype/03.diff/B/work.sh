ls */*|perl -e 'while(<>){chomp;$t=$_;print "$t\n";open IN,$t;while(<IN>){chomp;if(/,FDCSP,/ || /,FCER2,/ || /,CR1,/ || /,CR2,/ || /,CD35,/ || /MFGE8/ || /VCAM1/ || /CXCL13/ || /CNA42/){print "$_\n"}}}' > stat.b.csv

ls */*|perl -e 'while(<>){chomp;$t=$_;print "$t\n";open IN,$t;while(<IN>){chomp;@a=split;if($a[0]>5000){last};if(/,FDCSP,/ || /,FCER2,/ || /,CR1,/ || /,CR2,/ || /,CD35,/ || /MFGE8/ || /VCAM1/ || /CXCL13/ || /CNA42/){print "$_\n"}};print "\n"}' > stat.b.csv &

ls */*|perl -e 'while(<>){chomp;$t=$_;$o= ">$t\n";open IN,$t;$n=0;while(<IN>){chomp;@a=split;if($a[0]>1000){last};if(/,FDCSP,/ || /,FCER2,/ || /,CR1,/ || /,CR2,/ || /,CD35,/ || /MFGE8/ || /VCAM1/ || /CXCL13/ || /CNA42/){$n++; $o .= "$_\n"}};if($n>=3){print "$o\n"}}' > stat.b.csv &

ls 2000_50_1/* 2000_50_1.2/* 2000_50_1.5/* 2000_50_2/*|perl -e 'while(<>){chomp;$t=$_;$o= ">$t\n";open IN,$t;$n=0;while(<IN>){chomp;@a=split;if($a[0]>5000){last};if(/,FDCSP,/ || /,FCER2,/ || /,CR1,/ || /,CR2,/ || /,CD35,/ || /MFGE8/ || /VCAM1/ || /CXCL13/ || /CNA42/){$n++; $o .= "$_\n"}};if($n>=3){print "$o\n"}}' > stat.b.csv &

