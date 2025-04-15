perl -e '$t=<>;chomp $t;print "$t,stage\n";while(<>){chomp;@a=split(/,/);if($a[1] eq "Lymphoid1" || $a[1] eq "Lymphoid2"){print "$_,early\n"}elsif($a[1] eq "Lymphoid3"){print "$_,late\n"}}' l2m.stat.merge.csv.sn > l2m.stat.merge.csv.sn.stage

perl -e '$t=<>;chomp $t;print "$t,stage\n";while(<>){chomp;@a=split(/,/);if($a[1] eq "Lymphoid0" || "Lymphoid1" || $a[1] eq "Lymphoid2"){print "$_,early\n"}elsif($a[1] eq "Lymphoid3"){print "$_,late\n"}}' l2m.stat.merge.csv.sn > l2m.stat.merge.csv.sn.stage
