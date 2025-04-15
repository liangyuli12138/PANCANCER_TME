perl -e '$t=<>;while(<>){chomp;@a=split;if($a[3] ne "Unknown" && $a[4] eq "Immune"){print "$_\n"}}' merge.stat.csv > merge.stat.csv.filter

perl merge.stat.pl immune.list score.list ori.list > merge.stat.csv

perl -e '$t=<>;print "$t";while(<>){chomp;@a=split;if($a[4] eq "others" && $a[3] ne "Lymphoid_Plamsa_IGLC" && $a[3] ne "Lymphoid_Plamsa_IGKC"){}else{print "$_\n"}}' merge.stat.csv.xls > merge.all.stat.csv.xls

cut -f 6-11 merge.stat2.csv |perl -e 'while(<>){chomp;@a=split;if(@a==6){print "$_\n"}}' > for.excel.xls

