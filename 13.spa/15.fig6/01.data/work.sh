cut -f 1 sn.list |while read i;do echo perl get.group.pl new.group.list.at at100/$i.ex.at $i \> at100.group/$i.ex.at;done|sh &

perl -e 'while(<>){chomp;@a=split;$o=`cat tmp.py`;$o=~s/aaaa/$a[0]/g; print "$o"}' sn.list >> merge.py 

perl -e 'print "cell,new_cluster,merge_groups\n";<>;while(<>){chomp;@a=split(/,/);$b=$a[0];$b=~s/\_\d+\_\d+$//;if($a[-2] eq Lymphoid1 || $a[-2] eq Lymphoid2){print "$a[0],$b","$a[-3],Lymphoid_1_2\n"}elsif($a[-2] eq Lymphoid3){print "$a[0],$b","$a[-3],$a[-2]\n"}elsif($a[-2] eq Myeloid7 || $a[-2] eq Myeloid8){print "$a[0],$b","$a[-3],$a[-2]\n"}}' pancancer.icar.all.cell.obs > pancancer.icar.all.cell.obs.at

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];<IN1>;print "cell,merge_celltype\n";while(<IN1>){chomp;@a=split(/,/);if($ha{$a[-6]} eq "delete"){next};if(exists $ha{$a[-6]}){print "$a[0],$ha{$a[-6]}\n"}}' list pancancer.icar.all.cell.obs > pancancer.icar.all.cell.obs.at.at

