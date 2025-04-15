perl -e '$o=`cat head.R`;print "$o";while(<>){chomp;@a=split;$s=$a[0];$o=`cat tmp.R`;$o=~s/outfix_sample/$s/g;print "$o"}' sn.list > all.filter.R

perl -e 'print "cellbin\n";<>;while(<>){s/"//g;@a=split(/,/);print "$a[0]\n"}' SS200000929BL_D2_cellbin.filter.gene.list > cellbin.filter.list

ls ../result/*/*filter.gene.list|while read i;do echo perl getcell.pl $i \> $i.cellbin;done|sh &

