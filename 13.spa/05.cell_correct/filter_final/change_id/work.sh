rename final tmp ../../result/*/*

perl -e '$o=`cat head.py`;print "$o";while(<>){chomp;@a=split;$s=$a[0];$o=`cat tmp.py`;$o=~s/aaaa/$s/g;print "$o"}' sample_list > all.obs.py

for i in ../../result/*/*cellbin.tmp.obs;do echo perl center2id.pl $i \> $i.id;done|sh &

perl -e '$o=`cat head.py`;print "$o";while(<>){chomp;@a=split;$s=$a[0];$o=`cat tmp2.py`;$o=~s/aaaa/$s/g;print "$o"}' sample_list > all.final.py

