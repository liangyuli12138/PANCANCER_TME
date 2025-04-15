perl -e '$o=`cat head.py`;print "$o";while(<>){chomp;@a=split;$s=$a[0];$o=`cat tmp.py`;$o=~s/aaaa/$s/g;print "$o"}' sn.list > all.filter.py

