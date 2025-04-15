perl -e 'while(<>){chomp;@a=split;$o=`cat tmp.py`;$o=~s/aaaa/$a[0]/g; print "$o"}' sn.list >> merge.py 

