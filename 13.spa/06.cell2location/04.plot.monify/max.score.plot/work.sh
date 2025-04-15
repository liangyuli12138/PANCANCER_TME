perl -e '$o=`cat head.py`;print "$o";while(<>){chomp;@a=split;$o=`cat tmp.py`;$c="$a[1].$a[0]";$o=~s/aaaa/$c/g;print "$o"}' sn.list > plot.py

/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python plot.py &

