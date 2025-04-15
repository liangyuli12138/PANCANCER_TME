perl -e 'while(<>){chomp;@a=split;$o=`cat network.py`;$o=~s/aaaa/$a[0]/g;open OUT,">shell/$a[0].network.py";print OUT "$o";open OUT,">shell/$a[0].network.sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $a[0].network.py\n"}' sn.list 

