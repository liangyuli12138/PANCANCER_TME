perl -e 'while(<>){chomp;@a=split;$o=`cat network.100.py`;$o=~s/aaaa/$a[0]/g;open OUT,">shell/$a[0].100.network.py";print OUT "$o";open OUT,">shell/$a[0].100.network.sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $a[0].100.network.py\n"}' sn.list

perl -e 'while(<>){chomp;@a=split;$o=`cat network.200.py`;$o=~s/aaaa/$a[0]/g;open OUT,">shell/$a[0].200.network.py";print OUT "$o";open OUT,">shell/$a[0].200.network.sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $a[0].200.network.py\n"}' sn.list

perl -e 'while(<>){chomp;@a=split;$o=`cat network.300.py`;$o=~s/aaaa/$a[0]/g;open OUT,">shell/$a[0].300.network.py";print OUT "$o";open OUT,">shell/$a[0].300.network.sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $a[0].300.network.py\n"}' sn.list

