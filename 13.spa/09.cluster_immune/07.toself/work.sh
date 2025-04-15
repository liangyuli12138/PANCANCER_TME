perl -e 'while(<>){chomp;s/\.at//;$t=$_;$t=~s/\_Cl\S+//;$o=`cat spatial_neighbors.merge.global.py`;$o=~s/aaaa/$_/g;$o=~s/bbbb/$t/g;open OUT,">shell/$_.py";print OUT "$o";open OUT,">shell/$_.sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $_.py\n"}' cluster.list

perl get.stat.pl cluster.list > all.toself.stat.csv

perl -e 'while(<>){chomp;s/\.at//;$t=$_;$t=~s/\_Cl\S+//;$o=`cat spatial_neighbors.merge.abs.py`;$o=~s/bbbb/$_/g;$o=~s/aaaa/$t/g;open OUT,">shell/$_.py";print OUT "$o";open OUT,">shell/$_.sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $_.py\n"}' cluster.list

