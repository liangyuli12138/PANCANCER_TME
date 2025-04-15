perl -e 'while(<>){chomp;@a=split;$t="$a[1].$a[0]";$o=`cat spatial_neighbors.merge.py`;$o=~s/aaaa/$a[0]/g;$o=~s/bbbb/$t/g;open OUT,">shell1/$t.py";print OUT "$o";open OUT,">shell1/$t.sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $t.py\n"}' sn.list

perl -e 'while(<>){chomp;@a=split;$t="$a[1].$a[0]";$o=`cat spatial_neighbors.sub.py`;$o=~s/aaaa/$a[0]/g;$o=~s/bbbb/$t/g;open OUT,">shell2/$t.py";print OUT "$o";open OUT,">shell2/$t.sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $t.py\n"}' sn.list

perl -e 'while(<>){chomp;@a=split;$t="$a[1].$a[0]";$o=`cat spatial_neighbors.merge.global.py`;$o=~s/aaaa/$a[0]/g;$o=~s/bbbb/$t/g;open OUT,">shell/merge.global/$t.py";print OUT "$o";open OUT,">shell/merge.global/$t.sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $t.py\n"}' sn.list
perl -e 'while(<>){chomp;@a=split;$t="$a[1].$a[0]";$o=`cat spatial_neighbors.merge.local.py`;$o=~s/aaaa/$a[0]/g;$o=~s/bbbb/$t/g;open OUT,">shell/merge.local/$t.py";print OUT "$o";open OUT,">shell/merge.local/$t.sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $t.py\n"}' sn.list

perl stat.pl csv.list2 > nbor_pvals.out.merge.global.stat.xls

perl -e 'while(<>){chomp;@a=split;$t="$a[1].$a[0]";$o=`cat spatial_neighbors.sub.abs.py`;$o=~s/aaaa/$a[0]/g;$o=~s/bbbb/$t/g;open OUT,">shell/sub.abs/$t.py";print OUT "$o";open OUT,">shell/sub.abs/$t.sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $t.py\n"}' sn.list
perl -e 'while(<>){chomp;@a=split;$t="$a[1].$a[0]";$o=`cat spatial_neighbors.merge.abs.py`;$o=~s/aaaa/$a[0]/g;$o=~s/bbbb/$t/g;open OUT,">shell/merge.abs/$t.py";print OUT "$o";open OUT,">shell/merge.abs/$t.sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $t.py\n"}' sn.list

