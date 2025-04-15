perl -e 'while(<>){chomp;@a=split;$o=`cat get.mean.py`;$o=~s/aaaa/$a[0]/g;open OUT,">shell/$a[0].py";print OUT "$o";open OUT,">shell/$a[0].sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $a[0].py\n"}' sn.list

perl -e '$o=`cat header`;print "$o";while(<>){chomp;@a=split(/\t/);$t="out/$a[0].geneAmean.csv";open IN,$t;<IN>;while(<IN>){chomp;print "$a[0]","$_\n"}}' sn.list > merge.gene.csv


