perl -e 'while(<>){chomp;@a=split;print "ln -s /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/$a[0]/$a[0]","_cellbin.final.h5ad cluster\n"}' sn.list |sh

perl -e 'while(<>){chomp;@a=split;$o=`cat tmp.py`;$o=~s/aaaa/$a[0]/g;open OUT,">cluster/$a[0].py";print OUT "$o";open OUT,">cluster/$a[0].sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $a[0].py\n"}'  sn.list

