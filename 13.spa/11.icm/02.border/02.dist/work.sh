perl -e 'while(<>){chomp;@a=split;$o=`cat tmp.py`;$o=~s/aaaa/$a[0]/g;open OUT,">shell/$a[0].py";print OUT "$o";open OUT,">shell/$a[0].sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $a[0].py\n"}' all.sn.list 

perl -e 'while(<>){chomp;@a=split;print "perl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/11.icm/02.border/02.dist/get.dist.pl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/11.icm/01.plot/at/$a[0].all.at > /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/11.icm/02.border/02.dist/out/$a[0].dist.csv\n"}' all.sn.list > all.dist.sh

