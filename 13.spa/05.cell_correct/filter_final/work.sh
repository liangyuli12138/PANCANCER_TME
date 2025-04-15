perl -e '$o=`cat head.py`;print "$o";while(<>){chomp;@a=split;$s=$a[0];$o=`cat tmp.py`;$o=~s/aaaa/$s/g;print "$o"}' sample_list > all.filter.py

perl -e 'while(<>){chomp;@a=split;$s=$a[0];$h=`cat head.py`;$o=`cat tmp.py`;$o=~s/aaaa/$s/g;open OUT,">shell/$a[0].py";print OUT "$h\n$o";open OUT,">shell/$a[0].sh";print OUT "/jdfssz1/ST_TSCBI/P22Z10200N0433/USER/zhangzhao/software/anaconda3/bin/python $a[0].py\n"}' sample_list 

les sample_list |while read i;do echo -e "/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R422/bin/R --slave -f /jdfssz1/ST_TSCBI/P22Z10200N0433/USER/liuhui/script/h5ad_cov_rds.R '--args /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/$i/$i""_cellbin.final.h5ad /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/05.cell_correct/result/$i/$i""_cellbin.final.rds'";done > h5ad2rds.sh


