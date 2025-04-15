cat /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/03.TLS/region/*csv > all.immue.csv

perl get.at.pl merge.list all.immue.csv obs.list all.tls.windows.list

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split;$o=`cat tem.R`;open OUT,">shell/$a[0].R";$t=$a[0];$t=~s/\_TLS\S+//;$d=$1 if($a[0]=~/(TLS\S+)/);$o=~s/aaaa/$t/g;$o=~s/bbbb/$t\.$ha{$t}/g;$o=~s/cccc/$d/g;$o=~s/xxaa/$a[5]/g;$o=~s/xxbb/$a[6]/g;$o=~s/yyaa/$a[7]/g;$o=~s/yybb/$a[8]/g;print OUT "$o";open OUT,">shell/$a[0].sh";print OUT "/ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R411/bin/Rscript $a[0].R\n"}' sn.list all.tls.windows.list

