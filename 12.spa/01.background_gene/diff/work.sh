perl -e 'while(<>){chomp;@a=split(/\//);$a[-1]=~s/diff\.//;$a[-1]=~s/\.csv//;$n=0;open IN,$_;while(<IN>){chomp;@b=split(/\,/);$n++;if($n>100){last};$b[-2]=sprintf("%.2f",$b[-2]);$b[-1]=sprintf("%.2f",$b[-1]);$ha{$b[1]}.="{$a[-1]|$b[-2]|$b[-1]};"}};for $i(keys %ha){print "$i\t$ha{$i}\n"}' all.diff.list |les

les all.stat.list |while read i;do echo perl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/12.spa/01.background_gene/diff/add.diff.pl /zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/12.spa/01.background_gene/diff/all.diff.list $i \> $i.diff;done > all.sh

