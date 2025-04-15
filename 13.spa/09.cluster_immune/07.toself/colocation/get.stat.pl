while(<>){
chomp;s/\.at//;$a=$_;
open IN,"/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/out/$a.out.global.zscore.csv";
$n=-1;$m=0;
while(<IN>){
chomp;$n++;
@b=split(/,/);
if($b[0] eq "Lymphoid_B"){$b[$n]=$b[$n]==""?0:$b[$n];print "$a\t$b[$n]\t";$m=1};
}
open IN,"/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/out/$a.out.global.nbor_pvals.csv";
$n=-1;
while(<IN>){
chomp;$n++;
@b=split(/,/);
if($b[0] eq "Lymphoid_B"){$b[$n]=$b[$n]==""?1:$b[$n];print "$b[$n]\n"};
}

if($m==0){print "$a\t0\t1\n"}
}
