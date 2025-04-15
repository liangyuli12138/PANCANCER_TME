while(<>){
chomp;s/\.at//;$a=$_;
open IN,"/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/out/$a.obs.csv";
$n=-1;$m=0;$nn=0;$o="";
while(<IN>){
chomp;$n++;
@b=split(/,/);
if($b[0] eq "Lymphoid_B"){$b[$n]=$b[$n]==""?0:$b[$n];$o.= "$a\t$b[$n]\t";$x=$b[$n];$m=1};
}
open IN,"/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/out/$a.shuffled.csv";
$n=-1;
while(<IN>){
chomp;$n++;
@b=split(/,/);
if($b[0] eq "Lymphoid_B"){$b[$n]=$b[$n]==""?0:$b[$n];$o.= "$b[$n]\t";$y=$b[$n];};
}
open IN,"/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/09.cluster_immune/07.toself/at/$a.at";
while(<IN>){if(/Lymphoid_B/){$nn++}}
$o.= "$nn\t";
if($m>0){
$k=$x/$nn-$y/$nn;print "$o\t$k\n"}
if($m==0){print "$a\t0\t0\t0\t0\n"}
}
