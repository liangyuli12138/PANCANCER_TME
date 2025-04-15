open IN0,$ARGV[0];
while(<IN0>){
chomp;
$t=$_;
$t=~s/region\///;
$t=~s/_immue_region.csv//;
$hb{$t}=1;
open IN,$_;<IN>;
while(<IN>){
chomp;@a=split(/,/);
$ha{$t}{$a[0]}=$a[1];
}
}

open IN1,$ARGV[1];
while(<IN1>){chomp;
@a=split;
$x="/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/08.HE/01.Hexcell/region/$a[0]_region.csv";
#$y="/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/01.from_hainan/ori_file/h5ad/data/$a[0]_cellbin.final.celltype.obs.csv";
#undef %hb;undef %hc;
if(!exists $hb{$a[0]}){next};
open OUT0,">at/$a[0].input";
open OUT1,">at/$a[0].at";
print OUT0 "cell\n";
print OUT1 "cell,region\n";

open INX,$x;
<INX>;
while(<INX>){chomp;@b=split(/,/);
$ha{$a[0]}{$b[0]}=$ha{$a[0]}{$b[0]}?$ha{$a[0]}{$b[0]}:"others";
print OUT0 "$b[0]\n";print OUT1 "$b[0],$ha{$a[0]}{$b[0]}\n";
}
}
