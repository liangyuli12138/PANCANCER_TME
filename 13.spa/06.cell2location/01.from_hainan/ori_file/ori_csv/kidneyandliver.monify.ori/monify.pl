while(<>){
if(/^#/){next};
chomp;@a=split;
$p="Merge_"."$a[1].$a[0]"."_meta_ori.monify.csv";
$q="Merge_"."$a[1].$a[0]"."_meta_ori.monify.csv.list";
`mv ../$p $p`;`mv ../$q $q`;
$t="/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/06.cell2location/06.K_and_L/03.monify/$a[0]"."_cellbin.final.celltype.group.obs.monify";
open IN,$t;<IN>;
undef %ha;
while(<IN>){chomp;@c=split(/,/);$ha{$c[0]}=$c[-4]};

open IN,$p;open OUT,">../$p";
$x=<IN>;print OUT $x;
while(<IN>){chomp;@d=split(/,/);$id="$d[3]"."_"."$d[4]";if(!exists $ha{$id}){print OUT "$_\n"}else{$d[-2]="\"$ha{$id}\"";print OUT join (",",@d);print OUT "\n"}};

open IN,$q;open OUT,">../$q";
$x=<IN>;print OUT $x;
while(<IN>){chomp;@d=split(/,/);$d[1]=~s/\"//g;if(!exists $ha{$d[1]}){print OUT "$_\n"}else{$d[-3]="\"$ha{$d[1]}\"";print OUT join (",",@d);print OUT "\n"}}
}
