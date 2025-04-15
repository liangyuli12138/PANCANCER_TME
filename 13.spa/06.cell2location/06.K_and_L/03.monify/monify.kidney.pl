while(<>){
if(/^#/){next};
chomp;@a=split;@b=split(/,/,$a[1]);
$t="$a[0]"."_cellbin.final.celltype.group.obs";$p="$a[0]"."_cellbin.final.celltype.group.obs.monify";open OUT,">$p";
$q="Merge_Kidney.$a[0]"."_meta_ori.csv";
open IN,$q;
<IN>;
while(<IN>){chomp;s/\"POLYGON.+\)\)\"\,//;s/\"//g;@d=split(/,/);$id="$d[3]"."_"."$d[4]";$hd{$id}=$d[-1]}
open IN,$t;$f=<IN>;print OUT "$f";while(<IN>){chomp;@c=split(/,/);for($i=0;$i<@b;$i++){if($b[$i] == $c[-1] && $c[-4] eq "Epithelium_Normal"){s/Epithelium_Normal/$hd{$c[0]}/;}}print OUT "$_\n"}

}
