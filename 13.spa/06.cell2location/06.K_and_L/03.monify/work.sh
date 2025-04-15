cut -f 1 sn.list |while read i;do echo perl get.merge.pl $i""_cellbin.final.leiden.obs.csv $i""_cellbin.final.celltype.obs ">" $i""_cellbin.final.celltype.group.obs;done|sh

cut -f 1 sn.list |while read i;do echo perl stat.pl $i""_cellbin.final.celltype.group.obs \> $i""_cellbin.final.celltype.group.obs.stat.xls;done|sh

perl -e 'while(<>){chomp;@a=split;@b=split(/,/,$a[1]);$t="$a[0]"."_cellbin.final.celltype.group.obs";$p="$a[0]"."_cellbin.final.celltype.group.obs.monify";open OUT,">$p";open IN,$t;$f=<IN>;print OUT "$f";while(<IN>){chomp;@c=split(/,/);for($i=0;$i<@b;$i++){if($b[$i] == $c[-1] && $c[-4] eq "Epithelium_Normal"){s/Epithelium_Normal/Hepatocyte/;s/Epithelium/Hepatocyte/}};print OUT "$_\n"}}' monify.list.liver


perl monify.kidney.pl monify.list.kidney
