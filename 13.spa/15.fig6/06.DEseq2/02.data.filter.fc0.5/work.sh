perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;s/\,+$//;@a=split(/,/);for($i=1;$i<@a;$i++){$ha{$a[$i]}=$a[0]}};open IN1,$ARGV[1];$t=<IN1>;chomp $t;print "$t,GOPB\n";while(<IN1>){chomp;s/\"//g;@a=split(/,/);$ha{$a[0]}=$ha{$a[0]}?$ha{$a[0]}:"NA";print "$_,$ha{$a[0]}\n"}' gopb_gene.csv Lymphoid3up_vs_Lymphoid_1_2down_FC0.5_gene.res1_total.csv > Lymphoid3up_vs_Lymphoid_1_2down_FC0.5_gene.res1_total.csv.xls

perl -e 'print "Lymphoid3\tLymphoid_1_2\n";while(<>){chomp;s/\"//g;@a=split(/,/);if($a[1]>5){if($a[2]>1){$x.="$a[0]\t"}elsif($a[2]<-1){$y.="$a[0]\t"}}};print "$x\n$y"' Lymphoid3up_vs_Lymphoid_1_2down_FC0.5_gene.res1_total.csv > for.annotation.lym123.xls

