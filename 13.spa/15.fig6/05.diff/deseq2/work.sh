perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;s/\,+$//;@a=split(/,/);for($i=1;$i<@a;$i++){$ha{$a[$i]}=$a[0]}};open IN1,$ARGV[1];$t=<IN1>;chomp $t;print "$t,gopb\n";while(<IN1>){chomp;s/\"//g;@a=split(/,/);$ha{$a[0]}=$ha{$a[0]}?$ha{$a[0]}:"NA";print "$_,$ha{$a[0]}\n"}' ../diff/gopb_gene.csv DESeq2_lym3_vs_lym12.csv > DESeq2_lym3_vs_lym12.csv.gopb

