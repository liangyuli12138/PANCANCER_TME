KEGG work in software-0-1 
TMPDIR=/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/13.spa/15.fig6/14.cell_diff/all2/diff_merge/fun /ldfssz1/ST_OCEAN/USER/liaoshangfeng/software/anaconda3/envs/R411/bin/Rscript enrich.kegg.R


perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;s/\"//g;@a=split(/,/);$ha{$a[1]}=$a[0]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/ - /);if(@a>2){$a[0]="$a[0] - $a[1]"};$a[0]=~s/^\d+\.\s+//;s/^\d+\.\s+//;$hb{$a[0]}=$_;};open IN2,$ARGV[2];$t=<IN2>;chomp $t;print "$t\tChinese\n";while(<IN2>){chomp;@a=split(/\t/);@b=split(/\//,$a[-2]);$o="";for($i=0;$i<@b;$i++){$o.="$ha{$b[$i]}/"};$a[-2]=$o;print join "\t",@a;print "\t$hb{$a[2]}\n"}'  ENTREZID_SYMBOL.csv kegg.term.list.csv.kimi.csv result.enrich.KEGG.csv > result.enrich.KEGG.kimi.csv

mv result.enrich.GO.csv result.enrich.GO.txt

