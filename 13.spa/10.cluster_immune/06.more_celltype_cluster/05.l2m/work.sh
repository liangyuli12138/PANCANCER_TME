perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;s/SS200000929BL/SS200000929BL_D2/;@a=split(/,/);$ha{$a[1]}="$a[-3],$a[-2],$a[-1]"};open IN2,$ARGV[2];while(<IN2>){chomp;s/SS200000929BL/SS200000929BL_D2/;@a=split(/,/);$hb{$a[0]}="$a[-2]"};open IN1,$ARGV[1];$t=<IN1>;chomp $t;$t=~s/\t/,/g;print "$t,density,area,elongation,distance,Dominant\n";while(<IN1>){chomp;@a=split;s/\t/,/g;$t=$a[1]=~/Lymphoid/?"Lymphoid_Dominant_ICAR":"Myeloid_Dominant_ICAR";print "$_,$ha{$a[0]},$hb{$a[0]},$t\n"}' cell_group_meta.csv l2m.stat.csv cell_group_meta.border.csv > l2m.stat.merge.csv

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$_};open IN1,$ARGV[1];<IN1>;print "cell,groups,Lymphoid_percentage,Myeloid_percentage,Lym2Mye_ratio,log10_Lym2Mye_ratio,density,area,elongation,distance,Dominant\n";while(<IN1>){chomp;@a=split(/,/);print "$ha{$a[0]}\n"}' l2m.stat.merge.csv immune.cluster.11.r1.5.new.obs > l2m.stat.merge.csv.at

perl -e 'while(<>){chomp;print "sc.pl.umap(adata_concat,color=\"$_\",frameon=False, na_color=\"grey\",save=\".$_.gene.png\")\n"}' list >> umap.py 

perl -e 'while(<>){chomp;$a=$_;$o=`cat tmp.py`;$o=~s/aaaa/$a/g;print "$o"}' list >> plot.boxplot.py 

