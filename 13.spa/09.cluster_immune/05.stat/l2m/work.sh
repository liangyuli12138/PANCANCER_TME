perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;s/SS200000929BL/SS200000929BL_D2/;@a=split(/,/);$ha{$a[1]}="$a[-3],$a[-2],$a[-1]"};open IN1,$ARGV[1];$t=<IN1>;chomp $t;$t=~s/\t/,/g;print "$t,density,area,elongation\n";while(<IN1>){chomp;@a=split;s/\t/,/g;print "$_,$ha{$a[0]}\n"}' cell_group_meta.csv l2m.stat.csv > l2m.stat.merge.csv

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$_};open IN1,$ARGV[1];<IN1>;print "cell,group,Lymphoid_percentage,Myeloid_percentage,Lym/Mye_ratio,log10(Lym\/Mye_ratio),density,area,elongation\n";while(<IN1>){chomp;@a=split(/,/);print "$ha{$a[0]}\n"}' l2m.stat.merge.csv immune.cluster.r0.5.obs > l2m.stat.merge.csv.at

perl -e 'while(<>){chomp;print "sc.pl.umap(adata_concat,color=\"$_\",frameon=False, na_color=\"grey\",save=\".$_.gene.png\")\n"}' list >> umap.py 

