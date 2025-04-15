open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$n++;$ha{$a[0]}=$n};
open IN1,$ARGV[1];while(<IN1>){chomp;@a=split;s/^[^\t]+\t//;s/\t/\",\"/g;s/^/\"/;s/$/\"/;$hb{$a[0]}=$_};
open IN2,$ARGV[2];while(<IN2>){chomp;@a=split(/,/);$hc{$a[-1]}{$a[0]}=1};
$h=`cat head.py`;
print "$h";

for $i(keys %hc){
print "adata_concat = sc.read_h5ad(\"/zfssz2/ST_TSCBI/P22Z10200N0433/USER/wubin2/wubin2/pancnew/09.cluster_show/02.split/pancancer.split.$i.0905.h5ad\")\nmarker_genes_dict = {\n";
$s="";
for $j(sort {$ha{$a} <=> $ha{$b}} %ha){
if(exists $hc{$i}{$j}){
print "\"$j\":[$hb{$j}],\n";
$s.="\"$j\",";
}}
print "}\n";
print "sc.pl.dotplot(adata_concat, marker_genes_dict, groupby=\"groups_pri\",dot_max=0.8, categories_order=[$s],use_raw=False, colorbar_title=\"mean z-score\", vmin=-1, vmax=1, cmap=\"RdBu_r\",save=\".sub.cluster.markergene.$i.pdf\")\n\n\n";

}
