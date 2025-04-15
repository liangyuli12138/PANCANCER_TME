perl -e '<>;while(<>){chomp;@a=split(/,/);$ha{$a[-1]}++;$hb{$a[-1]}{$a[-2]}++};for $i(sort {$a<=>$b} keys %hb){print "$i";for $j(sort {$hb{$i}{$b} <=> $hb{$i}{$a}} keys %{$hb{$i}}){$p=sprintf "%.2f",$hb{$i}{$j}/$ha{$i};print "\t$j|$hb{$i}{$j}|$p"};print "\n"}' pancancer.ref.umap.0723.obs.txt > pancancer.ref.umap.0723.obs.txt.stat &

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);print "$a[0],$ha{$a[-1]}\n"}' pancancer.ref.umap.0723.obs.txt.stat.rela pancancer.ref.umap.0723.obs.txt > pancancer.ref.umap.0723.obs.txt.stat.rela.at

