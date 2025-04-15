cat ../result/* > merge.doublet.csv

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$a[-1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);$ha{$a[0]}=$ha{$a[0]}?$ha{$a[0]}:"Missing";$hb{$a[-1]}{$ha{$a[0]}}++};for $i(sort {$a <=> $b} keys %hb){$x="Doublet";$y="Singlet";$z="Missing";$hb{$i}{$x}=$hb{$i}{$x}?$hb{$i}{$x}:0;$hb{$i}{$y}=$hb{$i}{$y}?$hb{$i}{$y}:0;$hb{$i}{$z}=$hb{$i}{$z}?$hb{$i}{$z}:0;print "$i\t$hb{$i}{$x}\t$hb{$i}{$y}\t$hb{$i}{$z}\n"}' merge.doublet.csv pancancer.ref.0721.raw.obs.csv > pancancer.ref.0721.raw.obs.csv.check

les pancancer.ref.0721.raw.obs.csv.check|perl -e 'while(<>){chomp;@a=split(/\t/);$a+=$a[1];$b+=$a[2];$c+=$a[3]};print "$a\t$b\t$c\n"'
85593	1055616	1

