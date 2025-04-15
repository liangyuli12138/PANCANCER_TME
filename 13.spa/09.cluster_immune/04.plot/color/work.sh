perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;s/group//;@a=split(/,/);$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];<IN1>;while(<IN1>){chomp;s/Cluster/,/;@a=split(/,/);$hb{$a[0]}{$a[1]}=$ha{$a[2]}};for $i(keys %hb){open OUT,">data/$i.color";print OUT "Sub_cluster,Colour\n";for $j(sort {$a <=> $b} keys %{$hb{$i}}){print OUT "$j,$hb{$i}{$j}\n"}}' colour.group.list immune.cluster.r0.5.obs 

