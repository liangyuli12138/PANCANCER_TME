perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$hb{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);if(/Lymphoid/){s/Cluster\S+//;$ha{$_}++;$hc{$_}+=$hb{$a[0]}}};for $i(keys %ha){print "$i\t$ha{$i}\t$hc{$i}\n"}' all.icar.cell.csv immune.cluster.11.r1.5.new.obs > tls.stat

