open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);if($a[3]=~/Cluster/){$ha{$a[0]}=$a[3]}};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);if($a[3]=~/Cluster/){if($ha{$a[0]} ne $a[3]){print "$_\t$ha{$a[0]}\n"}}}

