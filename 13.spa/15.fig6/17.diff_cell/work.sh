perl -e 'while(<>){chomp;@a=split(/,/);print "$a[-5]\t$a[-3]\n"}' pancancer.icar.all.cell.filter.obs |sort -u > filter.all.cluster.list

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN2,$ARGV[2];while(<IN2>){chomp;@a=split;$hb{$a[1]}=1};open IN1,$ARGV[1];$t=<IN1>;chomp $t;print "cell,difftype\n";while(<IN1>){chomp;@a=split(/,/);if(($a[-5]=~/Lym/) && exists $hb{$a[-3]}){print "$a[0],$ha{$a[-7]}","_","$a[-5]\n"}else{print "$a[0],other\n"}}' merge.list pancancer.icar.all.cell.filter.obs filter.all.cluster.list > filter.all.cluster.list.at

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split;$ha{$a[0]}=$a[1]};open IN2,$ARGV[2];while(<IN2>){chomp;@a=split;$hb{$a[1]}=1};open IN1,$ARGV[1];$t=<IN1>;chomp $t;print "cell,difftype\n";while(<IN1>){chomp;@a=split(/,/);if(($a[-5]=~/Lym/) && exists $hb{$a[-3]} && exists $ha{$a[-7]}){print "$a[0],$ha{$a[-7]}","_","$a[-5]\n"}else{print "$a[0],other\n"}}' merge.list2 pancancer.icar.all.cell.filter.obs filter.all.cluster.list > filter.all.cluster.list.at2

