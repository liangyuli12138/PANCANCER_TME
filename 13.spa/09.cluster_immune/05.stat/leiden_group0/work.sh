perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);if($a[1] eq "0"){$ha{$a[0]}=1}};open IN1,$ARGV[1];$t=<IN1>;print "$t";while(<IN1>){chomp;@a=split(/,/);if(exists $ha{$a[0]}){print "$_\n"}}' immune.cluster.r0.5.obs cluster.immune.cell.csv > cluster.immune.cell.group0.csv

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);if($a[1] eq "0"){$ha{$a[0]}=1}};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);$hb{$a[0]}=$a[1]};open IN2,$ARGV[2];$t=<IN2>;print "group_pri,group_sub,$t";while(<IN2>){chomp;@a=split(/,/);$x=$hb{$a[0]}?$hb{$a[0]}:"NA";print "$ha{$a[0]},$x,$_\n"}'

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;@a=split(/,/);$ha{$a[0]}=$a[1]};open IN1,$ARGV[1];while(<IN1>){chomp;@a=split(/,/);$hb{$a[0]}=$a[1]};open IN2,$ARGV[2];$t=<IN2>;$t=~s/^,/id,/;print "group_pri,group_sub,$t";while(<IN2>){chomp;@a=split(/,/);if(!exists $hb{$a[0]}){$hb{$a[0]}="NA"};print "$ha{$a[0]},$hb{$a[0]},$_\n"}' immune.cluster.r0.5.obs immune.group0.r0.5.obs cluster.immune.cell.csv > cluster.immune.cell.csv.stat.xls

