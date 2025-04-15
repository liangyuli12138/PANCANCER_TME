perl -e '<>;while(<>){chomp;@a=split(/,/);$ha{$a[-3]}{$a[-1]}++;$hb{$a[-3]}=$a[-2];$hc{$a[-1]}=1;$hd{$a[-3]}++};print "groups\tcluster";for $i(sort {$a cmp $b} keys %hc){print "\t$i"};print "\n";for $j(sort {$a cmp $b} keys %hb){print "$hb{$j}\t$j";for $i(sort {$a cmp $b} keys %hc){$ha{$j}{$i}=$ha{$j}{$i}?$ha{$j}{$i}:0;$o=$ha{$j}{$i}/$hd{$j};print "\t$o"};print "\n"}' pancancer.icar.all.cell.obs > pancancer.icar.all.cell.obs.stat.xls

perl -e '<>;while(<>){chomp;@a=split(/,/);$ha{$a[-3]}{$a[-1]}++;$hb{$a[-3]}=$a[-2];$hc{$a[-1]}=1;$hd{$a[-3]}++};print "groups\tcluster";for $i(sort {$a cmp $b} keys %hc){print "\t$i"};print "\n";for $j(sort {$a cmp $b} keys %hb){print "$hb{$j}\t$j";for $i(sort {$a cmp $b} keys %hc){$ha{$j}{$i}=$ha{$j}{$i}?$ha{$j}{$i}:0;$o=$ha{$j}{$i}/$hd{$j};print "\t$ha{$j}{$i}"};print "\n"}' pancancer.icar.all.cell.obs > pancancer.icar.all.cell.obs.num.xls

perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;$ha{$_}=1};open IN1,$ARGV[1];$t=<IN1>;print "$t";while(<IN1>){chomp;@a=split;if(exists $ha{$a[1]}){print "$_\n"}}'  filter.list pancancer.icar.all.cell.obs.stat.xls > pancancer.icar.all.cell.obs.stat.filter.xls

