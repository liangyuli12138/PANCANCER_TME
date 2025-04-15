perl -e 'open IN0,$ARGV[0];while(<IN0>){chomp;s/\"//g;@a=split(/,/);if($a[2]=~/group/){$ha{$a[1]}=$a[2]}};open IN1,$ARGV[1];while(<IN1>){chomp;s/\"//g;@a=split(/,/);if(exists $ha{$a[1]}){print "$_,$ha{$a[1]}\n"}}' D01972C4.group.at D01972C4.ori.at > D01972C4.group.at.cell

perl -e 'while(<>){chomp;@a=split(/,/);$ha{$a[2]}{$a[4]}++;$hb{$a[4]}=1};for $i(sort {$a cmp $b} keys %hb){print ",$i"};print "\n";for $j(sort {$a cmp $b} keys %ha){print "$j";for $i(sort {$a cmp $b} keys %hb){$ha{$j}{$i}=$ha{$j}{$i}?$ha{$j}{$i}:0;print ",$ha{$j}{$i}"};print "\n"}' D01972C4.group.at.cell > D01972C4.group.at.cell.stat

grep Lymphoid_B_memory D01972C4.group.at.cell |cut -d "," -f 2 > Lymphoid_B_memory.index

